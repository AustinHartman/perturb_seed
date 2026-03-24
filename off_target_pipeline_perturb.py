#!/usr/bin/env python3
"""
Script for the nomination of off-target CRISPRi dCas9 events using perturb-seq data

this pipeline identifies off-target candidates by:
1. Clustering pseudobulked transcriptomes by guide and identifying guide neighborhoods
2. Identifying guides which knock down one of the genes targeted by a different guide in the neighborhood
3. Searching the guide sequence for seed matches in the promoter of the off-target gene

Input to this script is provided in a yaml file. See the README for details on input data and parameters.

python off_target_pipeline_perturb.py --params <dataset>/params.yaml
"""

import sys
import argparse
import time
from pathlib import Path
import gc
import re
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from scipy.sparse import issparse, csr_matrix as csr


def load_params(yaml_path: Path) -> dict:
    """Load pipeline parameters from a YAML file."""
    import yaml
    with open(yaml_path) as fh:
        p = yaml.safe_load(fh)

    base = yaml_path.parent

    def resolve(key):
        if p.get(key):
            v = Path(p[key])
            return v if v.is_absolute() else base / v
        return None

    p["adata_path"] = resolve("adata_path")
    p["sgrna_csv"] = resolve("sgrna_csv")
    p["output_csv"] = resolve("output_csv")
    p["output_dir"] = resolve("output_dir")
    p["tss_map_csv"] = resolve("tss_map_csv")
    p["genome_fa"] = resolve("genome_fa")
    p["genome_fai"] = resolve("genome_fai")

    p.setdefault("obs_guide_col", "guide_id")
    p.setdefault("obs_gene_col", "perturbed_gene_name")
    p.setdefault("nt_label", "NTC")
    p.setdefault("guide_filter_col", None)
    p.setdefault("guide_filter_val", None)
    p.setdefault("dataset_type", "single_guide")
    p.setdefault("min_kd_pairs", 1)
    p.setdefault("seed_length", 10)
    p.setdefault("n_top_features", 1_000)
    p.setdefault("n_neighbors", 30)
    p.setdefault("kd_lfc_thresh", 1.0)
    p.setdefault("log2fc_pseudocount", 0.1)
    p.setdefault("min_cells_per_guide_vector", 20)
    p.setdefault("promoter_upstream", 500)
    p.setdefault("promoter_downstream", 500)
    p.setdefault("n_seed_workers", 30)

    return p


P: dict = {}

CHUNK_SIZE = 250_000 # num cells per chunk when streaming h5ad into memory
NPCS = 100 # num PCs for guide neighborhood detection
RC_MAP = str.maketrans("ACGTacgt", "TGCAtgca") # dict for rev comping DNA sequences
_TSS_MAP: dict | None = None


def compute_pseudobulk():
    """
    Load h5ad in backed mode, normalize cells to 10k, log1p, and compute
    mean expression for every gene group and every guide group.

    If there is just one guide vector per gene, then the guide and gene pseudobulks will be the same.

    Returns
    -------
    pb_gene      : DataFrame  (n_target_genes, n_features)
    pb_guide     : DataFrame  (n_guides_passing_min_cells, n_features)
    var          : DataFrame  h5ad .var metadata (with 'std' column)
    obs          : DataFrame  filtered h5ad .obs (cells from passing guides only)
    excluded_df  : DataFrame  (guide_id, n_cells) for guides below min_cells
    """
    adata_path = P["adata_path"]
    obs_guide_col = P["obs_guide_col"]
    obs_gene_col = P["obs_gene_col"]
    guide_filter_col = P.get("guide_filter_col")
    guide_filter_val = P.get("guide_filter_val")
    nt_label = P["nt_label"]

    # could simplify this by loading whole thing into mem
    # but many GWPS datasets would require 256GB+ RAM
    adata_backed = sc.read_h5ad(adata_path, backed="r")
    n_total = adata_backed.n_obs # num cells
    n_genes = adata_backed.n_vars # num features
    var = adata_backed.var.copy()

    if guide_filter_col is None:
        sg_mask = np.ones(n_total, dtype=bool)  # use all cells
    else:
        sg_mask = (adata_backed.obs[guide_filter_col] == guide_filter_val).values
    obs = adata_backed.obs[sg_mask].copy()
    n_cells = int(sg_mask.sum())
    abs_idx = np.where(sg_mask)[0]
    sort_perm = np.argsort(abs_idx)
    abs_idx_sorted = abs_idx[sort_perm]

    gene_labels_sorted = obs[obs_gene_col].astype(str).values[sort_perm]
    guide_labels_sorted = obs[obs_guide_col].astype(str).values[sort_perm]

    gene_list = sorted(set(gene_labels_sorted))
    guide_list = sorted(set(guide_labels_sorted))
    gene_map = {g: i for i, g in enumerate(gene_list)}
    guide_map = {g: i for i, g in enumerate(guide_list)}
    n_gg = len(gene_list)
    n_pg = len(guide_list)

    gene_codes = np.array([gene_map[g] for g in gene_labels_sorted], dtype=np.int32)
    guide_codes = np.array([guide_map[g] for g in guide_labels_sorted], dtype=np.int32)
    nt_code = gene_map.get(nt_label, -1)

    gene_sum = np.zeros((n_gg, n_genes), dtype=np.float64)
    guide_sum = np.zeros((n_pg, n_genes), dtype=np.float64)
    gene_cnt = np.zeros(n_gg, dtype=np.int64)
    guide_cnt = np.zeros(n_pg, dtype=np.int64)
    nt_n, nt_s, nt_s2 = 0, np.zeros(n_genes), np.zeros(n_genes)

    n_chunks = (n_cells + CHUNK_SIZE - 1) // CHUNK_SIZE
    print(f"  Streaming {n_chunks} chunks of ≤{CHUNK_SIZE:,} cells…", flush=True)

    for ci in range(n_chunks):
        # figure out which cells to load for the chunk
        lo, hi = ci * CHUNK_SIZE, min((ci + 1) * CHUNK_SIZE, n_cells)
        chunk_abs = abs_idx_sorted[lo:hi]
        chunk_sz = hi - lo
        chunk_gcod = gene_codes[lo:hi]
        chunk_pcod = guide_codes[lo:hi]

        # load into mem
        X = adata_backed.X[chunk_abs, :]
        if issparse(X):
            X = X.toarray()
        X = X.astype(np.float32, copy=False)

        # normalize
        row_sums = X.sum(axis=1, keepdims=True)
        np.maximum(row_sums, 1e-10, out=row_sums)
        X *= 10_000.0 / row_sums
        np.log1p(X, out=X)

        # matrix math tricks to sum up the normalized value for each guide id
        rows = np.arange(chunk_sz, dtype=np.int32)
        gene_ind = csr((np.ones(chunk_sz, np.float32), (rows, chunk_gcod)),
                       shape=(chunk_sz, n_gg))
        guide_ind = csr((np.ones(chunk_sz, np.float32), (rows, chunk_pcod)),
                        shape=(chunk_sz, n_pg))
        gene_sum += gene_ind.T @ X
        guide_sum += guide_ind.T @ X
        gene_cnt  += np.bincount(chunk_gcod, minlength=n_gg)
        guide_cnt += np.bincount(chunk_pcod, minlength=n_pg)

        nt_sel = (chunk_gcod == nt_code)
        if nt_sel.any():
            Xnt = X[nt_sel].astype(np.float64)
            nt_n += Xnt.shape[0]
            nt_s += Xnt.sum(axis=0)
            nt_s2 += (Xnt ** 2).sum(axis=0)

        if (ci + 1) % max(1, n_chunks // 10) == 0 or (ci + 1) == n_chunks:
            print(f"    chunk {ci+1}/{n_chunks}  ({hi:,} cells)", flush=True)

        del X
        gc.collect()

    adata_backed.file.close()
    del adata_backed
    gc.collect()

    pb_gene = pd.DataFrame(
        (gene_sum / np.maximum(gene_cnt, 1)[:, None]).astype(np.float32),
        index=pd.Index(gene_list, name="gene"),
        columns=var.index,
    )
    pb_guide = pd.DataFrame(
        (guide_sum / np.maximum(guide_cnt, 1)[:, None]).astype(np.float32),
        index=pd.Index(guide_list, name=obs_guide_col),
        columns=var.index,
    )
    del gene_sum, guide_sum
    gc.collect()

    pb_gene["_n_cells"] = gene_cnt.tolist()
    pb_guide["_n_cells"] = guide_cnt.tolist()

    if nt_n > 1:
        nt_mean = nt_s / nt_n
        nt_var = np.maximum(nt_s2 / nt_n - nt_mean ** 2, 0) * nt_n / (nt_n - 1)
        nt_std = np.sqrt(nt_var)
    else:
        nt_std = np.ones(n_genes, dtype=np.float64)
    var["std"] = np.where(nt_std == 0, 1e-9, nt_std)

    # drop guides with too few cells — their pseudobulks are too noisy for
    # reliable fingerprinting and off-target log2FC estimation
    min_cells = P["min_cells_per_guide_vector"]
    low_cell_mask = pb_guide["_n_cells"] < min_cells
    excluded_df = (pb_guide.loc[low_cell_mask, "_n_cells"]
                   .rename("n_cells")
                   .reset_index()
                   .rename(columns={obs_guide_col: "guide_id"}))
    pb_guide = pb_guide.loc[~low_cell_mask].copy()

    kept_guides = set(pb_guide.index)
    obs = obs[obs[obs_guide_col].isin(kept_guides)].copy()

    n_excluded = int(low_cell_mask.sum())
    print(f"  Excluded {n_excluded:,} guides with < {min_cells} cells "
          f"({len(pb_guide):,} guides retained)", flush=True)

    return pb_gene, pb_guide, var, obs, excluded_df


def compute_neighbors_perturb(pb_guide: pd.DataFrame, pb_gene: pd.DataFrame, obs: pd.DataFrame, var: pd.DataFrame) -> dict:
    """
    Find N_NEIGHBORS nearest neighbor guides for each guide using
    transcriptional fingerprint correlation on pb_guide.

    Unlike the gene-level pipeline (which averages all guides per gene),
    this computes a fingerprint for each individual guide.  Guides that did
    not knock down their intended target will occupy a different region of
    fingerprint space than effective guides, and their off-target candidates
    will reflect their actual transcriptional neighborhood.

    Returns
    -------
    guide_neighbor_map : dict  guide_id -> [guide_id, ...]  (n_neighbors entries)
    """
    nt_label = P["nt_label"]
    obs_guide_col = P["obs_guide_col"]
    obs_gene_col = P["obs_gene_col"]
    n_top_features = P["n_top_features"]
    n_neighbors = P["n_neighbors"]
    output_csv = P["output_csv"]

    nt_guide_ids = set(
        obs.loc[obs[obs_gene_col].astype(str) == nt_label, obs_guide_col].astype(str)
    )
    feat_cols = [c for c in pb_guide.columns if not c.startswith("_")]
    nt_mean = pb_gene.loc[nt_label, feat_cols].values.astype(np.float64)
    feat_std = pb_gene[feat_cols].std(axis=0).values
    feat_std = np.where(feat_std == 0, 1e-9, feat_std)

    # fp stands for "fingerprint" borrowing from Grabski et al
    # intuition is that we are most interested in transcriptional shift from baseline
    # so the delta bw each perturbation and the non-targeting mean is taken
    fp_vals = ((pb_guide[feat_cols].values - nt_mean[np.newaxis, :])
               / feat_std[np.newaxis, :])
    fp_df = pd.DataFrame(fp_vals, index=pb_guide.index, columns=feat_cols)
    fp_df = fp_df.loc[~fp_df.index.isin(nt_guide_ids)].fillna(0)
    fp_gene_vals = ((pb_gene[feat_cols].values - nt_mean[np.newaxis, :])
                    / feat_std[np.newaxis, :])
    fp_gene_for_var = (pd.DataFrame(fp_gene_vals, index=pb_gene.index, columns=feat_cols)
                       .drop(index=nt_label, errors="ignore")
                       .fillna(0))
    feat_var = fp_gene_for_var.var(axis=0)
    top_feats = feat_var.nlargest(n_top_features).index
    fp_sub = fp_df[top_feats].values.astype(np.float64)

    guide_ids = fp_df.index.tolist()
    n_perturb = fp_sub.shape[0]

    # generate a low dim embedding to compute guide neighborhoods
    # n_components is a somewhat arbitrary choice here.
    # some heuristic where n_components ~ number of distinct transcriptionally detected pathways might be ideal
    pca = PCA(n_components=NPCS, random_state=42)
    fp_pca = pca.fit_transform(fp_sub)
    # save fp_pca to a DataFrame for easier debugging and potential downstream use
    fp_pca_df = pd.DataFrame(fp_pca, index=guide_ids, columns=[f"PC{i+1}" for i in range(fp_pca.shape[1])])
    guide_fingerprint_pca_path = output_csv.with_name(output_csv.stem + "_guide_fingerprint_pca.csv")
    fp_pca_df.to_csv(Path(P["output_dir"]) / guide_fingerprint_pca_path)
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1, metric="euclidean", n_jobs=-1) # add one bc nearest neighbor is self
    nn.fit(fp_pca)
    _, indices = nn.kneighbors(fp_pca)

    guide_neighbor_map = {}
    for i, g in enumerate(guide_ids):
        guide_neighbor_map[g] = [guide_ids[idx] for idx in indices[i, 1:]]

    return guide_neighbor_map


def compute_log2fc(pb_guide: pd.DataFrame, pb_gene: pd.DataFrame) -> pd.DataFrame:
    """
    log2fc[guide, gene] = log2((mean_guide + pc) / (mean_NT + pc))
    A pseudocount (log2fc_pseudocount) is added to both numerator and
    denominator to prevent extreme values when either expression is near zero.
    """
    nt_label = P["nt_label"]
    pc = float(P.get("log2fc_pseudocount", 0.1))
    feat_cols = [c for c in pb_guide.columns if not c.startswith("_")]

    nt_mean = pb_gene.loc[nt_label, feat_cols].values.astype(np.float64)
    pair_mat = pb_guide[feat_cols].values.astype(np.float64)
    log2fc = np.log2((pair_mat + pc) / (nt_mean[np.newaxis, :] + pc))

    return pd.DataFrame(log2fc, index=pb_guide.index, columns=feat_cols)


def compute_log2fc_alt(pb_guide: pd.DataFrame, pb_gene: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Computes log2 fold change and returns three DataFrames:
    1. Log2FC values
    2. Mean expression for the numerator (guide-level)
    3. Mean expression for the denominator (NT-level)
    """
    nt_label = P["nt_label"]
    pc = float(P.get("log2fc_pseudocount", 0.1))
    feat_cols = [c for c in pb_guide.columns if not c.startswith("_")]

    # Extract means
    nt_mean = pb_gene.loc[nt_label, feat_cols].values.astype(np.float64)
    pair_mat = pb_guide[feat_cols].values.astype(np.float64)
    
    # Perform the LFC calculation
    # numerator = (pair_mat + pc), denominator = (nt_mean + pc)
    log2fc = np.log2((pair_mat + pc) / (nt_mean[np.newaxis, :] + pc))

    # Construct the DataFrames
    df_lfc = pd.DataFrame(log2fc, index=pb_guide.index, columns=feat_cols)
    
    # Numerator DataFrame (guide mean expression)
    df_num = pd.DataFrame(pair_mat, index=pb_guide.index, columns=feat_cols)
    
    # Denominator DataFrame (NT mean expression broadcasted to match the guide index)
    # We use np.tile or broadcasting to ensure it matches the shape of the guides
    nt_broadcasted = np.tile(nt_mean, (pair_mat.shape[0], 1))
    df_den = pd.DataFrame(nt_broadcasted, index=pb_guide.index, columns=feat_cols)

    return df_lfc, df_num, df_den


def detect_knockdowns(
    pb_guide: pd.DataFrame,
    pb_gene: pd.DataFrame,
    obs:     pd.DataFrame,
    var:     pd.DataFrame):
    """
    For every guide, test whether it knocks down its intended target gene.
    A gene is 'validated' if >= MIN_KD_PAIRS of its guides show knockdown.

    Cache lives in the shared base directory so it is reused if the gene-level
    pipeline has already run.

    Returns
    -------
    kd_df           : DataFrame of all guide by target KD assessments
    validated_genes : set of gene names with validated knockdown
    log2fc_df       : full log2FC DataFrame (guide by features)
    """
    nt_label = P["nt_label"]
    obs_guide_col = P["obs_guide_col"]
    obs_gene_col = P["obs_gene_col"]
    kd_lfc_thresh = P["kd_lfc_thresh"]
    min_kd_pairs = P["min_kd_pairs"]

    log2fc_df, numerator_df, denominator_df = compute_log2fc_alt(pb_guide, pb_gene)

    gene_names = var["gene_name"] if "gene_name" in var.columns else var.index
    gene_name_to_ensembl = dict(zip(gene_names, var.index))

    pair_meta = (obs[[obs_guide_col, obs_gene_col]]
                 .drop_duplicates()
                 .set_index(obs_guide_col))

    n_cells_per_pair = pb_guide["_n_cells"]

    records = []
    for sgID, row in pair_meta.iterrows():
        target_gene = row[obs_gene_col]
        if str(target_gene) == nt_label:
            continue

        ensembl = gene_name_to_ensembl.get(str(target_gene))
        if ensembl is None or ensembl not in log2fc_df.columns:
            continue
        if sgID not in log2fc_df.index:
            continue

        lfc_target = log2fc_df.at[sgID, ensembl]
        num_target = numerator_df.at[sgID, ensembl]
        den_target = denominator_df.at[sgID, ensembl]
        records.append({
            "gene":           target_gene,
            "guide_id":       sgID,
            "target_ensembl": ensembl,
            "target_log2fc":  lfc_target,
            "ontarget_expression": num_target,
            "ntc_expression": den_target,
            "n_cells":        n_cells_per_pair.get(sgID, 0),
            "is_knockdown":   lfc_target < -kd_lfc_thresh,
        })

    kd_df = pd.DataFrame(records)

    kd_counts = (kd_df[kd_df["is_knockdown"]]
                       .groupby("gene")["guide_id"]
                       .nunique())
    validated_genes = set(kd_counts[kd_counts >= min_kd_pairs].index)
    kd_df["validated_target"] = kd_df["gene"].isin(validated_genes)

    return kd_df, validated_genes, log2fc_df


def find_offtarget_events_perturb(
    kd_df: pd.DataFrame,
    validated_genes: set,
    guide_neighbor_map: dict,
    log2fc_df: pd.DataFrame,
    var: pd.DataFrame,
    obs: pd.DataFrame) -> pd.DataFrame:
    """
    Identify guides that knock down a gene targeted by one of their guide-level transcriptional neighbors.

    Returns:
    events : DataFrame with one row per (guide by off-target gene) event
    """
    kd_lfc_thresh = P["kd_lfc_thresh"]
    obs_guide_col = P["obs_guide_col"]
    obs_gene_col = P["obs_gene_col"]

    gene_names = var["gene_name"] if "gene_name" in var.columns else var.index
    gene_name_to_ensembl = dict(zip(gene_names, var.index))

    # guide_id → intended target gene
    guide_to_gene = (obs[[obs_guide_col, obs_gene_col]]
                     .drop_duplicates()
                     .set_index(obs_guide_col)[obs_gene_col]
                     .to_dict())

    # look up 'reverse' neighbors using the original nn map
    from collections import defaultdict
    reverse_map: dict = defaultdict(list)
    for g, nbrs in guide_neighbor_map.items():
        for rank, nbr in enumerate(nbrs, 1):
            reverse_map[str(nbr)].append((str(g), rank))

    # loop over all of the guides and see if any of their neighbors target a gene that they knock down
    records = []
    candidate_guides = kd_df[kd_df["validated_target"]]
    for _, row in candidate_guides.iterrows():
        target_gene = row["gene"]
        guide_id = str(row["guide_id"])

        if guide_id not in log2fc_df.index:
            continue

        # forward = guide in its own neighbor list; reverse = guide in another guide's neighbor list
        fwd = [(nbr, rank, "forward")
               for rank, nbr in enumerate(guide_neighbor_map.get(guide_id, []), 1)]
        rev = [(g, rank, "reverse")
               for g, rank in reverse_map.get(guide_id, [])]

        # Accumulate best record per off-target gene, merging direction
        best: dict = {}
        for nbr_guide, rank, direction in fwd + rev:
            offtarget_gene = guide_to_gene.get(str(nbr_guide))
            if not offtarget_gene:
                continue
            if str(offtarget_gene) == str(target_gene): # make sure not to call self-knockdown an off-target event
                continue

            # TODO: probably rename; not actually using ENSEMBL IDs in a lot of cases
            offtarget_ensembl = gene_name_to_ensembl.get(str(offtarget_gene))
            ontarget_ensembl = gene_name_to_ensembl.get(str(target_gene))
            if not offtarget_ensembl or offtarget_ensembl not in log2fc_df.columns:
                continue

            # get log2fc of the on-target guide at the off-target gene and see if it passes threshold
            lfc_ontarget_at_offtarget = log2fc_df.at[guide_id, offtarget_ensembl]
            lfc_ontarget_at_ontarget = log2fc_df.at[guide_id, ontarget_ensembl]

            # DONT USE A CUTOFF IN THIS VERSION - we want to see all of the candidates and their LFCs
            # if lfc_ontarget_at_offtarget >= -kd_lfc_thresh:
            #     continue

            if offtarget_gene not in best:
                best[offtarget_gene] = {
                    "gene":              target_gene,
                    "guide_id":          guide_id,
                    "target_log2fc":     lfc_ontarget_at_ontarget,
                    "offtarget_gene":    offtarget_gene,
                    "offtarget_ensembl": offtarget_ensembl,
                    "offtarget_log2fc":  lfc_ontarget_at_offtarget,
                    "neighbor_rank":     rank,
                    "direction":         direction,
                }
            else: # many neighbors are in the fwd and rev neighborhoods. just save the rank of the closest neighbor on rank
                prev = best[offtarget_gene]
                if prev["direction"] != direction:
                    prev["direction"] = "both"
                if rank < prev["neighbor_rank"]:
                    prev["neighbor_rank"] = rank

        records.extend(best.values())

    events = pd.DataFrame(records)
    events = events.sort_values(["direction", "neighbor_rank"],
                                key=lambda s: s if s.name != "direction"
                                else s.map({"both": 0, "forward": 1, "reverse": 2}))

    print(f"  {len(events):,} off-target events "
          f"({events['gene'].nunique() if not events.empty else 0} target genes, "
          f"{events['guide_id'].nunique() if not events.empty else 0} guides)",
          flush=True)

    return events


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence"""
    return seq.translate(RC_MAP)[::-1]


def load_fai(fai_path: Path) -> dict:
    """Load a .fai index file for a reference genome FASTA"""
    fai = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            fai[chrom] = (int(parts[1]), int(parts[2]),
                          int(parts[3]), int(parts[4]))
    return fai


def extract_genome_region(fai_dict: dict, chrom: str,
                          start: int, end: int,
                          genome_path: Path) -> str:
    """Extract a genomic region from a reference genome FASTA using a .fai index"""
    if chrom not in fai_dict:
        return ""
    chrom_len, fai_offset, linebases, linewidth = fai_dict[chrom]
    start = max(0, start)
    end = min(chrom_len, end)
    if start >= end:
        return ""

    byte_start = fai_offset + (start // linebases) * linewidth + (start % linebases)
    byte_last = fai_offset + ((end - 1) // linebases) * linewidth + ((end - 1) % linebases)

    with open(genome_path, "rb") as fh:
        fh.seek(byte_start)
        raw = fh.read(byte_last - byte_start + 1)

    seq = raw.decode("ascii", errors="ignore").replace("\n", "").replace("\r", "")
    return seq[: end - start].upper()


def find_max_seed_match(full_spacer: str, seq: str,
                        min_len: int = 5) -> tuple[int, int | None, int | None, str | None, str | None, str | None]:
    """
    Find the maximum seed length (from the 3' end of the spacer) that produces
    a seed+NGG match in seq (forward or reverse complement strand).

    Also computes the Hamming distance.

    Returns
    (best_seed_len, hamming_distance, seed_local_pos, seed_strand, genomic_sequence, protospacer_sequence)
    seed_local_pos      : 0-based position within seq where the seed starts (fwd) or
                          where rc_seed starts on the fwd strand (rev); None if no match
    seed_strand         : "+" for a forward-strand match, "-" for a reverse-strand match;
                          None when best_seed_len == 0 (no match found).
    genomic_sequence    : spacer-length genomic sequence at the match site (None if no match)
    protospacer_sequence: the guide spacer sequence used for matching (None if no match)
    hamming_distance is None when best_seed_len == 0.
    """
    best = 0
    best_info = None   # ("fwd"|"rev", seed_start_pos)
    spacer_len = len(full_spacer)
    seq_up = seq.upper()

    for length in range(min_len, spacer_len + 1):
        seed = full_spacer[-length:].upper()
        rc_seed = revcomp(seed)
        m_fwd = re.search(re.escape(seed) + r"[ACGT]GG", seq_up)
        m_rev = re.search(r"CC[ACGT]" + re.escape(rc_seed), seq_up)
        if m_fwd:
            best = length
            best_info = ("fwd", m_fwd.start())
        elif m_rev:
            best = length
            best_info = ("rev", m_rev.start())

    if best == 0 or best_info is None:
        return 0, None, None, None, None, None

    strand, p = best_info
    if strand == "fwd":
        # Seed is at seq_up[p : p+best]; guide aligns from p-(spacer_len-best)
        g_start = p - (spacer_len - best)
        if g_start < 0:
            genomic = "N" * (-g_start) + seq_up[: p + best]
        else:
            genomic = seq_up[g_start : p + best]
        local_pos = p
        strand_char = "+"
    else:
        # Pattern CC[N]rc(seed) at fwd pos p.
        # Full spacer on rev strand = revcomp(seq_up[p+3 : p+3+spacer_len])
        s_end = p + 3 + spacer_len
        site_fwd = (seq_up[p + 3 : s_end] if s_end <= len(seq_up)
                    else seq_up[p + 3 :] + "N" * (s_end - len(seq_up)))
        genomic = revcomp(site_fwd)
        local_pos = p + 3   # 0-based start of rc_seed on the fwd strand
        strand_char = "-"

    hamming = sum(a != b for a, b in zip(full_spacer.upper(), genomic))
    return best, hamming, local_pos, strand_char, genomic, full_spacer.upper()


def get_tss_map() -> dict:
    global _TSS_MAP
    if _TSS_MAP is not None:
        return _TSS_MAP

    upstream = P["promoter_upstream"]
    downstream = P["promoter_downstream"]
    tss_map: dict = {}

    tss_map_csv = P.get("tss_map_csv")
    tss_df = pd.read_csv(tss_map_csv)
    missing_strand = "strand" not in tss_df.columns
    for _, row in tss_df.iterrows():
        gene = row["gene_name"]
        chrom = row["chrom"]
        tss = int(row["tss"])
        strand = row["strand"] if not missing_strand else "+"
        tss_0 = tss - 1
        if strand == "+":
            prom_start, prom_end = tss_0 - upstream, tss_0 + downstream
        else:
            prom_start, prom_end = tss_0 - downstream, tss_0 + upstream
        tss_map[gene] = (chrom, max(0, prom_start), prom_end)
    _TSS_MAP = tss_map
    return tss_map


def get_promoter_coords(gene_sym: str) -> tuple | None:
    return get_tss_map().get(gene_sym)


def load_sgrna_sequences(sgrna_csv: Path) -> dict:
    dataset_type = P.get("dataset_type", "single_guide")
    df = pd.read_csv(sgrna_csv)
    lookup = {}

    if dataset_type == "paired_guide":
        for _, row in df.iterrows():
            pair_id = f"{row['sgID_A']}|{row['sgID_B']}"
            lookup[pair_id] = {
                "seq_a": str(row.get("targeting sequence A", "") or ""),
                "seq_b": str(row.get("targeting sequence B", "") or ""),
                "gene":  row.get("gene", ""),
            }
    else:
        for _, row in df.iterrows():
            guide_id = row["sgRNA"]
            lookup[guide_id] = {
                "seq":  str(row.get("seq", "") or ""),
                "gene": str(row.get("designed_target_gene_name", "") or ""),
            }

    return lookup


# Module-level state populated in each worker process by _seed_worker_init.
# Using an initializer avoids pickling the large cache and lookup dicts on every task.
_W_PROMOTER_CACHE: dict = {}
_W_SGRNA_LKP: dict = {}
_W_DATASET_TYPE: str = "single_guide"


def _seed_worker_init(promoter_cache: dict, sgrna_lkp: dict, dataset_type: str) -> None:
    global _W_PROMOTER_CACHE, _W_SGRNA_LKP, _W_DATASET_TYPE
    _W_PROMOTER_CACHE = promoter_cache
    _W_SGRNA_LKP = sgrna_lkp
    _W_DATASET_TYPE = dataset_type


def _process_seed_row(args: tuple) -> dict:
    """Per-row worker: look up promoter sequence and run seed match."""
    guide_id, offtarget = args
    prom_data = _W_PROMOTER_CACHE.get(str(offtarget))
    entry = _W_SGRNA_LKP.get(guide_id)

    null = {"seed": None, "seed_match_len": 0, "guide_hamming_dist": None,
            "seed_match_chrom": None, "seed_match_start": None, "seed_match_strand": None,
            "genomic_sequence": None, "protospacer_sequence": None}

    if prom_data is None or entry is None:
        return null
    prom_seq, coords = prom_data
    if not prom_seq:
        return null

    if _W_DATASET_TYPE == "paired_guide":
        seq_a = entry.get("seq_a", "")
        seq_b = entry.get("seq_b", "")
        match_a, hamm_a, lpos_a, str_a, gseq_a, pseq_a = (find_max_seed_match(seq_a, prom_seq)
                                                             if seq_a else (0, None, None, None, None, None))
        match_b, hamm_b, lpos_b, str_b, gseq_b, pseq_b = (find_max_seed_match(seq_b, prom_seq)
                                                             if seq_b else (0, None, None, None, None, None))
        if match_a >= match_b:
            best_match, best_hamm, best_seq, best_lpos, best_str, best_gseq, best_pseq = (
                match_a, hamm_a, seq_a, lpos_a, str_a, gseq_a, pseq_a)
        else:
            best_match, best_hamm, best_seq, best_lpos, best_str, best_gseq, best_pseq = (
                match_b, hamm_b, seq_b, lpos_b, str_b, gseq_b, pseq_b)
        seed = best_seq[-best_match:].upper() if best_match > 0 else None
    else:
        seq = entry.get("seq", "")
        best_match, best_hamm, best_lpos, best_str, best_gseq, best_pseq = (
            find_max_seed_match(seq, prom_seq) if seq else (0, None, None, None, None, None))
        seed = seq[-best_match:].upper() if best_match > 0 else None

    if best_match > 0 and coords is not None and best_lpos is not None:
        chrom, prom_start, _ = coords
        seed_chrom  = chrom
        seed_start  = prom_start + best_lpos
        seed_strand = best_str
    else:
        seed_chrom = seed_start = seed_strand = None

    return {"seed": seed, "seed_match_len": best_match,
            "guide_hamming_dist": best_hamm,
            "seed_match_chrom": seed_chrom,
            "seed_match_start": seed_start,
            "seed_match_strand": seed_strand,
            "genomic_sequence": best_gseq,
            "protospacer_sequence": best_pseq}


def search_seeds_in_promoters(events: pd.DataFrame,
                               sgrna_csv: Path) -> pd.DataFrame:
    fai = load_fai(P["genome_fai"])
    sgrna_lkp = load_sgrna_sequences(sgrna_csv)
    dataset_type = P.get("dataset_type", "single_guide")
    n_workers = int(P.get("n_seed_workers", 30))
    events_r = events.reset_index(drop=True)
    n = len(events_r)

    # Pre-fetch all unique promoter sequences sequentially.
    # Doing this before forking workers avoids concurrent seeks on the same FASTA file.
    unique_genes = events_r["offtarget_gene"].unique()
    print(f"  Pre-fetching promoter sequences for {len(unique_genes):,} unique genes…", flush=True)
    tss_map = get_tss_map()
    promoter_cache: dict = {}
    for gene in unique_genes:
        coords = tss_map.get(str(gene))
        if coords is None:
            promoter_cache[str(gene)] = None
        else:
            chrom, pstart, pend = coords
            seq = extract_genome_region(fai, chrom, pstart, pend, P["genome_fa"])
            promoter_cache[str(gene)] = (seq, coords) if seq else None

    args_list = [(row["guide_id"], row["offtarget_gene"]) for _, row in events_r.iterrows()]
    chunksize = max(1, n // (n_workers * 8))
    print(f"  Seed search: {n:,} events, {n_workers} workers, chunksize={chunksize}", flush=True)

    with ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_seed_worker_init,
        initargs=(promoter_cache, sgrna_lkp, dataset_type),
    ) as pool:
        result_rows = list(pool.map(_process_seed_row, args_list, chunksize=chunksize))

    result_df = pd.DataFrame(result_rows)
    return pd.concat([events_r, result_df], axis=1)


def main():
    global P, _TSS_MAP

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--params", required=True,
                        help="Path to params.yaml")
    args = parser.parse_args()

    params_path = Path(args.params)
    if not params_path.exists():
        print(f"ERROR: params file not found: {params_path}", file=sys.stderr)
        sys.exit(1)

    P = load_params(params_path)
    _TSS_MAP = None

    perturb_dir = Path(P["output_dir"])
    perturb_dir.mkdir(parents=True, exist_ok=True)

    output_csv = perturb_dir / Path(P["output_csv"]).name

    print(f"=== Off-target pipeline (perturbation neighborhoods): "
          f"{P.get('dataset_name', params_path.stem)} ===", flush=True)
    print(f"  adata:     {P['adata_path']}", flush=True)
    print(f"  sgrna_csv: {P['sgrna_csv']}", flush=True)
    print(f"  output:    {output_csv}", flush=True)
    print(f"  nt_label:  {P['nt_label']}", flush=True)
    print(f"  dataset:   {P.get('dataset_type', 'single_guide')}", flush=True)
    print("", flush=True)

    t_total = time.time()

    print("Computing pseudobulk profiles per guide (Step 1)…")
    pb_gene, pb_guide, var, obs, excluded_df = compute_pseudobulk()
    print(f"  pb_gene: {pb_gene.shape}  pb_guide: {pb_guide.shape}\n", flush=True)

    excluded_csv = output_csv.with_name(output_csv.stem + "_excluded_guides.csv")
    excluded_df.to_csv(excluded_csv, index=False)
    print(f"  Excluded guides written to: {excluded_csv}\n", flush=True)

    print("Computing guide-level neighbors (Step 2)…")
    guide_neighbor_map = compute_neighbors_perturb(pb_guide, pb_gene, obs, var)
    print(f"  {len(guide_neighbor_map):,} guides mapped\n", flush=True)

    print("Detecting knockdowns per guide (Step 3)…")
    kd_df, validated_genes, log2fc_df = detect_knockdowns(pb_guide, pb_gene, obs, var)
    print(f"  Validated target genes: {len(validated_genes):,}\n", flush=True)

    ontarget_lfc_csv = output_csv.with_name(output_csv.stem + "_ontarget_lfcs.csv")
    kd_df.to_csv(ontarget_lfc_csv, index=False)
    print(f"  On-target LFCs written to: {ontarget_lfc_csv}\n", flush=True)

    print("Finding guide-level off-target events (Step 4)…")
    events = find_offtarget_events_perturb(kd_df, validated_genes, guide_neighbor_map, log2fc_df, var, obs)

    final = search_seeds_in_promoters(events, P["sgrna_csv"])
    cols_out = [
        "gene", "guide_id", "target_log2fc",
        "offtarget_gene", "offtarget_ensembl", "offtarget_log2fc",
        "neighbor_rank", "direction",
        "seed", "seed_match_len", "guide_hamming_dist",
        "seed_match_chrom", "seed_match_start", "seed_match_strand",
        "genomic_sequence", "protospacer_sequence",
    ]
    final = final[[c for c in cols_out if c in final.columns]]
    final = final.sort_values(["seed_match_len", "offtarget_log2fc"], ascending=[False, True])

    final.to_csv(output_csv, index=False)
    print(f"\n=== DONE in {time.time()-t_total:.0f}s ===", flush=True)
    print(f"Output: {output_csv}  ({len(final):,} rows)", flush=True)
    has_match = (final["seed_match_len"] >= 12).sum() if "seed_match_len" in final.columns else 0
    print(f"  Rows with seed match in off-target promoter: {has_match:,}", flush=True)
    print(final[final["seed_match_len"] >= 12].to_string())


if __name__ == "__main__":
    main()
