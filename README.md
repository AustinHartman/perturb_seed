# Perturb Seed

## Overview

This repository contains the pipeline for identifying seed-driven off-target effects in GWPS Perturb-seq experiments, as described in the preprint: **"Systematic identification of seed-driven off-target effects in Perturb-seq experiments"**

The pipeline nominates off-target CRISPR guide candidates by combining seed alignments at off-target promoters with measured transcriptional repression at off-target genes. Please read the preprint for more details on the method.

## Repository Structure

- `off_target_pipeline_perturb.py` — Main pipeline script for off-target candidate generation
- `tss_map.csv` — Transcription start sites coordinates for hg38
- `guide_libraries/` — sgRNA library CSVs for each analyzed dataset
  - `sgrna_bradu.csv` — K562 cells (Bradu & Blair dataset)
  - `sgrna_replogle.csv` — K562 cells (Replogle dataset)
  - `sgrna_xaira.csv` — HCT116 cells (Xaira dataset)
  - `sgrna_zhu.csv` — CD4+ T cells (Zhu & Dann dataset)
- `input_yamls/` — YAML configuration files for each dataset and condition. AnnData objects must be downloaded separately (links in data availability section of preprint).
- `environment.yml` - Conda environment used for analyses

## Output

`*_off_target_candidates.csv`: Sorted by seed match length and off-target gene log2 fold-change so higher likelihood off-targets will be found at the beginning of the file with very low likelihood off-target events at the bottom

#### Output Columns

- `gene` — Gene the guide is designed to target
- `guide_id` — Guide identifier
- `target_log2fc` — Pseudobulk log2 fold-change of `guide_id` at `gene`
- `offtarget_gene` — Gene the guide may have off-target activity at
- `offtarget_log2fc` — Pseudobulk log2 fold-change of `guide_id` at `offtarget_gene`. **NOTE:** currently all events are recorded in output and not filtered by log2 fold-change
- `neighbor_rank` — Lowest rank in KNN neighborhood (of the direct and asymmetrical neighbors -- see paper methods for more information on neighbor finding)
- `direction` — Direct neighbor, asymmetrical neighbor, or both
- `seed` — Guide seed sequence which aligns to the off-target site adjacent to an `NGG` PAM sequence
- `seed_match_len` — Length of the `seed` alignment
- `guide_hamming_dist` — Hamming distance of the full 20 bp alignment even though only the first `seed_match_len` bases match exactly (additional matches might increase the strength of association with an off-target locus)
- `seed_match_chrom` — Chromosome of seed match
- `seed_match_start` — Start position of seed match
- `seed_match_strand` — Strand of seed match

## Installation

Clone the repository and create the conda environment from `environment.yml`:

```bash
git clone https://github.com/AustinHartman/perturb_seed.git
cd perturb_seed
conda env create -f environment.yml
conda activate vcp
```

All analyses described in the preprint were performed on a GCP machine with 256 GB memory and 32 CPUs.

**Expected installation time:** Creating the conda environment takes ~10 minutes.

## Usage

```bash
python off_target_pipeline_perturb.py --yaml input_yamls/params_replogle.yaml
```

Each YAML file specifies the input h5ad AnnData object, sgRNA library, reference genome paths, and other parameters.

**Expected runtime:** On a machine comparable to the one described above (256 GB memory, 32 CPUs), a single dataset typically completes in ~10 minutes, though this depends on the on-disk matrix format of the input AnnData object (e.g. dense vs. sparse, CSR vs. CSC), which affects how quickly the data can be loaded and processed.

#### Parameters:

- `n_top_features` — Number of top variable genes for fingerprinting
- `n_neighbors` — KNN neighbors in PCA space
- `kd_lfc_thresh` — Log2FC threshold for on-target knockdown
- `min_cells_per_guide_vector` — Minimum cells per guide to retain
- `promoter_upstream` — Bases upstream of TSS to search for seed alignments (250-2,000 are seem to be reasonable values, but it likely also depends on the effector tethered to dCas)
- `promoter_downstream` — Bases downstream of TSS to search for seed alignments
- `dataset_type` — `single_guide` meaning each vector delivered to cells contains a single guide or `paired_guide` which means two guides targeting the same gene are delivered to cells together

The two described `dataset_type` values covered the experimental designs for the four datasets we analyzed in the preprint, but feel free to open an issue or email me if there are additional experimental designs I can incorporate.

## Resources

- 📄 [Read the preprint](https://www.biorxiv.org/content/10.64898/2026.03.27.714658v2)
- 🔧 [GitHub repository to reproduce analyses in the preprint](https://github.com/AustinHartman/perturb_seed_reproduce)
- 🌐 [Web-app to query guides for seed alignments near transcription start sites](https://crispr-seed-finder.vercel.app/)

## Data & Contact

Input data for each plot includes the AnnData object. See the preprint's data availability section for download links.

For additional output files or scripts feel free to email me at [hartmana@stanford.edu](mailto:hartmana@stanford.edu).
