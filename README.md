# aSBM-Louvain

Code to reproduce all experiments from the paper *"aSBM-Louvain: Attributed Network Clustering with Louvain-Based Block Structure Learning"*.

## Overview

This repository implements **aSBM-Louvain**, a method for clustering attributed networks that combines the attributed Stochastic Block Model (aSBM) with Louvain community detection. The Louvain step replaces the slow likelihood-based partition search of Pure-aSBM (Stanley et al., 2019), reducing per-run time from hours to under a second.

The experiments cover:

- **Synthetic benchmarks** — A scenarios (structure degradation) and B scenarios (attribute noise), comparing Louvain, KMeans, Pure-aSBM, and aSBM-Louvain.
- **Real-data benchmarks** — STARmap (1,207 cells, 16 cell types) and MERFISH (5,202 cells, 19 cell types) spatial transcriptomics datasets, adding SpaGCN and BACT as baselines.

## Repository Structure

```
├── asbmlou.R                  # Main aSBM-Louvain algorithm
├── FitAttribute.R             # Pure aSBM baseline (Stanley et al., 2019)
├── run_simulations.R          # Unified synthetic experiment runner
├── plot_results.py            # Generate publication figures from CSV results
├── starmap_PCA_scores_50PCs.csv       # STARmap precomputed PCA (50 PCs)
├── coord_and_pc_19.RData              # MERFISH precomputed PCs + coordinates
│
└── realdata_experiments/
    ├── realdata_submission_helpers.R           # Shared R: data loading, KNN, output
    ├── spagcn_submission_helpers.py            # Shared Python: SpaGCN utilities
    ├── starmap_asbmlouvain_pca_coord_input_10runs.R    # aSBM-Louvain
    ├── starmap_louvain_asbm_pca_coord_input_10runs.R   # Louvain + Pure-aSBM
    ├── starmap_bact_pca_coord_input_10runs.R           # BACT
    ├── starmap_spagcn_expression_coord_input_10runs.py # SpaGCN
    ├── merfish_asbmlouvain_pc_coord_input_10runs.R     # aSBM-Louvain
    ├── merfish_louvain_asbm_pc_coord_input_10runs.R    # Louvain + Pure-aSBM
    ├── merfish_bact_pc_coord_input_10runs.R            # BACT
    ├── merfish_spagcn_expression_coord_input_10runs.py # SpaGCN
    ├── scanpy/                              # Local scanpy overrides for SpaGCN
    └── results/                             # Output directory (created at runtime)
```

## Dependencies

### R packages

```r
install.packages(c("igraph", "mvtnorm", "mclust", "FNN"))
```

- `igraph` (>= 1.6.0)
- `mvtnorm`
- `mclust`
- `FNN`
- `parallel` (built-in)
- `BACT` — only required for BACT baseline scripts; install from [BACT](https://github.com/username/BACT)

### Python packages

```bash
pip install matplotlib numpy pandas scanpy SpaGCN torch scikit-learn anndata
```

- `SpaGCN` — only required for SpaGCN baseline scripts
- `anndata` — only required for MERFISH SpaGCN script
- `scanpy` — the repo includes local overrides in `realdata_experiments/scanpy/`

## Synthetic Experiments

All synthetic experiments are driven by `run_simulations.R` from the repository root.

### Quick Start

```bash
# Run all experiments (expect several days for Pure-aSBM)
Rscript run_simulations.R all

# Run only fast baselines (Louvain, KMeans, aSBM-Louvain) -- finishes in minutes
Rscript run_simulations.R baselines A 3.0    # A scenarios, sigma=3.0
Rscript run_simulations.R baselines B         # B scenarios, r=4.0

# Run Pure-aSBM (slow -- hours to days)
Rscript run_simulations.R pure A 3.0          # A scenarios, sigma=3.0
Rscript run_simulations.R pure B 4.0          # B scenarios, r=4.0
```

### Experiment Settings

**A Scenarios (Structure Degradation)**
- N = 1000 (A1/A2) or 5000 (A3)
- K = 5 (A1/A3) or 10 (A2)
- D = 50, target degree ~ 20
- sigma = 3.0 (fixed attribute noise)
- r = 7.5, 6.0, 4.5, 3.0, 1.5, 1.0 (decreasing structure strength)
- A3 only runs r >= 4.5 due to computational infeasibility at low r

**B Scenarios (Attribute Noise)**
- N = 1000 (B1/B2) or 5000 (B3)
- K = 5 (B1/B3) or 10 (B2)
- D = 50, target degree ~ 20
- r = 4.0 (fixed structure strength)
- sigma = 0.5, 1.0, ..., 5.0 (increasing attribute noise)

### Methods Compared

| Method | Description |
|---|---|
| Louvain | Graph topology only |
| KMeans | Node attributes only |
| Pure aSBM | Likelihood-based attributed SBM (Stanley et al., 2019) |
| aSBM-Louvain | Proposed method |

### Performance & Runtime

- All methods: 5 independent repetitions per configuration
- Pure-aSBM: 10-hour timeout per single run; timed-out runs report current partition
- Metric: Adjusted Rand Index (ARI)

| Experiment | Approximate Time |
|---|---|
| Baselines (A+B, all scenarios) | ~5 minutes |
| Pure-aSBM A1/A2 (sigma=3.0) | ~6-12 hours |
| Pure-aSBM A3 (high r only) | ~12-24 hours |
| Pure-aSBM B1/B2 (r=4.0) | ~1-3 hours |
| Pure-aSBM B3 (r=4.0) | ~2-4 days |

### Output

Each run produces:
- `*_results.csv`: per-repetition raw results
- `*_point_summary.csv`: mean +/- SD per parameter point

## Real-Data Experiments

### Data Preparation

The repository includes two preprocessed data files:

- `starmap_PCA_scores_50PCs.csv` — STARmap PCA scores (precomputed 50 PCs)
- `coord_and_pc_19.RData` — MERFISH principal components and spatial coordinates

**Additional data files** must be downloaded from the original publications and placed in the repository root:

**STARmap** (Wang et al., 2018, *Science*):
- `starmap_spatial.csv` — X, Y coordinates per cell
- `Annotation_STARmap_20180505_BY3_1k.txt` — ground truth cell-type annotation (tab-separated, `Annotation` column)
- `starmap_counts.csv` — gene expression count matrix (SpaGCN only)

**MERFISH** (Moffitt et al., 2018, *Science*):
- `cell_type_annotation_0.19.csv` — ground truth cell-type annotation (first column)
- `merfish19_counts.csv` — gene expression count matrix (SpaGCN only)
- `merfish19_spatial.csv` — spatial coordinates (SpaGCN only)

**Important:** All data files should be placed in the repository root (same directory as `asbmlou.R`). Scripts load data from the project root directory.

### Running the Experiments

All scripts must be run **from the repository root** so that `source()` calls resolve correctly.

**aSBM-Louvain:**

```bash
Rscript realdata_experiments/starmap_asbmlouvain_pca_coord_input_10runs.R
Rscript realdata_experiments/merfish_asbmlouvain_pc_coord_input_10runs.R
```

**Louvain + Pure-aSBM:**

```bash
Rscript realdata_experiments/starmap_louvain_asbm_pca_coord_input_10runs.R
Rscript realdata_experiments/merfish_louvain_asbm_pc_coord_input_10runs.R
```

**BACT:**

```bash
Rscript realdata_experiments/starmap_bact_pca_coord_input_10runs.R
Rscript realdata_experiments/merfish_bact_pc_coord_input_10runs.R
```

**SpaGCN:**

```bash
python realdata_experiments/starmap_spagcn_expression_coord_input_10runs.py
python realdata_experiments/merfish_spagcn_expression_coord_input_10runs.py
```

### Methods x Datasets

| Method | STARmap | MERFISH | Input |
|---|---|---|---|
| Louvain | 10 runs | 10 runs | PCA/PC KNN graph |
| aSBM (Pure) | 10 runs | 10 runs | PCA/PC KNN graph + spatial attributes |
| aSBM-Louvain | 10 runs | 10 runs | PCA/PC KNN graph + spatial attributes |
| BACT | 10 runs | 10 runs | PCA/PC matrix + spatial coordinates |
| SpaGCN | 10 runs | 10 runs | Expression matrix + spatial coordinates |

### Output

Each script writes results to `realdata_experiments/results/`:
- `{prefix}_metrics.csv` — per-run ARI, NMI, K_Found
- `{prefix}_labels.csv` — cluster labels per cell per run
- `{prefix}_summary.csv` — mean, SD, median across runs

## Figures

After obtaining result CSV files (either by re-running experiments or using pre-computed results), generate publication figures:

```bash
python plot_results.py
```

Figures are saved as PNG files at 300 DPI.

## Pre-computed Results

Pre-computed result CSV files for the synthetic experiments can be provided upon request. Place them in the repository root and run `plot_results.py` directly without re-running the experiments.

The expected CSV files follow the naming convention:
- `A_baselines_sigma3.0_results.csv`, `A_baselines_sigma3.0_point_summary.csv`
- `B_baselines_r4_results.csv`, `B_baselines_r4_point_summary.csv`
- `pure_asbm_A_sigma3.0_results.csv`, `pure_asbm_A_sigma3.0_point_summary.csv`
- `pure_asbm_B_r4_results.csv`, `pure_asbm_B_r4_point_summary.csv`

## Citation

If you use this code in your research, please cite:

```
[Paper citation to be added]
```

## References

- Stanley, N., et al. (2019). Stochastic Block Models with Continuous Attributes. *Applied Network Science*.
- Wang, X., et al. (2018). Three-dimensional intact-tissue sequencing of single-cell transcriptional states. *Science*.
- Moffitt, J. R., et al. (2018). Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. *Science*.
