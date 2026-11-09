# aSBM-Louvain

Code for the experiments in *aSBM-Louvain: Attributed Network Clustering with Louvain-Based Block Structure Learning*. The method combines an attributed stochastic block model with Louvain community detection.

## Run aSBM-Louvain

From the repository root, load the algorithm and fit it to an undirected graph.
`A` is a symmetric, binary adjacency matrix with a zero diagonal and at least
one edge. `X` is a numeric matrix with one row per node in the same order. [Algorithm source](asbm_louvain.R)

```r
install.packages(c("igraph", "mvtnorm", "Matrix"))
source("asbm_louvain.R")
fit <- aSBM_Louvain_revised(A, X)
labels <- fit$partition
```

For sparse adjacency matrices, load the `Matrix` package before fitting.

## Reproducing the experiments

Run the commands in the experiment guide from `reproduction/`. The guide lists dependencies, data sources, and experiment commands. [Setup and run instructions](reproduction/README.md)

- The balanced simulations cover scenarios A1–A3 and B1–B3, with scripts for generating data, running the methods, and plotting the results. [Balanced simulations](reproduction/README.md#balanced-simulations)
- The supplementary experiments are indexed by section number. [Supplementary experiments](reproduction/README.md#supplementary-experiments)
- The STARmap and MERFISH workflow prepares the input data and runs aSBM-Louvain, vectorized Pure-aSBM, Louvain, and K-means. The experiment guide includes the matched timing commands. [Real-data experiments](reproduction/experiments/realdata_current/README.md) · [Timing experiments](reproduction/README.md#s10s12-real-data-experiments-and-graph-bootstrap)
- The BACT and SpaGCN scripts use their own input and dependency setup. [Baseline scripts](realdata_experiments/) · [BACT and SpaGCN setup](reproduction/experiments/realdata_current/README.md#bact-and-spagcn)
