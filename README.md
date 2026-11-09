# aSBM-Louvain

R implementation and experiment code for *aSBM-Louvain: Attributed Network Clustering with Louvain-Based Block Structure Learning*.

From the repository root, supply a symmetric, binary adjacency matrix `A` with a zero diagonal and at least one edge. The numeric attribute matrix `X` must have one row per node in the same order.

```r
install.packages(c("igraph", "mvtnorm", "Matrix"))
library(Matrix)
source("asbm_louvain.R")
fit <- aSBM_Louvain_revised(A, X)
labels <- fit$partition
```

[Algorithm source](asbm_louvain.R) · [Experiment commands](reproduction/README.md) · [STARmap, MERFISH, and baselines](reproduction/experiments/realdata_current/README.md)
