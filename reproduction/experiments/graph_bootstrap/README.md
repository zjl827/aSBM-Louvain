# S12 conditional graph bootstrap

From the release root, install R packages `igraph`, `Matrix`, `mclust` and `mvtnorm` and Python package `numpy`. Prepare the two primary-orientation real inputs using the [real-data suite](../realdata_current/README.md). The eight synthetic sources are regenerated directly using the first fixed Gaussian replicate in each source condition. The smoke run is synthetic and can run before preparing the real inputs.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript experiments/graph_bootstrap/run.R smoke
Rscript experiments/graph_bootstrap/run.R
python3 experiments/graph_bootstrap/summarize.py
```

The full run fits four methods on ten source inputs and on 50 paired draws per source. Pure-aSBM uses the vectorized R implementation. The two arms share a sampled graph. Graph-only holds the original attributes fixed, and joint draws attributes from the fitted within-community Gaussians. Each source's Proposed partition anchors the plug-in generator and agreement metrics. The resulting statistics characterize conditional stability under that generator. Pairwise comparisons among bootstrap partitions are descriptive because pairs share partitions.

The full output is `results/graph_bootstrap/full`. The smoke output is `results/graph_bootstrap/smoke` and contains one synthetic source and one pair of draws. A `preflight` run fits the original-data anchors and checks the generator. Set `REVISION_BOOTSTRAP_OUT` to a fresh output directory when the code or inputs differ from saved checkpoints. Use `REVISION_BOOTSTRAP_WORKERS=1` for serial execution and `REVISION_REALDATA_DIR` for prepared real-data inputs elsewhere. Run `summarize.py --input-dir DIR --expected-sources 1 --expected-b 1` to validate and summarize a smoke run.
