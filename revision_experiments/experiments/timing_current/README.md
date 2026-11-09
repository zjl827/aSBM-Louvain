# S10 matched runtime

From the release root, install R packages `igraph`, `Matrix`, `mclust` and `mvtnorm`. Prepare the two primary-orientation input RDS files using the [real-data suite](../realdata_current/README.md). The smoke run uses a synthetic input and can run independently of those files.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript experiments/timing_current/run.R smoke
Rscript experiments/timing_current/run.R
```

The formal run uses 70 matched cases, with ten cases at each of five simulation sizes and ten repeats on each of the two fixed real-data inputs. The bundled task grid records all original seeds. The runner executes serially, rotates the order of the two methods, warms both paths separately, shares estimated-K K-means labels and downstream RNG, and includes the shared initialization time in each method's total. The timed fits use the revised aSBM–Louvain and the exact R-only Pure-aSBM. Time measurement excludes simulation, loading and warmup. The summaries include successful, uncapped fits and the raw table retains all attempted runs and stop statuses. Elapsed times depend on the host.

Results and restart checkpoints are written under `results/timing_current/<stage>`. Use `REVISION_TIMING_OUT` for an alternate output path and `REVISION_REALDATA_DIR` for prepared real-data inputs in a different directory.
