# Paper experiments

## Setup

Run the commands from the `reproduction/` directory. The experiments used R 4.5.2, igraph 2.2.2, Matrix 1.7-4, mvtnorm 1.3-5, and mclust 6.1.2. Matching R and igraph versions helps reproduce seeded Louvain partitions.

```r
install.packages(c("igraph", "Matrix", "mvtnorm", "mclust"))
```

```sh
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Use separate output directories for partial runs and full runs. Use a fresh directory when the code or inputs differ from saved checkpoints. For runners that use forked R workers, set the worker count to one on Windows.

## Balanced simulations

The [task grid](experiments/main_simulations/tasks.csv) contains 48 settings with ten independent data replicates each. The two attributed fits share an initialization. The standalone K-means reference uses the planted community count.

```sh
Rscript experiments/main_simulations/run.R
python3 experiments/main_simulations/plot_results.py
```

Append case IDs such as `S01_rep01` to run a subset and set `MAIN_SIM_OUTPUT` to a separate directory. The runner defaults to `outputs/main_simulations` and resumes completed cases. To plot a full run from another location, set `MAIN_SIM_OUTPUT` for the runner and pass its `raw_metrics.csv` to the plotter with `--metrics`.

## Supplementary experiments

```sh
# S2. Unequal community sizes
Rscript experiments/unequal_sizes/run.R results/unequal_sizes

# S3. Joint signal strength and community-dependent noise
Rscript experiments/joint_signals/run.R
python3 experiments/joint_signals/analyze.py

# S4 and S8. Non-Gaussian attributes and attribute likelihood weights
Rscript experiments/weight_misspecification/run.R
python3 experiments/weight_misspecification/analyze.py

# S5. Attribute dimension. The extension reads the pilot cache and results.
Rscript experiments/dimension_pilot/run.R
Rscript experiments/dimension_pilot/analyze.R
Rscript experiments/dimension_extension/run.R
Rscript experiments/dimension_extension/analyze.R

# S6. Initialization stability
python3 experiments/initialization_stability/prepare.py
INIT_STAGE=full Rscript experiments/initialization_stability/run.R
python3 experiments/initialization_stability/analyze.py

# S7. Prior feedback. Use balanced for simulations, real or all after data setup below.
Rscript experiments/prior_feedback/run.R results/prior_feedback balanced
```

### S9. DANCer and I-Louvain

Use the [benchmark repository at commit 515bf25](https://github.com/changliu94/attributed-community-detection/tree/515bf251458be543cee43f27eb234bcef544ac4f). Keep its graphs and I-Louvain implementation in a separate checkout.

```sh
python3 experiments/ilouvain_current/prepare.py --source-dir /path/to/benchmark/dataset --output-dir data/dancer
Rscript experiments/ilouvain_current/run_proposed.R
python3 experiments/ilouvain_current/run_baseline.py --ilouvain-code /path/to/benchmark/src/IL.py
Rscript experiments/ilouvain_current/summarize.R
```

### S10–S12. Real-data experiments and graph bootstrap

Prepare the [STARmap and MERFISH inputs](experiments/realdata_current/README.md) before running the real-data part of S7 or any of S10–S12. S11 also uses the ordered neighbor lists from that preparation.

```sh
# S10. Graph construction, real-data accuracy, and matched runtime
Rscript experiments/realdata_current/run.R
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript experiments/timing_current/run.R

# S11. Resolution and parameter sensitivity
Rscript experiments/parameter_sensitivity/run.R results/parameter_sensitivity

# S12. Conditional graph bootstrap
Rscript experiments/graph_bootstrap/run.R
python3 experiments/graph_bootstrap/summarize.py
```

The timing and bootstrap runners accept `smoke` without real-data inputs. For bootstrap smoke summaries, pass `--input-dir results/graph_bootstrap/smoke --expected-sources 1 --expected-b 1` to `summarize.py`.

### Additional checks

```sh
Rscript experiments/ilouvain_synthetic/run.R /path/to/benchmark/src/IL.py results/ilouvain_synthetic
Rscript experiments/potential_diagnostics/run.R results/potential_diagnostics
```
