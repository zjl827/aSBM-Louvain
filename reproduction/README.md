# aSBM-Louvain experiment guide

This directory contains the scripts for the balanced simulations and the experiments in the Supplementary Material. Run the commands below from this directory.

## Setup

The experiments used R 4.5.2 with igraph 2.2.2, Matrix 1.7-4, mvtnorm 1.3-5, and mclust 6.1.2. Matching R and igraph versions matters because several designs use seeded Louvain partitions. Install the R dependencies and the Python packages used by the analysis scripts.

```r
install.packages(c("igraph", "Matrix", "mvtnorm", "mclust"))
```

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

## Implementation

The experiment scripts load the repository's `asbm_louvain.R` and call `aSBM_Louvain_revised()`. The standard configuration uses attribute K-means initialization with a Louvain estimate of the number of groups, prior feedback, and the joint stopping rule. Each runner sets the parameters specific to its experiment. [Algorithm source](../asbm_louvain.R)

The aSBM baseline in `code/pure_asbm_vectorized.R` preserves the per-pair denominator floor using the vectorized R implementation. For a fresh run, use an output directory without checkpoints from another version.

## Balanced simulations

`experiments/main_simulations/tasks.csv` records 48 settings across A1–A3 and B1–B3, each with 10 data replicates. Each data replicate regenerates the network and attributes, including community means. The two attributed models start from the same partition and post-initialization random-number state. The standalone K-means reference uses the planted number of groups.

```bash
# Quick check of one case, with four ARI values
Rscript experiments/main_simulations/run.R S01_rep01

# Full design and figures
Rscript experiments/main_simulations/run.R
python3 experiments/main_simulations/plot_results.py
```

The accuracy-only aSBM fits can take substantially longer than the other methods. The runner saves one checkpoint per case under `outputs/main_simulations/runs/` and resumes from completed checkpoints. Set `MAIN_SIM_OUTPUT` to use a separate output directory and pass the resulting CSV to the plotting script with `--metrics`. The analysis writes a table and the six ARI panels using means and sample standard deviations over the data replicates.

## Supplementary experiments

| Section | Experiment | Directory |
| --- | --- | --- |
| S2 | Unequal community sizes | [unequal_sizes](experiments/unequal_sizes/) |
| S3 | Joint signal strength and community-dependent noise | [joint_signals](experiments/joint_signals/) |
| S4, S8 | Non-Gaussian noise and attribute likelihood weighting | [weight_misspecification](experiments/weight_misspecification/) |
| S5 | Attribute dimension | [dimension_pilot](experiments/dimension_pilot/), [dimension_extension](experiments/dimension_extension/) |
| S6 | Initialization stability | [initialization_stability](experiments/initialization_stability/) |
| S7 | Prior feedback | [prior_feedback](experiments/prior_feedback/) |
| S9 | DANCer comparison with I-Louvain | [ilouvain_current](experiments/ilouvain_current/) |
| S10 | Graph construction, real-data accuracy, and runtime | [realdata_current](experiments/realdata_current/), [timing_current](experiments/timing_current/) |
| S11 | Resolution and parameter sensitivity | [parameter_sensitivity](experiments/parameter_sensitivity/) |
| S12 | Conditional graph bootstrap | [graph_bootstrap](experiments/graph_bootstrap/) |

The joint-signal, noise, dimension, and initialization experiments generate their inputs from the included designs and seeds. The other directories have their own setup and run commands. Start with the real-data preparation before running the real-data branches of S7, S10, S11, and S12.

```bash
Rscript experiments/joint_signals/run.R
python3 experiments/joint_signals/analyze.py

# Run once for both S4 and S8.
Rscript experiments/weight_misspecification/run.R
python3 experiments/weight_misspecification/analyze.py

# The dimension extension uses the pilot inputs and results.
Rscript experiments/dimension_pilot/run.R
Rscript experiments/dimension_pilot/validate.R
Rscript experiments/dimension_pilot/analyze.R
Rscript experiments/dimension_extension/run.R
Rscript experiments/dimension_extension/analyze.R

python3 experiments/initialization_stability/prepare.py
INIT_STAGE=full Rscript experiments/initialization_stability/run.R
python3 experiments/initialization_stability/analyze.py
```

For a quick run, use a separate output directory and set `PILOT_SMOKE=1`, `WM_SMOKE=1`, `DIM_SMOKE=1`, or `INIT_STAGE=smoke` for the corresponding runner. The full-design analysis scripts expect all replicates. Where an R runner uses forked workers, set its worker count to one on Windows.

The [synthetic I-Louvain probe](experiments/ilouvain_synthetic/) and [finite potential diagnostics](experiments/potential_diagnostics/) reproduce additional checks of clustering and potential trajectories.

## Data and outputs

The synthetic experiments generate their data locally. For the fixed datasets, follow the corresponding directory's instructions to prepare the inputs using the existing repository's preprocessed files and the public data sources. Generated outputs, additional data downloads, and compiled libraries are excluded from version control.

An independent data replicate is the unit for the simulation ARI summaries. On a fixed graph and attribute matrix, variation across repeated runs reflects the algorithm's randomness. See the individual experiment protocols for paired comparisons and graph-bootstrap anchors.
