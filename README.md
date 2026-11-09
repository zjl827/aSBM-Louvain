# aSBM-Louvain

Code for the experiments in *aSBM-Louvain: Attributed Network Clustering with Louvain-Based Block Structure Learning*. The method combines an attributed stochastic block model with Louvain community detection.

## Reproducing the experiments

Run the commands in the experiment guide from `revision_experiments/`. The guide covers dependencies, input data, and the outputs of each experiment. [Setup and run instructions](revision_experiments/README.md)

- The balanced simulations cover scenarios A1–A3 and B1–B3, with scripts for generating data, running the methods, and plotting the results. [Balanced simulations](revision_experiments/README.md#balanced-simulations)
- The supplementary experiments are indexed by section number. [Supplementary experiments](revision_experiments/README.md#supplementary-experiments)
- The STARmap and MERFISH instructions show how to prepare the input data and run the accuracy and timing comparisons. [Real-data experiments](revision_experiments/experiments/realdata_current/README.md) · [Timing experiments](revision_experiments/experiments/timing_current/README.md)
- The BACT and SpaGCN scripts use their own input and dependency setup. [Baseline scripts](realdata_experiments/) · [Baseline setup](https://github.com/zjl827/aSBM-Louvain/blob/bf6b721df90ad321efc8e1026d26ab27a6a60e57/README.md#real-data-experiments)

The core R entry point is `aSBM_Louvain_revised()`. [Source files and loading order](revision_experiments/README.md#implementation)
