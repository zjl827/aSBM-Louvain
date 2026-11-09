# Finite potential diagnostics

This runner reproduces the finite trajectory check reported in the reply to Reviewer 1, Major Comment 1. It uses the 480 seeded inputs from the balanced main study and adds a read-only observer to the shared wrapper at runtime. Each traced fit is checked against an ordinary fit from the same initial labels and random state.

Run from the `revision_experiments` directory.

```sh
Rscript experiments/potential_diagnostics/run.R results/potential_diagnostics
```

Append case IDs from `experiments/main_simulations/tasks.csv` to run a subset. Use a separate output directory for subsets. The output records the four increments for eligible fixed-dimension transitions, with the diagnostic coefficient fixed at one. It counts decreases using the absolute threshold of (10^{-8}). These finite traces describe observed trajectories and cannot establish summability along an infinite continuation.
