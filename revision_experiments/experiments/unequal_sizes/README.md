# Unequal community sizes

This runner reproduces the four methods plotted in Supplementary Section S2 on 360 seeded synthetic inputs. The task grid covers both size ratios, both one-factor sweeps, and the moderate-signal stress conditions. The proposed method and Pure-aSBM share attribute K-means labels and the random state after initialization. The standalone K-means reference uses the planted community count.

Run the full grid from the `revision_experiments` directory.

```sh
Rscript experiments/unequal_sizes/run.R results/unequal_sizes
```

Append case IDs from `tasks.csv` to run only those cases. For example, use `S49_rep01 U3_C_r3.5_s3_rep01 U5_A_r1_s3_rep01` to select three cases across the size ratios and the stress grid. Use a separate output directory for subsets. The output includes per-case scores, a combined metrics file, and mean ARI with sample SD by condition and method.
