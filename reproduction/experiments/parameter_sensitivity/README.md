# Resolution and parameter sensitivity

The 340 tasks in Supplementary Section S11 vary one parameter at a time on fixed STARmap and MERFISH data. Each dataset has one shared default configuration and four nondefault resolution settings. The remaining settings vary the expression-graph neighbor count, cooling factor, initial annealing strength, or stopping tolerances. Resolution is applied both during the initial community-count estimate and subsequent Louvain updates.

Prepare the fixed datasets with `experiments/realdata_current/prepare.R`, then run the sweep from the `reproduction` directory.

```sh
Rscript experiments/parameter_sensitivity/run.R results/parameter_sensitivity
```

The input directory defaults to `results/realdata_current/inputs`. It contains the expression-graph and spatial-attribute RDS files and the ordered `starmap_expression_knn12.csv` and `merfish_expression_knn12.csv` files. Set `REVISION_REALDATA_DIR` to use a different directory. If the neighbor CSVs are in a separate directory, set `REVISION_KNN_DIR` as well.

Append case IDs from `tasks.csv` to run a subset. Use a separate output directory for subsets. The output includes per-run scores and counts, and summaries by dataset and configuration.
