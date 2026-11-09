# Partition-prior feedback

This runner reproduces the paired comparison in Supplementary Section S7. Both branches start from the same labels and post-initialization random state. The Uniform branch replaces the partition-dependent prior and keeps the entropy term. The 480 balanced cases use the same design as `main_simulations`. The fixed STARmap and MERFISH inputs contribute ten matched random runs each.

Run from the `revision_experiments` directory.

```sh
Rscript experiments/prior_feedback/run.R results/prior_feedback balanced
```

After preparing the fixed datasets with `experiments/realdata_current/prepare.R`, use `real` for the 20 fixed-data pairs or `all` for all 500 pairs. Set `REVISION_REALDATA_DIR` if the prepared files are outside `results/realdata_current/inputs`. You can append case IDs such as `S01_rep01` or `STARmap__expression_graph_spatial_attributes__rep01` to select a subset. Use a separate output directory for subsets. The output gives paired scores and their mean and sample SD by condition.
