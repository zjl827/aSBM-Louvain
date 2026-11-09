# Extension to 400 and 800 attribute dimensions

This part of Supplementary Section S5 extends the [50-, 100-, and 200-dimensional experiment](../dimension_pilot/protocol.md) to 400 and 800 dimensions at noise SD 1.5 and 4. It uses the same ten independent replicates and pairs all five dimensions. The extension adds 160 method records across four conditions and ten replicates.

The graph, 50-dimensional community means, and seeds come from the original replicate. Using the same residual seed, the runner generates a $1{,}000\times800$ matrix of standard normal residuals and checks that the first 200 columns match the earlier input exactly. Population community means are zero in coordinates 51 through 800. Methods, stopping settings, seeds, and initialization follow the first part. The extension also logs the aSBM objective values and records whether its returned partition equals the K-means initialization.

Run the first part through its analysis from the same checkout before running `run.R` and `analyze.R` in this directory. The extension reads its input cache and baseline metrics through `DIM_PILOT_OUTPUT` (default `outputs/dimension_pilot`). Its output defaults to `outputs/dimension_extension`; use `DIM_OUTPUT` to override it. The combined analysis writes the five-dimension summary and paired ARI differences relative to dimension 50. `validate.R` is an optional independent check of full runs and nested inputs.

`DIM_WORKERS` sets the fitting worker count (default 2). With `DIM_SMOKE=1`, use a matching one-replicate pilot cache and separate output directories. The validation and combined analysis expect the full grid and ten replicates.
