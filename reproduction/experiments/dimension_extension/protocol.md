# Extension to 400 and 800 attribute dimensions

This part of Supplementary Section S5 extends the [50-, 100-, and 200-dimensional experiment](../dimension_pilot/protocol.md) to 400 and 800 dimensions at noise SD 1.5 and 4. It uses the same ten independent replicates and pairs all five dimensions. The extension adds 160 method records across four conditions and ten replicates.

The graph, 50-dimensional community means, and seeds come from the original replicate. Using the same residual seed, the runner generates a $1{,}000\times800$ matrix of standard normal residuals and checks that the first 200 columns match the earlier input exactly. Population community means are zero in coordinates 51 through 800. Methods, stopping settings, seeds, and initialization follow the first part. The extension also logs the aSBM objective values and records whether its returned partition equals the K-means initialization.

Run the first part from the same checkout before running `run.R` and `analyze.R` in this directory. The extension uses the first part's input cache and raw metrics via `DIM_PILOT_OUTPUT` (default `outputs/dimension_pilot`). Its output defaults to `outputs/dimension_extension`. Set `DIM_OUTPUT` to use another directory. The combined analysis writes the five-dimension summary and paired ARI differences relative to dimension 50.

`DIM_WORKERS` sets the fitting worker count (default 2). With `DIM_SMOKE=1`, use a matching one-replicate pilot cache and separate output directories. The validation and combined analysis expect the full grid and ten replicates.
