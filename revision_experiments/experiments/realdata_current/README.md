# S10 real-data graph and attribute orientation

From the release root, install the R packages `igraph`, `Matrix`, `mclust` and `mvtnorm`. Input preparation also needs Python with `h5py`, `numpy`, `pandas` and `scikit-learn`. Set `REVISION_PYTHON` if that environment uses a command other than `python3`.

The two processed files already in the parent project checkout are `starmap_PCA_scores_50PCs.csv` and `coord_and_pc_19.RData`. The default processed input directory is `..`. Use `--processed-dir` to point to another copy.

Place these separately obtained files in `data/realdata/` or pass their directory as `--input-dir`.

- [STARmap H5AD](http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/STARmap_20180505_BY3_1k.h5ad), saved as `STARmap_20180505_BY3_1k.h5ad`.
- [STARmap cell-type annotations](https://drive.google.com/drive/folders/1I1nxheWlc2RXSdiv24dex3YRaEh78Omy?usp=sharing), saved as `Annotation_STARmap_20180505_BY3_1k.txt`. This tab-delimited file supplies 16 cell types via its `Annotation` column and the aligned X/Y coordinates. The H5AD supplies the coordinates for the spatial graph. The cell-type scores use the annotation file.
- [MERFISH H5AD](http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/MERFISH_0.19.h5ad), saved as `MERFISH_0.19.h5ad`. The 5,803 cells are aligned to the processed PCs by ID. The reference uses the H5AD's `cell_class` field (15 types).

```sh
Rscript experiments/realdata_current/prepare.R
Rscript experiments/realdata_current/run.R
```

Preparation writes four input RDS files and two nested expression kNN lists under `results/realdata_current/inputs`. They are also consumed by the real-data S7, S10 timing, S11 and S12 runners. Set `REVISION_REALDATA_DIR` or use `--output-dir` when preparing an alternate location. Each RDS contains `adjacency`, `attributes`, `truth`, `cell_ids`, and `resolution` along with the input definitions. The preparation uses symmetrized Euclidean 6NN graphs in unscaled PCA or spatial coordinates and column-standardized attributes in the opposite modality. The top-12 CSVs use the same PCA cell order, exclude self and index neighbors from one.

The full run fits four methods with ten algorithmic repeats on each of the four fixed inputs. Proposed and Pure-aSBM share estimated-K attribute K-means initialization (5 starts, 10 iterations), seeds `20260901 + rep`, and the post-initialization RNG stream. The ordinary Louvain baseline uses the original sorted-cell order, resolution 1 and seeds `123 + rep`. Proposed uses resolution 1.5 for STARmap and 0.5 for MERFISH. Pure-aSBM uses the vectorized R implementation. The `smoke` positional option runs the first STARmap primary-orientation repeat. Results and restart checkpoints go in `results/realdata_current` or `REVISION_REALDATA_OUT`. Set `REVISION_REALDATA_OUT` to a fresh directory for a smoke run or when the code or inputs differ from saved checkpoints.
