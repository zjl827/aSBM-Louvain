# STARmap and MERFISH

Install the [experiment dependencies](../../README.md#setup). The repository root contains `starmap_PCA_scores_50PCs.csv` and `coord_and_pc_19.RData`. Put these additional files in `reproduction/data/realdata` or pass `--input-dir` when running `prepare.R`.

- [STARmap H5AD](http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/STARmap_20180505_BY3_1k.h5ad), saved as `STARmap_20180505_BY3_1k.h5ad`.
- [STARmap cell-type annotations](https://drive.google.com/drive/folders/1I1nxheWlc2RXSdiv24dex3YRaEh780my?usp=sharing), saved as `Annotation_STARmap_20180505_BY3_1k.txt`.
- [MERFISH H5AD](http://sdmbench.drai.cn/tcm/download/?file_path=/mnt/JINGD/data/file/sdmbench/db/MERFISH_0.19.h5ad), saved as `MERFISH_0.19.h5ad`. The reference labels use its `cell_class` field.

Run from `reproduction/`.

```sh
Rscript experiments/realdata_current/prepare.R
Rscript experiments/realdata_current/run.R
```

Preparation writes the four graph/attribute inputs and two expression neighbor lists to `results/realdata_current/inputs`. Set `REVISION_REALDATA_DIR` in the runners if you use `--output-dir` to prepare the inputs elsewhere. The run script accepts `smoke` for one STARmap repeat. Use `REVISION_REALDATA_OUT` for a separate results directory.

## BACT and SpaGCN

The [baseline scripts](../../../realdata_experiments/) read inputs from the repository root. Install [BACT](https://github.com/yinqiaoyan/BACT) for the R scripts. The Python scripts use [SpaGCN](https://github.com/jianhuupenn/SpaGCN), PyTorch, anndata, and the Python analysis dependencies.

| Dataset | BACT inputs in addition to the bundled PCs | SpaGCN inputs |
| --- | --- | --- |
| STARmap | `starmap_spatial.csv`, `Annotation_STARmap_20180505_BY3_1k.txt` | `starmap_spatial.csv`, `starmap_counts.csv`, `Annotation_STARmap_20180505_BY3_1k.txt` |
| MERFISH | `cell_type_annotation_0.19.csv` | `cell_type_annotation_0.19.csv`, `merfish19_spatial.csv`, `merfish19_counts.csv` |

Run from the repository root:

```sh
Rscript realdata_experiments/starmap_bact_pca_coord_input_10runs.R
Rscript realdata_experiments/merfish_bact_pc_coord_input_10runs.R
python3 realdata_experiments/starmap_spagcn_expression_coord_input_10runs.py
python3 realdata_experiments/merfish_spagcn_expression_coord_input_10runs.py
```
