# S9 DANCer comparison with I-Louvain

From the release root, install R packages `igraph`, `mclust` and `mvtnorm` and Python packages `numpy`, `scipy` and `networkx`. Obtain the [benchmark repository](https://github.com/changliu94/attributed-community-detection) at commit `515bf251458be543cee43f27eb234bcef544ac4f` separately. Use its `dataset/` directory containing `toy_0.graph`, `toy_2.graph` and `toy_3.graph`, and the external `src/IL.py`.

```sh
python3 experiments/ilouvain_current/prepare.py --source-dir /path/to/benchmark/dataset --output-dir data/dancer
Rscript experiments/ilouvain_current/run_proposed.R
python3 experiments/ilouvain_current/run_baseline.py --ilouvain-code /path/to/benchmark/src/IL.py
Rscript experiments/ilouvain_current/summarize.R
```

The preparation recreates the 3 graph × 3 attribute combinations, using the original DANCer attributes for the weak condition. The pinned external `IL.py` is adapted at load time for Python 3 and NetworkX APIs. Its nine ARI/NMI results match the validated vectorized port used in the paper. The implementation and input graphs stay in the external checkout. The baseline runner records the supplied source's SHA256.

Each fixed condition has ten Proposed repeats and one deterministic I-Louvain fit. The Proposed seeds and all settings match the paper, including the 5-start, 10-iteration K-means initializer. Runs returning at the outer cap remain in the table. Use `smoke` with the R runner and `--smoke` with the Python runner to fit the first setting only. The default data path is `data/dancer` and the result path is `results/ilouvain_current`. Override them using `REVISION_DANCER_DIR` and `REVISION_DANCER_OUT` for R or `--data-dir` and `--output` for Python. Use separate smoke and formal output directories.
