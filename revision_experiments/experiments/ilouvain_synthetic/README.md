# I-Louvain on the original synthetic generator

This probe covers the 480 balanced and 360 unequal-size inputs described in the reply to Reviewer 1 and Supplementary Section S9. It regenerates graphs and attributes from the two checked-in simulation manifests and evaluates I-Louvain once per input.

Provide an external `IL.py` as described in the [DANCer comparison instructions](../ilouvain_current/). The historical sweep used the Python 3 port with the optimized initial Gini calculation. The supplied implementation remains external, and its SHA-256 and optimization flag are recorded. The compatible upstream constructor uses the full initial Gini calculation, which can take substantially more memory for the larger graphs.

Run from the `revision_experiments` directory.

```sh
Rscript experiments/ilouvain_synthetic/run.R /path/to/IL.py results/ilouvain_synthetic
```

Set `PYTHON` if your Python interpreter is outside the default `python3` path. Append case IDs from either simulation manifest to run a subset, such as `S01_rep01 S49_rep01 U5_A_r1_s3_rep01`. Use a separate output directory for subsets. The output includes ARI and the number of communities returned for every evaluated input.
