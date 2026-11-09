"""Prepare H5AD-aligned metadata and historical Euclidean neighbor lists."""

from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors


def decode(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def read_labels(values):
    if isinstance(values, h5py.Dataset):
        return [decode(value) for value in values[:]]
    categories = [decode(value) for value in values["categories"][:]]
    codes = values["codes"][:]
    if np.any(codes < 0):
        raise ValueError("Missing reference labels in the H5AD")
    return [categories[code] for code in codes]


def read_metadata(path, label_key):
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        cells = [decode(value) for value in obs[decode(obs.attrs["_index"])][:]]
        coords = handle["obsm"]["spatial"][:]
        metadata = pd.DataFrame(coords, index=pd.Index(cells, name="Cell"), columns=["X", "Y"])
        metadata["truth"] = read_labels(obs[label_key])
    if not metadata.index.is_unique:
        raise ValueError(f"Duplicate H5AD cell IDs in {path}")
    return metadata


def neighbors(features, k, query_self):
    matrix = features.to_numpy(dtype=float)
    model = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", metric="euclidean").fit(matrix)
    candidates = model.kneighbors(matrix if query_self else None, return_distance=False)
    rows = []
    for i, candidates_i in enumerate(candidates):
        chosen = candidates_i[candidates_i != i][:k]
        if len(chosen) != k:
            raise ValueError(f"Too few non-self neighbors for row {i + 1}")
        rows.append(chosen + 1)
    return pd.DataFrame(rows, index=pd.Index(features.index, name="Cell"),
                        columns=[f"Neighbor{i}" for i in range(1, k + 1)])


def main():
    input_dir, pca_star, pca_mer, tmp, output = map(Path, sys.argv[1:])
    for stem, filename, label_key, pca_file in (
        ("starmap", "STARmap_20180505_BY3_1k.h5ad", "ground_truth", pca_star),
        ("merfish", "MERFISH_0.19.h5ad", "cell_class", pca_mer),
    ):
        pcs = pd.read_csv(pca_file, index_col=0)
        pcs.index = pcs.index.astype(str)
        metadata = read_metadata(input_dir / filename, label_key).loc[pcs.index]
        # The original spatial kNN script read its H5AD metadata back as CSV.
        metadata.to_csv(tmp / f"{stem}_metadata.csv")
        coords = pd.read_csv(tmp / f"{stem}_metadata.csv", index_col=0)[["X", "Y"]]
        expr6 = neighbors(pcs, 6, query_self=True)
        expr12 = neighbors(pcs, 12, query_self=False)
        if not expr12.iloc[:, :6].equals(expr6):
            raise ValueError(f"{stem}: nested expression neighbor list does not match the 6NN list")
        expr6.to_csv(tmp / f"{stem}_expression_knn6.csv")
        neighbors(coords, 6, query_self=True).to_csv(tmp / f"{stem}_spatial_knn6.csv")
        expr12.to_csv(output / f"{stem}_expression_knn12.csv")


if __name__ == "__main__":
    main()
