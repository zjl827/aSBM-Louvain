from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


def _to_dense(x):
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def filter_cells(x, min_genes=None, max_genes=None, min_counts=None, max_counts=None):
    arr = _to_dense(x)
    gene_counts = (arr > 0).sum(axis=1)
    total_counts = arr.sum(axis=1)
    mask = np.ones(arr.shape[0], dtype=bool)
    if min_genes is not None:
        mask &= gene_counts >= min_genes
    if max_genes is not None:
        mask &= gene_counts <= max_genes
    if min_counts is not None:
        mask &= total_counts >= min_counts
    if max_counts is not None:
        mask &= total_counts <= max_counts
    return mask, total_counts


def filter_genes(x, min_cells=None, max_cells=None, min_counts=None, max_counts=None):
    arr = _to_dense(x)
    cell_counts = (arr > 0).sum(axis=0)
    total_counts = arr.sum(axis=0)
    mask = np.ones(arr.shape[1], dtype=bool)
    if min_cells is not None:
        mask &= cell_counts >= min_cells
    if max_cells is not None:
        mask &= cell_counts <= max_cells
    if min_counts is not None:
        mask &= total_counts >= min_counts
    if max_counts is not None:
        mask &= total_counts <= max_counts
    return mask, total_counts


def _normalize_inplace(adata, target_sum):
    x = adata.X
    if sparse.issparse(x):
        totals = np.asarray(x.sum(axis=1)).ravel()
        scale = np.divide(target_sum, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        adata.X = sparse.diags(scale).dot(x)
    else:
        totals = x.sum(axis=1)
        scale = np.divide(target_sum, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        adata.X = x * scale[:, None]
    adata.obs["size_factors"] = scale


def normalize_total(adata, target_sum=1e4):
    _normalize_inplace(adata, target_sum)


def normalize_per_cell(adata, counts_per_cell_after=None):
    if counts_per_cell_after is None:
        x = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        totals = x.sum(axis=1)
        positive = totals[totals > 0]
        counts_per_cell_after = float(np.median(positive)) if positive.size else 1.0
    _normalize_inplace(adata, counts_per_cell_after)


def log1p(adata, copy=False):
    target = adata.copy() if copy else adata
    if sparse.issparse(target.X):
        target.X = target.X.copy()
        target.X.data = np.log1p(target.X.data)
    else:
        target.X = np.log1p(np.asarray(target.X))
    return target if copy else None


def neighbors(adata, n_neighbors=10):
    x = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    if x.shape[0] <= 1:
        conn = sparse.csr_matrix((x.shape[0], x.shape[0]))
    else:
        k = min(max(1, n_neighbors), x.shape[0] - 1)
        nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")
        nbrs.fit(x)
        distances, indices = nbrs.kneighbors(x)
        rows = np.repeat(np.arange(x.shape[0]), k)
        cols = indices[:, 1 : k + 1].reshape(-1)
        vals = 1.0 / (1.0 + distances[:, 1 : k + 1].reshape(-1))
        conn = sparse.csr_matrix((vals, (rows, cols)), shape=(x.shape[0], x.shape[0]))
        conn = conn.maximum(conn.transpose())

    adata.obsp["connectivities"] = conn
    adata.uns["neighbors"] = {"connectivities_key": "connectivities", "n_neighbors": n_neighbors}
