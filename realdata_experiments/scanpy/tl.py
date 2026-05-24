from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
import igraph as ig


def louvain(adata, resolution=1.0):
    conn = adata.obsp.get("connectivities")
    if conn is None:
        raise ValueError("Call scanpy.pp.neighbors() before scanpy.tl.louvain().")

    matrix = conn if sparse.issparse(conn) else sparse.csr_matrix(conn)
    matrix = sparse.triu(matrix, k=1).tocsr()
    if matrix.shape[0] == 0:
        labels = []
    elif matrix.nnz == 0:
        labels = [0] * matrix.shape[0]
    else:
        rows, cols = matrix.nonzero()
        weights = matrix.data.astype(float).tolist()
        graph = ig.Graph(n=matrix.shape[0], edges=list(zip(rows.tolist(), cols.tolist())), directed=False)
        clustering = graph.community_multilevel(weights=weights, resolution=max(float(resolution), 1e-6))
        labels = clustering.membership

    adata.obs["louvain"] = pd.Categorical(np.asarray(labels, dtype=int).astype(str))


def rank_genes_groups(*args, **kwargs):
    raise NotImplementedError("rank_genes_groups is not implemented in the lightweight benchmark scanpy shim.")
