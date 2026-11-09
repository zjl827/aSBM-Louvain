#!/usr/bin/env python3
"""Fit an external I-Louvain implementation on one generated network."""

import argparse
import inspect
import sys
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ilouvain_current"))
from run_baseline import load_ilouvain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ilouvain-code", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--attributes", type=Path, required=True)
    parser.add_argument("--labels-out", type=Path, required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--d", type=int, required=True)
    parser.add_argument("--array-attributes", action="store_true")
    args = parser.parse_args()

    model, digest = load_ilouvain(args.ilouvain_code)
    supports_optimized = "optimized_initial_gini" in inspect.signature(model).parameters
    graph = nx.Graph()
    graph.add_nodes_from(range(args.n))
    edges = np.fromfile(args.edges, dtype="<i4").reshape((-1, 2))
    graph.add_edges_from((int(u), int(v)) for u, v in edges)
    x = np.fromfile(args.attributes, dtype="<f8").reshape((args.n, args.d))
    if args.array_attributes and supports_optimized:
        attrs = dict(enumerate(x))
    else:
        attrs = {node: dict(enumerate(map(float, row))) for node, row in enumerate(x)}
    kwargs = {"optimized_initial_gini": True} if supports_optimized else {}
    index = {node: node for node in range(args.n)}
    partition = model(graph, attrs, index, **kwargs).findPartition()
    np.asarray([partition[node] for node in range(args.n)], dtype="<i4").tofile(args.labels_out)
    print(f"external_code_sha256={digest} optimized_initial_gini={supports_optimized}")


if __name__ == "__main__":
    main()
