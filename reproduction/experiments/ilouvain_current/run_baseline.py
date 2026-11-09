#!/usr/bin/env python3
"""Run the mechanically ported released I-Louvain code on DANCer CSV inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import types
import math
import time
from pathlib import Path

import numpy as np


GRAPHS = ("G_strong", "G_medium", "G_weak")
CONDITIONS = ("strong_attributes", "medium_attributes", "weak_attributes")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("data/dancer"))
    parser.add_argument("--ilouvain-code", type=Path, required=True, help="Path to an external IL.py")
    parser.add_argument("--output", type=Path, default=Path("results/ilouvain_current/ilouvain_results.csv"))
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--conditions", nargs="+", choices=CONDITIONS, default=CONDITIONS)
    return parser.parse_args()


def load_column(path: Path, name: str) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as stream:
        return np.asarray([int(row[name]) for row in csv.DictReader(stream)], dtype=int)


def load_edges(path: Path) -> list[tuple[int, int]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [
            (int(row["source"]) - 1, int(row["target"]) - 1)
            for row in csv.DictReader(stream)
        ]


def load_attributes(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as stream:
        return np.asarray(
            [
                (float(row["attribute_1"]), float(row["attribute_2"]))
                for row in csv.DictReader(stream)
            ],
            dtype=float,
        )


def contingency(truth: np.ndarray, predicted: np.ndarray) -> np.ndarray:
    _, truth_index = np.unique(truth, return_inverse=True)
    _, predicted_index = np.unique(predicted, return_inverse=True)
    table = np.zeros((truth_index.max() + 1, predicted_index.max() + 1), dtype=np.int64)
    np.add.at(table, (truth_index, predicted_index), 1)
    return table


def nmi_geometric(truth: np.ndarray, predicted: np.ndarray) -> float:
    table = contingency(truth, predicted).astype(float)
    joint = table / table.sum()
    row_probability = joint.sum(axis=1)
    column_probability = joint.sum(axis=0)
    expected = np.outer(row_probability, column_probability)
    nonzero = joint > 0
    mutual_information = np.sum(joint[nonzero] * np.log(joint[nonzero] / expected[nonzero]))
    truth_entropy = -np.sum(row_probability[row_probability > 0] * np.log(row_probability[row_probability > 0]))
    predicted_entropy = -np.sum(column_probability[column_probability > 0] * np.log(column_probability[column_probability > 0]))
    if truth_entropy == 0 or predicted_entropy == 0:
        return 0.0
    return float(mutual_information / math.sqrt(truth_entropy * predicted_entropy))


def ari(truth: np.ndarray, predicted: np.ndarray) -> float:
    table = contingency(truth, predicted)
    choose_two = lambda values: np.sum(values * (values - 1) // 2)
    pair_count = truth.size * (truth.size - 1) // 2
    observed = choose_two(table)
    truth_pairs = choose_two(table.sum(axis=1))
    predicted_pairs = choose_two(table.sum(axis=0))
    expected = truth_pairs * predicted_pairs / pair_count
    maximum = 0.5 * (truth_pairs + predicted_pairs)
    return float((observed - expected) / (maximum - expected))


def load_ilouvain(path):
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    source = raw.decode("utf-8")
    # Apply only Python 3 and NetworkX API compatibility to the verified
    # external upstream version. The supplied implementation stays external.
    if digest == "7fed5c6dc562268d7bec27707061b18743b421651953fb35998b1379fd4ff233":
        source = source.replace(".iteritems()", ".items()").replace(".edges_iter(", ".edges(")
        source = source.replace('print "IS THAT POSSIBLE ???? (best_com not in neigh_communities)"',
                                'print("IS THAT POSSIBLE ???? (best_com not in neigh_communities)")')
    module = types.ModuleType("external_ilouvain")
    exec(compile(source, str(path), "exec"), module.__dict__)
    return module.ML2, digest


def main() -> None:
    args = parse_args()
    import networkx as nx

    ML2, source_sha256 = load_ilouvain(args.ilouvain_code)
    vectorized_port = "optimized_initial_gini" in inspect.signature(ML2).parameters
    args.output.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for graph_name in (GRAPHS[:1] if args.smoke else GRAPHS):
        truth = load_column(args.data_dir / f"{graph_name}_labels.csv", "label")
        graph = nx.Graph()
        graph.add_nodes_from(range(truth.size))
        graph.add_edges_from(load_edges(args.data_dir / f"{graph_name}_edges.csv"))

        for condition in (args.conditions[:1] if args.smoke else args.conditions):
            attribute_matrix = load_attributes(args.data_dir / f"{graph_name}_{condition}.csv")
            attributes = {node: attribute_matrix[node] for node in range(truth.size)}
            if not vectorized_port:
                attributes = {node: dict(enumerate(row)) for node, row in attributes.items()}
            author_index = {node: node for node in range(truth.size)}
            start = time.perf_counter()
            kwargs = {"optimized_initial_gini": True} if vectorized_port else {}
            partition = ML2(graph, attributes, author_index, **kwargs).findPartition()
            elapsed = time.perf_counter() - start
            predicted = np.asarray([partition[node] for node in range(truth.size)], dtype=int)
            row = {
                "graph": graph_name,
                "attribute_condition": condition,
                "N": truth.size,
                "M": graph.number_of_edges(),
                "K_true": np.unique(truth).size,
                "K_found": np.unique(predicted).size,
                "ARI": ari(truth, predicted),
                "NMI": nmi_geometric(truth, predicted),
                "runtime_seconds": elapsed,
                "external_code_sha256": source_sha256,
                "vectorized_port": vectorized_port,
            }
            rows.append(row)
            print(
                f"[{graph_name}/{condition}] ARI={row['ARI']:.4f} "
                f"NMI={row['NMI']:.4f} K={row['K_found']} time={elapsed:.2f}s",
                flush=True,
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
