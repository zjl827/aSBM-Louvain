#!/usr/bin/env python3
"""Prepare the exact public DANCer graphs used in Liu et al. (2020)."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


GRAPH_SCENARIOS = {
    "G_strong": "toy_0.graph",
    "G_medium": "toy_2.graph",
    "G_weak": "toy_3.graph",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_dancer_graph(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    markers = [index for index, line in enumerate(lines) if line == "#"]
    if len(markers) < 2:
        raise ValueError(f"Malformed DANCer graph: {path}")

    vertex_lines = lines[markers[0] + 2 : markers[1]]
    edge_lines = lines[markers[1] + 2 :]
    attributes = np.zeros((len(vertex_lines), 2), dtype=float)
    labels = np.zeros(len(vertex_lines), dtype=int)
    for line in vertex_lines:
        node_text, x_text, y_text, label_text = line.split(";")
        node = int(node_text)
        attributes[node] = (float(x_text), float(y_text))
        labels[node] = int(label_text)

    undirected_edges: set[tuple[int, int]] = set()
    for line in edge_lines:
        source_text, target_text = line.split(";")[:2]
        source, target = int(source_text), int(target_text)
        if source == target:
            continue
        undirected_edges.add((min(source, target), max(source, target)))
    edges = np.asarray(sorted(undirected_edges), dtype=int)
    return edges, attributes, labels


def redistributed_attributes(labels: np.ndarray, covariance_scale: float) -> np.ndarray:
    # This reproduces utils.redistribute_attribute(..., random_state=19910420):
    # the seed is reset for each graph and attribute-quality condition.
    np.random.seed(19910420)
    community_count = len(set(labels.tolist()))
    centers = np.random.multivariate_normal(
        mean=np.zeros(2), cov=5 * np.identity(2), size=community_count
    )
    attributes = np.zeros((labels.size, 2), dtype=float)
    covariance = covariance_scale * np.identity(2)
    for node, community in enumerate(labels):
        attributes[node] = np.random.multivariate_normal(
            mean=centers[community], cov=covariance
        )
    return attributes


def write_rows(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for scenario, filename in GRAPH_SCENARIOS.items():
        edges, original_attributes, labels = read_dancer_graph(args.source_dir / filename)
        write_rows(
            args.output_dir / f"{scenario}_edges.csv",
            ["source", "target"],
            [[int(source) + 1, int(target) + 1] for source, target in edges],
        )
        write_rows(
            args.output_dir / f"{scenario}_labels.csv",
            ["node", "label"],
            [[node + 1, int(label) + 1] for node, label in enumerate(labels)],
        )

        attribute_conditions = {
            "strong_attributes": redistributed_attributes(labels, 0.5),
            "medium_attributes": redistributed_attributes(labels, 1.5),
            "weak_attributes": original_attributes,
        }
        for condition, attributes in attribute_conditions.items():
            write_rows(
                args.output_dir / f"{scenario}_{condition}.csv",
                ["node", "attribute_1", "attribute_2"],
                [
                    [node + 1, float(row[0]), float(row[1])]
                    for node, row in enumerate(attributes)
                ],
            )

        print(
            f"Prepared {scenario}: N={labels.size}, M={edges.shape[0]}, "
            f"K={len(set(labels.tolist()))}"
        )


if __name__ == "__main__":
    main()
