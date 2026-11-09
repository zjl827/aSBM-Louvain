"""Summarize and plot ARI for the six balanced simulation scenarios."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


# Same methods, colors, and line styles as the original repository figures.
STYLES = {
    "Proposed": ("aSBM-Louvain", "#c65d08", "o", "-"),
    "KMeans_knownK": ("K-means", "#d99600", "^", (0, (5, 5))),
    "Louvain": ("Louvain", "#8f8f8f", "x", (0, (1.5, 4))),
    "Pure_aSBM": ("aSBM", "#5dade2", "s", (0, (10, 5))),
}
SCENARIOS = ("A1", "A2", "A3", "B1", "B2", "B3")


def draw_panel(ax, summary, scenario):
    is_a = scenario.startswith("A")
    x = "r" if is_a else "sigma"
    panel = summary.loc[summary.scenario.eq(scenario)].sort_values(x)
    ax.set_axisbelow(True)
    ax.grid(color="#e6e6e6")
    for method, (label, color, marker, linestyle) in STYLES.items():
        data = panel.loc[panel.method.eq(method)]
        ax.plot(data[x], data.ARI_mean, color=color, marker=marker,
                linestyle=linestyle, markerfacecolor="white", label=label)
    for method in ("KMeans_knownK", "Louvain", "Pure_aSBM", "Proposed"):
        _, color, _, _ = STYLES[method]
        data = panel.loc[panel.method.eq(method)]
        ax.errorbar(data[x], data.ARI_mean, yerr=data.ARI_SD, fmt="none",
                    ecolor=color, capsize=3, zorder=4)
    ax.set(xlabel="Structure strength (r)" if is_a else "Attribute noise (σ)",
           ylabel="ARI", ylim=(-0.05, 1.05), yticks=[0, 0.25, 0.5, 0.75, 1])
    if is_a:
        ax.invert_xaxis()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path,
                        default=Path("outputs/main_simulations/raw_metrics.csv"))
    parser.add_argument("--out", type=Path, default=Path("outputs/main_simulations/figures"))
    args = parser.parse_args()
    raw = pd.read_csv(args.metrics)
    raw = raw.loc[raw.scope.eq("balanced")]
    counts = raw.groupby(["point_id", "method"]).ARI.count()
    if (len(counts) != 192 or not counts.eq(10).all()
            or raw.duplicated(["point_id", "method", "rep"]).any()):
        parser.error("Expected 10 finite ARIs per scenario point and method.")
    summary = (
        raw.groupby(["point_id", "scenario", "N", "K_true", "r", "sigma", "method"])
        .ARI.agg(ARI_mean="mean", ARI_SD="std", n="count").reset_index()
    )
    args.out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out / "condition_method_summary.csv", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    for ax, scenario in zip(axes.flat, SCENARIOS):
        draw_panel(ax, summary, scenario)
        ax.set_title(scenario)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncols=4, frameon=False)
    fig.savefig(args.out / "balanced_simulations.pdf")
    plt.close(fig)
    for scenario in SCENARIOS:
        fig, ax = plt.subplots(figsize=(8, 4.8), layout="constrained")
        draw_panel(ax, summary, scenario)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncols=4, frameon=False)
        fig.savefig(args.out / f"scenario_{scenario}.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()
