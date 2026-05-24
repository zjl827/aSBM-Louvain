"""
Generate publication figures for aSBM-Louvain simulation experiments.
Reads result CSV files and plots ARI curves for all scenarios.
"""

import matplotlib.pyplot as plt
import numpy as np

# ---- Common style ----
COLORS = {
    "aSBM-Louvain": "#c65d08",
    "KMeans":       "#d99600",
    "Louvain":      "#8f8f8f",
    "Pure-aSBM":    "#5dade2",
}
MARKERS = {
    "aSBM-Louvain": "o",
    "KMeans":       "^",
    "Louvain":      "x",
    "Pure-aSBM":    "s",
}
LINESTYLES = {
    "aSBM-Louvain": (0, ()),
    "KMeans":       (0, (5, 5)),
    "Louvain":      (0, (1.5, 4)),
    "Pure-aSBM":    (0, (10, 5)),
}


def make_figure(x, data, xlabel, filename, xlim=None, xticks=None):
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.grid(True, color="#e6e6e6", linewidth=1)

    order = ["aSBM-Louvain", "KMeans", "Louvain", "Pure-aSBM"]
    for method in order:
        y = data[method]
        ax.plot(x, y, color=COLORS[method], marker=MARKERS[method],
                linestyle=LINESTYLES[method], linewidth=1.6,
                markersize=6 if method in ("aSBM-Louvain", "Pure-aSBM") else 7,
                label=method)

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel("ARI", fontsize=16)
    ax.set_ylim(-0.05, 1.05)
    ax.tick_params(axis='both', labelsize=14)

    if xlim:
        ax.set_xlim(xlim)
    if xticks:
        ax.set_xticks(xticks)

    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4,
              frameon=False, fontsize=15, markerscale=1.8)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved: {filename}")


# ---- Data ----
x_a = np.array([7.5, 6.0, 4.5, 3.0, 1.5, 1.0])
x_b = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])

# A1: N=1000, K=5, sigma=3.0
make_figure(x_a, {
    "Louvain":      [0.997, 0.981, 0.929, 0.071, 0.001, -0.000],
    "KMeans":       [0.606, 0.598, 0.639, 0.632, 0.607, 0.642],
    "Pure-aSBM":    [0.991, 0.989, 0.929, 0.729, 0.751, 0.647],
    "aSBM-Louvain": [0.952, 0.904, 0.901, 0.738, 0.007, -0.001],
}, xlabel="Structure Strength (r)", filename="scenario_A1_structure_degradation.png",
   xlim=(7.8, 0.7), xticks=[6, 4, 2])

# A2: N=1000, K=10, sigma=3.0
make_figure(x_a, {
    "Louvain":      [0.932, 0.657, 0.070, 0.013, 0.002, -0.001],
    "KMeans":       [0.430, 0.512, 0.401, 0.444, 0.442, 0.419],
    "Pure-aSBM":    [0.923, 0.917, 0.893, 0.817, 0.689, 0.650],
    "aSBM-Louvain": [0.889, 0.738, 0.478, 0.084, 0.004, -0.001],
}, xlabel="Structure Strength (r)", filename="scenario_A2_structure_degradation.png",
   xlim=(7.8, 0.7), xticks=[6, 4, 2])

# A3: N=5000, K=5, sigma=3.0
make_figure(x_a, {
    "Louvain":      [0.926, 0.634, 0.453, 0.007, 0.001, 0.000],
    "KMeans":       [0.670, 0.641, 0.663, 0.683, 0.686, 0.659],
    "Pure-aSBM":    [0.984, 0.972, 0.945, None, None, None],
    "aSBM-Louvain": [0.997, 0.993, 0.970, 0.878, 0.005, 0.001],
}, xlabel="Structure Strength (r)", filename="scenario_A3_structure_degradation.png",
   xlim=(7.8, 0.7), xticks=[6, 4, 2])

# B1: N=1000, K=5, r=4.0
make_figure(x_b, {
    "Louvain":      [0.833, 0.679, 0.706, 0.816, 0.818, 0.872, 0.820, 0.675, 0.777, 0.794],
    "KMeans":       [1.000, 1.000, 0.989, 0.917, 0.777, 0.647, 0.433, 0.317, 0.288, 0.192],
    "Pure-aSBM":    [0.979, 1.000, 0.988, 0.920, 0.812, 0.658, 0.465, 0.363, 0.248, 0.142],
    "aSBM-Louvain": [0.998, 0.999, 0.996, 0.973, 0.908, 0.895, 0.876, 0.848, 0.882, 0.830],
}, xlabel="Sigma", filename="scenario_B1_attribute_noise.png")

# B2: N=1000, K=10, r=4.0
make_figure(x_b, {
    "Louvain":      [0.114, 0.077, 0.084, 0.099, 0.069, 0.116, 0.077, 0.080, 0.094, 0.072],
    "KMeans":       [0.510, 0.535, 0.510, 0.393, 0.313, 0.249, 0.150, 0.095, 0.051, 0.043],
    "Pure-aSBM":    [0.915, 0.972, 0.973, 0.811, 0.617, 0.374, 0.212, 0.082, 0.062, 0.032],
    "aSBM-Louvain": [0.894, 0.948, 0.915, 0.868, 0.544, 0.353, 0.153, 0.064, 0.055, 0.052],
}, xlabel="Sigma", filename="scenario_B2_attribute_noise.png")

# B3: N=5000, K=5, r=4.0
make_figure(x_b, {
    "Louvain":      [0.292, 0.285, 0.391, 0.425, 0.262, 0.247, 0.397, 0.289, 0.295, 0.134],
    "KMeans":       [1.000, 1.000, 0.995, 0.915, 0.774, 0.631, 0.462, 0.366, 0.191, 0.111],
    "Pure-aSBM":    [0.716, 0.657, 0.725, 0.720, 0.676, 0.675, 0.696, 0.474, 0.430, 0.351],
    "aSBM-Louvain": [0.999, 0.999, 0.996, 0.980, 0.937, 0.866, 0.659, 0.206, 0.059, 0.020],
}, xlabel="Sigma", filename="scenario_B3_attribute_noise.png")

print("\nAll 6 figures saved.")
