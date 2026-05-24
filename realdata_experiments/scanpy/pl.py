from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def scatter(adata, x, y, color, alpha=1.0, title=None, color_map=None, show=False, size=20):
    fig, ax = plt.subplots()
    x_vals = np.asarray(adata.obs[x])
    y_vals = np.asarray(adata.obs[y])
    c_vals = np.asarray(adata.obs[color])
    scatter_kwargs = {"alpha": alpha, "s": size}

    if np.issubdtype(c_vals.dtype, np.number):
        scatter_kwargs["c"] = c_vals
        if color_map is not None:
            scatter_kwargs["cmap"] = color_map
    else:
        categories, codes = np.unique(c_vals.astype(str), return_inverse=True)
        scatter_kwargs["c"] = codes
        if color_map is not None:
            scatter_kwargs["cmap"] = color_map

    ax.scatter(x_vals, y_vals, **scatter_kwargs)
    if title is not None:
        ax.set_title(title)
    if show:
        plt.show()
    return ax
