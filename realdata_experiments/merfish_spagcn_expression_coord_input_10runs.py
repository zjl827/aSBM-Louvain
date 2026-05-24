"""MERFISH SpaGCN real-data script using expression plus coordinate inputs.

This script follows the native SpaGCN workflow on MERFISH expression and
spatial inputs. It preserves SpaGCN's own preprocessing and
hyperparameter-search pipeline instead of converting the method to a
PCA-input workflow.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time
import warnings

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import anndata
import numpy as np
import pandas as pd
import SpaGCN as spg
import torch
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from spagcn_submission_helpers import (
    build_summary_frame,
    capture_stdout_with_timestamps,
    extract_total_epoch,
    records_to_frame,
    resolve_paths,
)


warnings.filterwarnings("ignore")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MERFISH SpaGCN script using expression and coordinate inputs.")
    parser.add_argument("--n-runs", type=int, default=10)
    parser.add_argument("--seed-base", type=int, default=100)
    parser.add_argument("--search-seed", type=int, default=100)
    parser.add_argument("--n-clusters", type=int, default=15)
    parser.add_argument("--beta", type=float, default=49.0)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--p", type=float, default=0.5)
    parser.add_argument("--search-res-start", type=float, default=0.7)
    parser.add_argument("--search-res-step", type=float, default=0.1)
    parser.add_argument("--search-res-max-epochs", type=int, default=20)
    parser.add_argument("--train-lr", type=float, default=0.05)
    parser.add_argument("--train-max-epochs", type=int, default=200)
    parser.add_argument("--train-tol", type=float, default=5e-3)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = resolve_paths(__file__)
    results_dir = paths["results_dir"]
    prefix = "merfish_spagcn_expression_coord_input_10runs"

    print("Preparing MERFISH inputs for clean SpaGCN submission runs...")

    load_start = time.perf_counter()
    coord = pd.read_csv(paths["data_root"] / "merfish19_spatial.csv", index_col=0)
    coord.columns = ["coord_x", "coord_y"]
    gene_expr = pd.read_csv(paths["data_root"] / "merfish19_counts.csv", index_col=0).T
    annotation = pd.read_csv(paths["data_root"] / "cell_type_annotation_0.19.csv", index_col=0)
    common_cells = gene_expr.index.intersection(coord.index).intersection(annotation.index)
    gene_expr = gene_expr.loc[common_cells].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    coord = coord.loc[common_cells]
    true_labels = annotation.loc[common_cells, annotation.columns[0]]
    adata = anndata.AnnData(gene_expr)
    adata.var_names_make_unique()
    adata.obs["coord_x"] = coord["coord_x"].to_numpy()
    adata.obs["coord_y"] = coord["coord_y"].to_numpy()
    adata.obs["true_label"] = true_labels.to_numpy()
    load_align_seconds = time.perf_counter() - load_start

    x_array = adata.obs["coord_x"].to_numpy()
    y_array = adata.obs["coord_y"].to_numpy()

    adjacency_start = time.perf_counter()
    adj = spg.calculate_adj_matrix(
        x=x_array,
        y=y_array,
        x_pixel=x_array,
        y_pixel=y_array,
        image=0,
        beta=args.beta,
        alpha=args.alpha,
        histology=False,
    )
    adjacency_seconds = time.perf_counter() - adjacency_start

    prefilter_start = time.perf_counter()
    spg.prefilter_genes(adata, min_cells=3)
    spg.prefilter_specialgenes(adata)
    prefilter_seconds = time.perf_counter() - prefilter_start

    search_seed = args.search_seed
    np.random.seed(search_seed)
    torch.manual_seed(search_seed)

    search_l_start = time.perf_counter()
    with capture_stdout_with_timestamps() as search_l_capture:
        l_value = spg.search_l(args.p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    search_l_seconds = time.perf_counter() - search_l_start
    search_l_log = records_to_frame(search_l_capture.records, "search_l", search_l_start)

    search_res_start = time.perf_counter()
    with capture_stdout_with_timestamps() as search_res_capture:
        res_value = spg.search_res(
            adata,
            adj,
            l_value,
            args.n_clusters,
            start=args.search_res_start,
            step=args.search_res_step,
            tol=args.train_tol,
            lr=args.train_lr,
            max_epochs=args.search_res_max_epochs,
            r_seed=search_seed,
            t_seed=search_seed,
            n_seed=search_seed,
        )
    search_res_seconds = time.perf_counter() - search_res_start
    search_res_log = records_to_frame(search_res_capture.records, "search_res", search_res_start)

    shared_timings = pd.DataFrame(
        [
            {"phase": "load_align", "seconds": load_align_seconds, "notes": f"Loaded MERFISH expression, spatial coordinates, and annotation for {len(common_cells)} cells."},
            {"phase": "adjacency", "seconds": adjacency_seconds, "notes": "Calculated the SpaGCN adjacency matrix once."},
            {"phase": "prefilter", "seconds": prefilter_seconds, "notes": "Applied prefilter_genes plus prefilter_specialgenes once as part of the native SpaGCN workflow."},
            {"phase": "search_l", "seconds": search_l_seconds, "notes": f"Dataset-level l search; selected l={l_value}."},
            {"phase": "search_res", "seconds": search_res_seconds, "notes": f"Dataset-level resolution search; selected res={res_value}."},
        ]
    )

    phase_logs = [search_l_log, search_res_log]
    run_metrics: list[dict[str, float | int]] = []
    labels_df = pd.DataFrame(
        {
            "Cell": common_cells,
            "coord_x": coord["coord_x"].to_numpy(),
            "coord_y": coord["coord_y"].to_numpy(),
            "true_label": true_labels.to_numpy(),
        }
    )
    shared_total = float(shared_timings["seconds"].sum())

    for run_id in range(1, args.n_runs + 1):
        seed = args.seed_base + run_id - 1
        print(f"\n[MERFISH][SpaGCN] Run {run_id} / {args.n_runs} with seed {seed}")

        np.random.seed(seed)
        torch.manual_seed(seed)

        clf = spg.SpaGCN()
        clf.set_l(l_value)

        train_start = time.perf_counter()
        with capture_stdout_with_timestamps() as train_capture:
            clf.train(
                adata,
                adj,
                init_spa=True,
                init="louvain",
                res=res_value,
                tol=args.train_tol,
                lr=args.train_lr,
                max_epochs=args.train_max_epochs,
            )
        train_seconds = time.perf_counter() - train_start
        train_log = records_to_frame(train_capture.records, "train", train_start, run_id=run_id)
        phase_logs.append(train_log)
        epochs_used = extract_total_epoch(train_capture.records, args.train_max_epochs)

        predict_start = time.perf_counter()
        y_pred, _ = clf.predict()
        predict_seconds = time.perf_counter() - predict_start

        metric_start = time.perf_counter()
        ari = adjusted_rand_score(true_labels, y_pred)
        nmi = normalized_mutual_info_score(true_labels, y_pred)
        metric_seconds = time.perf_counter() - metric_start

        k_found = int(pd.Series(y_pred).nunique())
        total_seconds = train_seconds + predict_seconds + metric_seconds
        shared_setup_seconds_per_run = shared_total / args.n_runs

        run_metrics.append(
            {
                "Run": run_id,
                "Seed": seed,
                "K_Found": k_found,
                "ARI": ari,
                "NMI": nmi,
                "l_value": l_value,
                "resolution": res_value,
                "epochs_used": epochs_used,
                "fit_seconds": train_seconds,
                "predict_seconds": predict_seconds,
                "metric_seconds": metric_seconds,
                "total_seconds": total_seconds,
                "seconds_per_epoch": train_seconds / max(epochs_used, 1),
                "seconds_per_100_epochs": train_seconds * 100 / max(epochs_used, 1),
                "shared_setup_seconds_per_run": shared_setup_seconds_per_run,
                "amortized_total_seconds": total_seconds + shared_setup_seconds_per_run,
            }
        )
        labels_df[f"Run{run_id}"] = y_pred.astype(int)

    run_metrics_df = pd.DataFrame(run_metrics)
    phase_log_df = pd.concat(phase_logs, ignore_index=True) if phase_logs else pd.DataFrame(columns=["run", "phase", "event_index", "elapsed_seconds", "message"])
    summary_df = build_summary_frame(
        run_metrics_df,
        shared_timings,
        extra={"selected_l": l_value, "selected_resolution": res_value},
    )

    run_metrics_df.to_csv(results_dir / f"{prefix}_metrics.csv", index=False)
    labels_df.to_csv(results_dir / f"{prefix}_labels.csv", index=False)
    shared_timings.to_csv(results_dir / f"{prefix}_shared_timings.csv", index=False)
    phase_log_df.to_csv(results_dir / f"{prefix}_phase_log.csv", index=False)
    summary_df.to_csv(results_dir / f"{prefix}_summary.csv", index=False)

    print("Finished MERFISH SpaGCN clean submission runs.")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    main()
