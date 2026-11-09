#!/usr/bin/env python3
"""Audit and summarize the 21-condition paired data-generation experiment.

All uncertainty is across independently generated replicate indices. Conditions
within a replicate are paired, and are never treated as independent replicates.
The script writes incomplete/error counts before returning a nonzero exit code.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

METHODS = ["Proposed", "Pure-aSBM", "Louvain", "K-means"]
COLORS = {"Proposed": "#235f9c", "Pure-aSBM": "#c56334", "Louvain": "#478466", "K-means": "#875ba7"}
NA_TEXT = {"", "na", "nan", "none", "null"}
METRICS = ["ARI", "NMI", "initial_K", "final_K", "iterations", "elapsed_sec",
           "high_F1", "low_F1", "worst_F1", "high_precision", "low_precision",
           "worst_precision", "high_recall", "low_recall", "worst_recall"]
ALIASES = {
    "replicate": "rep", "replicate_id": "rep", "ari": "ARI", "nmi": "NMI",
    "init_K": "initial_K", "K_initial": "initial_K", "K_final": "final_K",
    "runtime_sec": "elapsed_sec", "elapsed_seconds": "elapsed_sec",
    "iterations_run": "iterations", "f1": "F1", "is_high": "is_high_noise",
    "true_community": "community_id", "reached_stopping_rule": "met_stopping_rule",
    "algorithm": "method", "warning": "warnings",
}


def as_bool(value):
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "t"}:
        return True
    if text in {"false", "0", "no", "f"}:
        return False
    return None


def present_text(value):
    return not pd.isna(value) and str(value).strip().lower() not in NA_TEXT


def read_csv(path):
    frame = pd.read_csv(path)
    frame = frame.rename(columns={a: b for a, b in ALIASES.items() if a in frame and b not in frame})
    if "method" in frame:
        frame["method"] = frame["method"].replace({
            "aSBM-Louvain": "Proposed", "aSBM--Louvain": "Proposed", "aSBM_Louvain": "Proposed",
            "Pure_ASBM": "Pure-aSBM", "Pure_aSBM": "Pure-aSBM", "Pure-aSBM": "Pure-aSBM", "KMeans": "K-means", "KMeans_estK": "K-means",
            "k-means": "K-means", "Kmeans": "K-means", "kmeans": "K-means",
        })
    if "rep" in frame:
        frame["rep"] = pd.to_numeric(frame["rep"], errors="raise").astype(int)
    return frame


def stats(values):
    clean = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    clean = clean[np.isfinite(clean)]
    count = len(clean)
    mean = float(np.mean(clean)) if count else np.nan
    sd = float(np.std(clean, ddof=1)) if count > 1 else np.nan
    mcse = sd / math.sqrt(count) if count > 1 else np.nan
    return {"n": count, "mean": mean, "sd": sd, "mcse": mcse}


def error_mask(frame):
    err = frame["error"].map(present_text) if "error" in frame else pd.Series(False, index=frame.index)
    if "status" in frame:
        err |= frame["status"].astype(str).str.lower().isin({"error", "failed", "failure", "task_error"})
    return err


def add_local_metrics(raw, local, issues):
    if local is None:
        if not all(m in raw for m in ["high_F1", "low_F1", "worst_F1"]):
            issues.append("Community results are absent; high/low/worst-community F1 unavailable.")
        return raw, pd.DataFrame()
    required = {"condition_id", "rep", "method", "community_id", "is_high_noise", "F1", "precision", "recall"}
    missing = required.difference(local.columns)
    if missing:
        issues.append("Community results missing columns: " + ", ".join(sorted(missing)))
        return raw, pd.DataFrame()
    keys = ["condition_id", "rep", "method"]
    if local.duplicated(keys + ["community_id"]).any():
        issues.append("Duplicate community-result keys; local metrics were not aggregated.")
        return raw, pd.DataFrame()
    local = local.copy()
    local["is_high_noise"] = local["is_high_noise"].map(as_bool)
    records = []
    for key, group in local.groupby(keys, sort=True):
        row = dict(zip(keys, key))
        row["local_communities"] = len(group)
        row["local_high_communities"] = int(group["is_high_noise"].eq(True).sum())
        row["local_low_communities"] = int(group["is_high_noise"].eq(False).sum())
        for metric in ["F1", "precision", "recall"]:
            vals = pd.to_numeric(group[metric], errors="coerce")
            # Never silently average a partial set of true communities.
            complete = len(group) == 5 and vals.notna().all() and np.isfinite(vals).all()
            if complete and ((vals < 0) | (vals > 1)).any():
                issues.append(f"Out-of-range {metric} in community results for {key}.")
                complete = False
            for category, expected_n, flag in [("high", 2, True), ("low", 3, False)]:
                subset = vals[group["is_high_noise"].eq(flag)]
                row[f"{category}_{metric}"] = float(subset.mean()) if complete and len(subset) == expected_n else np.nan
            row[f"worst_{metric}"] = float(vals.min()) if complete else np.nan
        records.append(row)
    aggregates = pd.DataFrame(records)
    overlaps = [m for m in METRICS if m in raw and m in aggregates]
    if overlaps:
        raw = raw.drop(columns=overlaps)
    return raw.merge(aggregates, on=keys, how="left", validate="one_to_one"), aggregates


def summarize(raw, grid, repetitions):
    rows = []
    for _, condition in grid.iterrows():
        for method in METHODS:
            block = raw[(raw.condition_id == condition.condition_id) & (raw.method == method)]
            failures = error_mask(block)
            usable = block.loc[~failures]
            row = {"condition_id": condition.condition_id, "method": method,
                   "r": condition.r, "sigma_rms": condition.sigma_rms,
                   "high_low_sd_ratio": condition.high_low_sd_ratio,
                   "in_joint_grid": as_bool(condition.in_joint_grid),
                   "in_heterogeneous_noise": as_bool(condition.in_heterogeneous_noise),
                   "n_expected": repetitions, "n_present": len(block),
                   "n_failed": int(failures.sum()), "n_missing": max(0, repetitions-len(block)),
                   "failure_fraction_of_expected": float(failures.sum())/repetitions,
                   "n_warnings": int(block["warnings"].map(present_text).sum()) if "warnings" in block else 0}
            if "met_stopping_rule" in block:
                flags = block.met_stopping_rule.map(as_bool)
                row.update(n_stopping_assessed=int(flags.notna().sum()),
                           n_met_stopping_rule=int(flags.eq(True).sum()),
                           stopping_fraction_assessed=float(flags.dropna().mean()) if flags.notna().any() else np.nan)
            if "stop_reason" in block:
                row["stop_reason_counts"] = json.dumps(block.stop_reason.fillna("not_recorded").value_counts().to_dict(), sort_keys=True)
            for metric in METRICS:
                if metric in usable:
                    row.update({f"{metric}_{name}": value for name, value in stats(usable[metric]).items()})
            rows.append(row)
    return pd.DataFrame(rows)


def pair_differences(raw, grid, repetitions):
    good = raw.loc[~error_mask(raw)]
    pairs = good.loc[good.method.isin(["Proposed", "Pure-aSBM"])].pivot(index=["condition_id", "rep"], columns="method", values="ARI")
    for method in ["Proposed", "Pure-aSBM"]:
        if method not in pairs:
            pairs[method] = np.nan
    pairs["difference"] = pairs["Proposed"] - pairs["Pure-aSBM"]
    records = []
    for _, condition in grid.iterrows():
        block = pairs.loc[pairs.index.get_level_values("condition_id") == condition.condition_id]
        values = block.difference.dropna()
        values = values[np.isfinite(values)]
        detail = stats(values)
        half = student_t.ppf(.975, detail["n"] - 1) * detail["mcse"] if detail["n"] > 1 else np.nan
        records.append({"condition_id": condition.condition_id, "r": condition.r,
                        "sigma_rms": condition.sigma_rms, "high_low_sd_ratio": condition.high_low_sd_ratio,
                        "n_expected_pairs": repetitions, "n_pairs": detail["n"],
                        "n_missing_or_failed_pairs": repetitions-detail["n"],
                        "difference_mean": detail["mean"], "difference_sd": detail["sd"],
                        "difference_mcse": detail["mcse"], "difference_ci95_low": detail["mean"]-half,
                        "difference_ci95_high": detail["mean"]+half,
                        "proposed_wins": int((values > 1e-12).sum()),
                        "ties": int((values.abs() <= 1e-12).sum()),
                        "pure_wins": int((values < -1e-12).sum()),
                        "interval_method": "two-sided Student-t across paired independent replicate indices"})
    return pd.DataFrame(records), pairs.reset_index()


def audit(raw, grid, repetitions, issues):
    keys = ["condition_id", "rep", "method"]
    expected = pd.MultiIndex.from_product([grid.condition_id, range(1, repetitions+1), METHODS], names=keys)
    actual = pd.MultiIndex.from_frame(raw[keys])
    missing = expected.difference(actual)
    extra = actual.difference(expected)
    duplicates = raw.loc[raw.duplicated(keys, keep=False), keys]
    if len(missing): issues.append(f"Missing {len(missing)} of {len(expected)} expected method-result rows.")
    if len(extra): issues.append(f"Unexpected condition/replicate/method combinations: {len(extra)}.")
    if len(duplicates): issues.append(f"Duplicate method-result rows: {len(duplicates)}.")
    for field, bounds in {"ARI": (-1, 1), "NMI": (0, 1)}.items():
        if field in raw:
            vals = pd.to_numeric(raw[field], errors="coerce")
            if ((vals.dropna() < bounds[0]-1e-10) | (vals.dropna() > bounds[1]+1e-10)).any():
                issues.append(f"Out-of-range {field} values found.")
            absent = ~np.isfinite(vals) & ~error_mask(raw)
            if absent.any():
                issues.append(f"{int(absent.sum())} non-error method rows have no finite {field}.")
    checks = []
    for column in ["data_seed", "graph_seed", "means_seed", "residual_seed", "assignment_seed", "algorithm_seed", "initial_K", "initial_partition_hash", "initialization_hash", "postinit_rng_hash", "post_init_rng_hash", "post_initialization_rng_hash", "shared_initial_K"]:
        if column in raw:
            pp = raw.loc[raw.method.isin(["Proposed", "Pure-aSBM"])]
            mismatch = int((pp.groupby(["condition_id", "rep"])[column].nunique(dropna=False) > 1).sum())
            checks.append({"check": f"Proposed/Pure shared {column}", "mismatched_groups": mismatch})
            if mismatch: issues.append(f"Proposed/Pure pairing mismatch in {column}: {mismatch} groups.")
    for column in ["graph_hash", "attribute_hash", "noise_sd", "high_groups"]:
        if column in raw:
            mismatch = int((raw.groupby(["condition_id", "rep"])[column].nunique(dropna=False) > 1).sum())
            checks.append({"check": f"All methods share condition/rep {column}", "mismatched_groups": mismatch})
            if mismatch: issues.append(f"Methods received inconsistent {column}: {mismatch} groups.")
    merged = raw.merge(grid[["condition_id", "r"]], on="condition_id", how="left", suffixes=("_recorded", ""))
    for column in ["graph_hash", "A_hash"]:
        if column in merged:
            mismatch = int((merged.groupby(["r", "rep"])[column].nunique(dropna=False) > 1).sum())
            checks.append({"check": f"Graph reuse by r/rep: {column}", "mismatched_groups": mismatch})
            if mismatch: issues.append(f"Graph pairing mismatch across noise levels: {mismatch} groups.")
    for column in ["high_noise_ids", "high_noise_community_ids", "high_groups"]:
        if column in raw:
            mismatch = int((raw.groupby("rep")[column].nunique(dropna=False) > 1).sum())
            checks.append({"check": f"High-noise groups fixed by rep: {column}", "mismatched_groups": mismatch})
            if mismatch: issues.append(f"High-noise group identities differ within {mismatch} replicates.")
    versions = sorted(raw.algorithm_version.dropna().astype(str).unique()) if "algorithm_version" in raw else []
    if len(versions)>1:issues.append("Mixed algorithm versions in one pilot analysis: "+", ".join(versions))
    warnings = raw.loc[raw.warnings.map(present_text)] if "warnings" in raw else raw.iloc[:0]
    init_warnings = warnings.loc[warnings.method=="K-means"]
    proposed = raw.loc[raw.method=="Proposed"]
    report = {"expected_conditions": len(grid), "expected_instances": len(grid)*repetitions,
              "repetitions_per_condition": repetitions,
              "expected_method_rows": len(expected), "observed_method_rows": len(raw),
              "observed_unique_instances": len(raw[["condition_id", "rep"]].drop_duplicates()),
              "missing_method_rows": len(missing), "unexpected_method_rows": len(extra),
              "duplicate_rows": len(duplicates), "failed_method_rows": int(error_mask(raw).sum()),
              "algorithm_versions": versions,
              "warning_method_rows": len(warnings),
              "warning_instances": len(warnings[["condition_id","rep"]].drop_duplicates()),
              "initialization_warning_instances": len(init_warnings),
              "warning_selected_kmeans_iterations_above_10": int((pd.to_numeric(init_warnings.iterations,errors="coerce")>10).sum()) if "iterations" in init_warnings else 0,
              "proposed_iteration_median": float(proposed.iterations.median()) if "iterations" in proposed else None,
              "proposed_iteration_max": float(proposed.iterations.max()) if "iterations" in proposed else None,
              "pairing_checks": checks, "issues": issues,
              "missing_keys": [list(x) for x in missing], "unexpected_keys": [list(x) for x in extra]}
    return report


def save_fig(fig, out, name):
    fig.savefig(out / f"{name}.png", dpi=190, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def uncertainty_limits(frame, metrics):
    lows, highs = [-.04], [1.04]
    for metric in metrics:
        if f"{metric}_mean" not in frame or f"{metric}_sd" not in frame:
            continue
        mean, sd = frame[f"{metric}_mean"], frame[f"{metric}_sd"].fillna(0)
        if (mean-sd).notna().any(): lows.append(float((mean-sd).min())-.025)
        if (mean+sd).notna().any(): highs.append(float((mean+sd).max())+.025)
    return min(lows),max(highs)


def plots(summary, differences, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    joint = summary.loc[summary.in_joint_grid.eq(True)]
    r_values = sorted(joint.r.unique())
    sigmas = sorted(joint.sigma_rms.unique())
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.2), constrained_layout=True)
    joint_ids = set(joint.condition_id)
    for ax, method in zip(axes[:2], ["Proposed", "Pure-aSBM"]):
        matrix = joint.loc[joint.method == method].pivot(index="r", columns="sigma_rms", values="ARI_mean").reindex(index=r_values, columns=sigmas)
        im = ax.imshow(matrix, vmin=0, vmax=1, cmap="YlGnBu", aspect="auto")
        for (i, j), value in np.ndenumerate(matrix.to_numpy()):
            ax.text(j, i, f"{value:.3f}" if np.isfinite(value) else "NA", ha="center", va="center", color="white" if value > .65 else "#202020")
        ax.set_title(f"{method}: mean ARI")
        fig.colorbar(im, ax=ax, shrink=.85)
    matrix = differences.loc[differences.condition_id.isin(joint_ids)].pivot(index="r", columns="sigma_rms", values="difference_mean").reindex(index=r_values, columns=sigmas)
    limit = max(.1, float(np.nanmax(np.abs(matrix.to_numpy())))) if np.isfinite(matrix.to_numpy()).any() else 1
    im = axes[2].imshow(matrix, vmin=-limit, vmax=limit, cmap="RdBu", aspect="auto")
    for (i, j), value in np.ndenumerate(matrix.to_numpy()):
        axes[2].text(j, i, f"{value:+.3f}" if np.isfinite(value) else "NA", ha="center", va="center", color="white" if abs(value) > .68*limit else "#202020")
    axes[2].set_title("Paired ARI: Proposed minus Pure")
    fig.colorbar(im, ax=axes[2], shrink=.85)
    for ax in axes:
        ax.set_xticks(range(len(sigmas)), [str(x) for x in sigmas]); ax.set_xlabel("Attribute noise sigma")
        ax.set_yticks(range(len(r_values)), [str(x) for x in r_values]); ax.set_ylabel("Structural ratio r")
    fig.suptitle(f"{int(summary.n_expected.iloc[0])}-replicate study of joint graph and attribute signal")
    save_fig(fig, out, "joint_ari_heatmaps")
    hetero = summary.loc[summary.in_heterogeneous_noise.eq(True)]
    rs = sorted(hetero.r.unique())
    fig, axes = plt.subplots(1, len(rs), figsize=(12.5, 4), sharey=True, constrained_layout=True)
    for ax, r in zip(np.atleast_1d(axes), rs):
        for method in METHODS:
            block = hetero.loc[(hetero.r == r) & (hetero.method == method)].sort_values("high_low_sd_ratio")
            ax.errorbar(block.high_low_sd_ratio, block.ARI_mean, yerr=block.ARI_sd,
                        color=COLORS[method], marker="o", markersize=4, capsize=3, label=method)
        ax.set_title(f"r = {r:g}"); ax.set_xlabel("High/low noise SD ratio h")
        ax.set_xticks([1, 1.5, 2, 3]); ax.set_ylim(*uncertainty_limits(hetero,["ARI"])); ax.grid(axis="y", alpha=.2)
    np.atleast_1d(axes)[0].set_ylabel("ARI (mean +/- one SD)")
    np.atleast_1d(axes)[-1].legend(fontsize=9, loc="best")
    fig.suptitle("Heterogeneous Gaussian noise; population-average variance fixed at 9")
    save_fig(fig, out, "heterogeneous_ari")
    if all(f"{name}_F1_mean" in hetero for name in ["high", "low", "worst"]):
        fig, axes = plt.subplots(2, len(rs), figsize=(12.5, 7.4), sharex=True, sharey=True, constrained_layout=True)
        for i, method in enumerate(["Proposed", "Pure-aSBM"]):
            for j, r in enumerate(rs):
                ax = axes[i,j]
                block = hetero.loc[(hetero.r == r) & (hetero.method == method)].sort_values("high_low_sd_ratio")
                for category, color, marker in [("high", "#be4b3c", "o"), ("low", "#27689c", "s"), ("worst", "#555555", "^")]:
                    ax.errorbar(block.high_low_sd_ratio, block[f"{category}_F1_mean"], yerr=block[f"{category}_F1_sd"],
                                color=color, marker=marker, markersize=4, capsize=3, label=f"{category} noise" if category != "worst" else "worst community")
                ax.set_title(f"{method}, r = {r:g}"); ax.set_xticks([1,1.5,2,3]); ax.set_ylim(*uncertainty_limits(hetero,["high_F1","low_F1","worst_F1"])); ax.grid(axis="y",alpha=.2)
                if i == 1: ax.set_xlabel("High/low noise SD ratio h")
                if j == 0: ax.set_ylabel("Matched community F1 +/- one SD")
        axes[0,-1].legend(fontsize=9,loc="best")
        fig.suptitle("Local recovery: fixed high-noise identities; at h=1 both groups have equal noise")
        save_fig(fig,out,"heterogeneous_community_f1")
    # All four methods remain visible in an additional joint-grid curve figure.
    fig, axes = plt.subplots(1,len(r_values),figsize=(14,3.8),sharey=True,constrained_layout=True)
    for ax,r in zip(np.atleast_1d(axes),r_values):
        for method in METHODS:
            block=joint.loc[(joint.r==r)&(joint.method==method)].sort_values("sigma_rms")
            ax.errorbar(block.sigma_rms,block.ARI_mean,yerr=block.ARI_sd,color=COLORS[method],marker="o",markersize=4,capsize=3,label=method)
        ax.set_title(f"r = {r:g}");ax.set_xlabel("Attribute noise sigma");ax.set_xticks(sigmas);ax.grid(axis="y",alpha=.2);ax.set_ylim(*uncertainty_limits(joint,["ARI"]))
    np.atleast_1d(axes)[0].set_ylabel("ARI (mean +/- one SD)")
    np.atleast_1d(axes)[-1].legend(fontsize=8,loc="best")
    fig.suptitle("Joint graph and attribute signal")
    save_fig(fig,out,"joint_ari_curves")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path)
    parser.add_argument("--community-results",type=Path)
    parser.add_argument("--grid",type=Path,default=Path(__file__).resolve().with_name("grid.csv"))
    parser.add_argument("--out",type=Path,default=Path(os.environ.get("PILOT_OUTPUT", Path(__file__).resolve().parents[2] / "outputs/joint_signals")))
    parser.add_argument("--repetitions",type=int,default=10)
    parser.add_argument("--allow-incomplete",action="store_true")
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    raw=read_csv(args.results or args.out / "raw_metrics.csv");grid=read_csv(args.grid);issues=[]
    required={"condition_id","rep","method","ARI"}
    if required.difference(raw):
        raise ValueError("Missing required raw columns: "+", ".join(sorted(required.difference(raw))))
    for metric in METRICS:
        if metric in raw:raw[metric]=pd.to_numeric(raw[metric],errors="coerce")
    report=audit(raw,grid,args.repetitions,issues)
    # Duplicate keys are ambiguous and must never be silently de-duplicated.
    if report["duplicate_rows"]:
        (args.out/"analysis_audit.json").write_text(json.dumps(report,indent=2)+"\n")
        raise SystemExit("Duplicate result keys; see analysis_audit.json")
    local_path=args.community_results or args.out / "community_metrics.csv"
    local=read_csv(local_path) if local_path.exists() else None
    raw,local_aggregates=add_local_metrics(raw,local,issues)
    if local is not None and "high_F1" in raw:
        local_missing = raw["high_F1"].isna() & ~error_mask(raw)
        if local_missing.any():issues.append(f"{int(local_missing.sum())} non-error runs lack complete matched-community F1 metrics.")
    summary=summarize(raw,grid,args.repetitions)
    differences,pairs=pair_differences(raw,grid,args.repetitions)
    summary.to_csv(args.out/"condition_method_summary.csv",index=False)
    differences.to_csv(args.out/"paired_proposed_minus_pure.csv",index=False)
    pairs.to_csv(args.out/"paired_run_differences.csv",index=False)
    if len(local_aggregates):local_aggregates.to_csv(args.out/"community_group_metrics_per_run.csv",index=False)
    raw.loc[error_mask(raw)].to_csv(args.out/"failure_records.csv",index=False)
    raw.loc[(raw.method=="K-means") & raw.warnings.map(present_text)].to_csv(args.out/"initialization_warning_records.csv",index=False)
    report["community_metric_rows"]=len(local) if local is not None else None
    report["issues"]=issues
    (args.out/"analysis_audit.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    plots(summary,differences,args.out)
    print(json.dumps({k:v for k,v in report.items() if k not in {"missing_keys","unexpected_keys"}},ensure_ascii=False,indent=2))
    if issues and not args.allow_incomplete:raise SystemExit(2)


if __name__=="__main__":
    main()
