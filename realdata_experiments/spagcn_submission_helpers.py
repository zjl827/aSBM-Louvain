"""Shared helpers for the clean SpaGCN real-data submission scripts."""

from __future__ import annotations

import contextlib
import io
import pathlib
import re
import sys
import time
from dataclasses import dataclass

import pandas as pd


def resolve_paths(current_file: str) -> dict[str, pathlib.Path]:
    script_file = pathlib.Path(current_file).resolve()
    submission_dir = script_file.parent
    project_root = submission_dir.parent
    data_root = project_root
    results_dir = submission_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    return {
        "script_file": script_file,
        "submission_dir": submission_dir,
        "project_root": project_root,
        "data_root": data_root,
        "results_dir": results_dir,
    }


@dataclass
class TimestampedLine:
    timestamp: float
    message: str


class TimestampedStdout(io.TextIOBase):
    def __init__(self, mirror: io.TextIOBase | None = None) -> None:
        super().__init__()
        self.mirror = mirror
        self._buffer = ""
        self.records: list[TimestampedLine] = []

    def write(self, s: str) -> int:
        if self.mirror is not None:
            self.mirror.write(s)
        self._buffer += s
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            line = line.rstrip("\r")
            if line:
                self.records.append(TimestampedLine(time.perf_counter(), line))
        return len(s)

    def flush(self) -> None:
        if self.mirror is not None:
            self.mirror.flush()

    def finalize(self) -> None:
        if self._buffer.strip():
            self.records.append(TimestampedLine(time.perf_counter(), self._buffer.strip()))
        self._buffer = ""


@contextlib.contextmanager
def capture_stdout_with_timestamps() -> TimestampedStdout:
    capture = TimestampedStdout(mirror=sys.stdout)
    with contextlib.redirect_stdout(capture):
        yield capture
    capture.finalize()


def records_to_frame(records: list[TimestampedLine], phase: str, phase_start: float, run_id: int | None = None) -> pd.DataFrame:
    rows = []
    for idx, record in enumerate(records, start=1):
        rows.append(
            {
                "run": run_id,
                "phase": phase,
                "event_index": idx,
                "elapsed_seconds": record.timestamp - phase_start,
                "message": record.message,
            }
        )
    return pd.DataFrame(rows)


def extract_total_epoch(records: list[TimestampedLine], default_max_epochs: int) -> int:
    for record in reversed(records):
        match = re.search(r"Total epoch:\s*(\d+)", record.message)
        if match:
            return int(match.group(1))
    return int(default_max_epochs)


def build_summary_frame(run_metrics: pd.DataFrame, shared_timings: pd.DataFrame, extra: dict[str, object] | None = None) -> pd.DataFrame:
    shared_total = float(shared_timings["seconds"].sum()) if not shared_timings.empty else 0.0
    n_runs = max(len(run_metrics), 1)
    amortized_total = run_metrics["total_seconds"] + shared_total / n_runs

    summary = {
        "runs": len(run_metrics),
        "mean_ARI": run_metrics["ARI"].mean(),
        "sd_ARI": run_metrics["ARI"].std(ddof=1),
        "median_ARI": run_metrics["ARI"].median(),
        "mean_NMI": run_metrics["NMI"].mean(),
        "sd_NMI": run_metrics["NMI"].std(ddof=1),
        "median_NMI": run_metrics["NMI"].median(),
        "mean_fit_seconds": run_metrics["fit_seconds"].mean(),
        "median_fit_seconds": run_metrics["fit_seconds"].median(),
        "mean_total_seconds": run_metrics["total_seconds"].mean(),
        "median_total_seconds": run_metrics["total_seconds"].median(),
        "shared_setup_seconds": shared_total,
        "mean_amortized_total_seconds": amortized_total.mean(),
        "median_amortized_total_seconds": amortized_total.median(),
    }
    if extra:
        summary.update(extra)
    return pd.DataFrame([summary])
