#!/usr/bin/env python3
"""
Concatenate downstream collation tables across multiple runs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd


def _read_if_exists(path: Path) -> pd.DataFrame:
    if path.exists() and path.is_file():
        return pd.read_csv(path)
    return pd.DataFrame()


def collate(run_dirs: List[Path], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics_all = []
    summary_all = []

    for run_dir in run_dirs:
        coll_dir = run_dir / "collation"
        metrics = _read_if_exists(coll_dir / "metrics_long.csv")
        summary = _read_if_exists(coll_dir / "metrics_summary_long.csv")

        if not metrics.empty:
            metrics_all.append(metrics)
        if not summary.empty:
            summary_all.append(summary)

    if metrics_all:
        pd.concat(metrics_all, ignore_index=True).to_csv(output_dir / "metrics_long.csv", index=False)
    else:
        pd.DataFrame().to_csv(output_dir / "metrics_long.csv", index=False)

    if summary_all:
        pd.concat(summary_all, ignore_index=True).to_csv(
            output_dir / "metrics_summary_long.csv", index=False
        )
    else:
        pd.DataFrame().to_csv(output_dir / "metrics_summary_long.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Collate downstream run tables")
    parser.add_argument("--runs", nargs="+", required=True, help="Downstream output directories")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    args = parser.parse_args()

    run_dirs = [Path(p) for p in args.runs]
    collate(run_dirs, Path(args.output))
