#!/usr/bin/env python3
"""
Batch backfill per-run availability/error reporting artifacts across completed runs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Sequence

from tool_availability_reporting import discover_run_dirs, write_run_report


REQUIRED_REPORT_FILES = [
    "process_status.tsv",
    "log_events.tsv",
    "tool_availability_per_run.tsv",
    "error_summary.csv",
    "error_summary.html",
]


def _is_fresh_enough(run_dir: Path) -> bool:
    pipeline_info = run_dir / "pipeline_info"
    return all((pipeline_info / name).exists() for name in REQUIRED_REPORT_FILES)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Backfill per-run availability/error reports across completed run directories"
    )
    parser.add_argument("--runs-root", required=True, help="Root containing results_*x run directories")
    parser.add_argument(
        "--run-glob",
        default="results_*x",
        help="Glob pattern under --runs-root used to discover run directories",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite reporting files even if the expected outputs already exist",
    )
    args = parser.parse_args(argv)

    runs_root = Path(args.runs_root)
    if not runs_root.exists():
        parser.error(f"runs root does not exist: {runs_root}")
        return 2
    if not runs_root.is_dir():
        parser.error(f"runs root is not a directory: {runs_root}")
        return 2

    run_dirs = discover_run_dirs(runs_root, args.run_glob)
    if not run_dirs:
        print("No run directories matched the requested pattern.")
        return 0

    processed = 0
    skipped = 0
    failures: List[str] = []

    for run_dir in run_dirs:
        if not args.force and _is_fresh_enough(run_dir):
            skipped += 1
            print(f"SKIP {run_dir}: reporting artifacts already present")
            continue
        try:
            write_run_report(run_dir)
            processed += 1
            print(f"OK   {run_dir}")
        except Exception as exc:
            failures.append(f"{run_dir}: {exc}")
            print(f"FAIL {run_dir}: {exc}")

    print(
        "Backfill complete: discovered=%d processed=%d skipped=%d failed=%d"
        % (len(run_dirs), processed, skipped, len(failures))
    )
    if failures:
        for line in failures:
            print(f"  - {line}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
