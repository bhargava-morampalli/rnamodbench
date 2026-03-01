#!/usr/bin/env python3
"""
Generate per-run availability and error reporting artifacts for a completed run.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from tool_availability_reporting import write_run_report


def _resolve_run_dir(run_dir: Optional[str], outdir: Optional[str]) -> Path:
    selected = run_dir or outdir
    if not selected:
        raise ValueError("One of --run-dir or --outdir must be provided.")
    path = Path(selected)
    if not path.exists():
        raise ValueError(f"Run directory does not exist: {path}")
    if not path.is_dir():
        raise ValueError(f"Run directory is not a directory: {path}")
    return path


def _resolve_pipeline_info_dir(pipeline_info_dir: Optional[str]) -> Optional[Path]:
    if not pipeline_info_dir:
        return None
    path = Path(pipeline_info_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"Pipeline info directory is not a directory: {path}")
    return path


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Backfill per-run availability and error reports for a completed run"
    )
    parser.add_argument(
        "--run-dir",
        help="Completed run directory containing pipeline_info/, logs/, and modifications/",
    )
    parser.add_argument(
        "--outdir",
        help="Backward-compatible alias for --run-dir",
    )
    parser.add_argument(
        "--output",
        default="error_summary",
        help="Prefix for the HTML/CSV wrapper reports in pipeline_info/ (default: error_summary)",
    )
    parser.add_argument(
        "--pipeline-info-dir",
        help="Optional override for where report files are written (default: <run-dir>/pipeline_info)",
    )
    args = parser.parse_args(argv)

    try:
        run_dir = _resolve_run_dir(args.run_dir, args.outdir)
        pipeline_info_dir = _resolve_pipeline_info_dir(args.pipeline_info_dir)
    except ValueError as exc:
        parser.error(str(exc))
        return 2

    report = write_run_report(
        run_dir,
        output_prefix=args.output,
        pipeline_info_dir=pipeline_info_dir,
    )
    paths = report["paths"]
    availability_df = report["tool_availability"]
    events_df = report["log_events"]

    print(f"Run directory: {run_dir}")
    print(f"Process status TSV: {paths['process_status']}")
    print(f"Log events TSV: {paths['log_events']}")
    print(f"Tool availability TSV: {paths['tool_availability']}")
    print(f"Summary CSV: {paths['summary_csv']}")
    print(f"Summary HTML: {paths['summary_html']}")
    print(
        "Summary: rows=%d usable=%d error_events=%d warning_events=%d"
        % (
            len(availability_df),
            int((availability_df["output_state"] == "parsed_nonempty").sum()),
            int((events_df["event_level"] == "ERROR").sum()) if not events_df.empty else 0,
            int((events_df["event_level"] == "WARNING").sum()) if not events_df.empty else 0,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
