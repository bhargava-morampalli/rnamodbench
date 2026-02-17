#!/usr/bin/env python3
"""
Collate downstream collation tables across runs with robust diagnostics.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd
from pandas.errors import EmptyDataError


DEFAULT_TABLES = [
    "metrics_long.csv",
    "metrics_summary_long.csv",
    "window_metrics_long.csv",
    "window_metrics_summary_long.csv",
    "lag_metrics_long.csv",
    "lag_metrics_summary_long.csv",
]

DEFAULT_RUN_GLOB = "results_*/downstream_analysis"
DEFAULT_REPORT_PREFIX = "collate_report"

# Known schemas used when writing header-only files for empty tables.
DEFAULT_TABLE_SCHEMAS: Dict[str, List[str]] = {
    "metrics_long.csv": [
        "run_id",
        "coverage_label",
        "quality_label",
        "tool",
        "reference",
        "replicate",
        "auprc",
        "auroc",
        "f1_optimal",
        "precision_optimal",
        "recall_optimal",
        "optimal_threshold",
        "n_true_positives",
        "n_false_positives",
        "n_false_negatives",
        "n_true_negatives",
        "n_ground_truth",
        "n_predictions",
        "auprc_reported",
        "auroc_reported",
        "f1_reported",
        "metrics_valid",
        "invalid_score_fraction",
        "invalid_score_n",
        "invalid_score_total",
        "metric_scope_note",
        "n_universe",
        "n_reported",
        "n_no_call",
        "no_call_rate",
        "n_gt_total",
        "n_gt_reported",
        "gt_recall_raw",
    ],
    "metrics_summary_long.csv": [
        "run_id",
        "coverage_label",
        "quality_label",
        "reference",
        "tool",
        "auprc_mean",
        "auprc_std",
        "auprc_median",
        "auroc_mean",
        "auroc_std",
        "auroc_median",
        "auprc_reported_mean",
        "auprc_reported_std",
        "auprc_reported_median",
        "auroc_reported_mean",
        "auroc_reported_std",
        "auroc_reported_median",
        "f1_optimal_mean",
        "f1_optimal_std",
        "f1_optimal_median",
        "f1_reported_mean",
        "f1_reported_std",
        "f1_reported_median",
        "precision_optimal_mean",
        "recall_optimal_mean",
        "metrics_valid_fraction",
        "no_call_rate_mean",
        "gt_recall_raw_mean",
    ],
    "window_metrics_long.csv": [
        "run_id",
        "coverage_label",
        "quality_label",
        "auprc",
        "auroc",
        "f1_optimal",
        "precision_optimal",
        "recall_optimal",
        "n_true_positives",
        "n_false_positives",
        "n_false_negatives",
        "n_true_negatives",
        "optimal_threshold",
        "n_positive",
        "n_universe",
        "prevalence",
        "metric_scope_note",
        "reference",
        "tool",
        "replicate",
        "scope",
        "window_size",
        "metrics_valid",
        "invalid_score_fraction",
        "invalid_score_n",
        "invalid_score_total",
    ],
    "window_metrics_summary_long.csv": [
        "run_id",
        "coverage_label",
        "quality_label",
        "reference",
        "tool",
        "scope",
        "window_size",
        "auprc_mean",
        "auprc_std",
        "auprc_median",
        "auroc_mean",
        "auroc_std",
        "auroc_median",
        "f1_optimal_mean",
        "f1_optimal_std",
        "f1_optimal_median",
        "precision_optimal_mean",
        "recall_optimal_mean",
        "prevalence_mean",
        "n_positive_mean",
        "n_universe_mean",
        "metrics_valid_fraction",
    ],
    "lag_metrics_long.csv": [
        "run_id",
        "coverage_label",
        "quality_label",
        "auprc",
        "auroc",
        "f1_optimal",
        "precision_optimal",
        "recall_optimal",
        "n_true_positives",
        "n_false_positives",
        "n_false_negatives",
        "n_true_negatives",
        "optimal_threshold",
        "n_positive",
        "n_universe",
        "prevalence",
        "metric_scope_note",
        "reference",
        "tool",
        "replicate",
        "scope",
        "delta",
        "metrics_valid",
        "invalid_score_fraction",
        "invalid_score_n",
        "invalid_score_total",
    ],
    "lag_metrics_summary_long.csv": [
        "run_id",
        "coverage_label",
        "quality_label",
        "reference",
        "tool",
        "scope",
        "delta",
        "auprc_mean",
        "auprc_std",
        "auprc_median",
        "auroc_mean",
        "auroc_std",
        "auroc_median",
        "f1_optimal_mean",
        "f1_optimal_std",
        "f1_optimal_median",
        "precision_optimal_mean",
        "recall_optimal_mean",
        "prevalence_mean",
        "n_positive_mean",
        "n_universe_mean",
        "metrics_valid_fraction",
    ],
}


class HardCollateError(RuntimeError):
    """Raised for hard failures where processing should stop."""


@dataclass
class IssueRecord:
    severity: str
    code: str
    message: str
    run: Optional[str] = None
    table: Optional[str] = None
    path: Optional[str] = None


@dataclass
class RunRecord:
    source: str
    input_path: str
    normalized_path: Optional[str]
    run_name: Optional[str]
    status: str
    reason: str = ""


@dataclass
class TableRecord:
    table: str
    output_path: Optional[str]
    rows_out: int
    cols_out: int
    requested_runs: int
    participating_runs: int
    missing_files: int
    empty_files: int
    parse_errors: int
    schema_drift: bool
    schema_variants: int
    columns_union: List[str] = field(default_factory=list)
    columns_intersection: List[str] = field(default_factory=list)
    drift_details: Dict[str, Dict[str, List[str]]] = field(default_factory=dict)


@dataclass
class CollateReport:
    generated_at: str
    output_dir: str
    strict: bool
    run_glob: str
    tables_requested: List[str]
    runs_requested: int = 0
    runs_discovered: int = 0
    runs_total: int = 0
    runs_usable: int = 0
    runs_unusable: int = 0
    run_records: List[RunRecord] = field(default_factory=list)
    table_records: List[TableRecord] = field(default_factory=list)
    issues: List[IssueRecord] = field(default_factory=list)
    output_files: Dict[str, str] = field(default_factory=dict)

    @property
    def issue_counts(self) -> Dict[str, int]:
        counts = {"error": 0, "warning": 0, "info": 0}
        for issue in self.issues:
            counts[issue.severity] = counts.get(issue.severity, 0) + 1
        return counts

    @property
    def soft_issue_count(self) -> int:
        counts = self.issue_counts
        return counts.get("warning", 0) + counts.get("error", 0)

    def to_dict(self) -> Dict[str, object]:
        data = asdict(self)
        data["issue_counts"] = self.issue_counts
        data["soft_issue_count"] = self.soft_issue_count
        return data


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _print_issue(issue: IssueRecord, quiet: bool) -> None:
    if quiet:
        return
    details = []
    if issue.run:
        details.append(f"run={issue.run}")
    if issue.table:
        details.append(f"table={issue.table}")
    if issue.path:
        details.append(f"path={issue.path}")
    suffix = f" ({', '.join(details)})" if details else ""
    print(f"[{issue.severity.upper()}] {issue.code}: {issue.message}{suffix}", file=sys.stderr)


def _add_issue(
    report: CollateReport,
    quiet: bool,
    severity: str,
    code: str,
    message: str,
    run: Optional[str] = None,
    table: Optional[str] = None,
    path: Optional[Path] = None,
) -> None:
    issue = IssueRecord(
        severity=severity,
        code=code,
        message=message,
        run=run,
        table=table,
        path=str(path) if path else None,
    )
    report.issues.append(issue)
    _print_issue(issue, quiet=quiet)


def _run_name_from_downstream(downstream_dir: Path) -> str:
    if downstream_dir.name == "downstream_analysis":
        return downstream_dir.parent.name
    return downstream_dir.name


def _normalize_input_path(input_path: Path) -> Tuple[Optional[Path], Optional[Path], str]:
    """
    Normalize run inputs to a downstream_analysis directory and collation directory.
    """
    if not input_path.exists():
        return None, None, "Input path does not exist"
    if not input_path.is_dir():
        return None, None, "Input path is not a directory"

    if input_path.name == "collation":
        downstream_dir = input_path.parent
        collation_dir = input_path
    elif input_path.name == "downstream_analysis":
        downstream_dir = input_path
        collation_dir = input_path / "collation"
    elif (input_path / "downstream_analysis").is_dir():
        downstream_dir = input_path / "downstream_analysis"
        collation_dir = downstream_dir / "collation"
    elif (input_path / "collation").is_dir():
        downstream_dir = input_path
        collation_dir = input_path / "collation"
    else:
        return None, None, "Could not find downstream_analysis/collation under input path"

    if not collation_dir.exists() or not collation_dir.is_dir():
        return downstream_dir, None, "Collation directory not found"

    return downstream_dir, collation_dir, ""


def _discover_runs(runs_root: Path, run_glob: str) -> List[Path]:
    return sorted(p for p in runs_root.glob(run_glob))


def _ordered_union(columns_by_run: Dict[str, List[str]]) -> List[str]:
    seen = set()
    ordered = []
    for cols in columns_by_run.values():
        for col in cols:
            if col not in seen:
                seen.add(col)
                ordered.append(col)
    return ordered


def _ordered_intersection(columns_by_run: Dict[str, List[str]], union_cols: List[str]) -> List[str]:
    if not columns_by_run:
        return []
    run_sets = [set(cols) for cols in columns_by_run.values()]
    return [col for col in union_cols if all(col in run_set for run_set in run_sets)]


def _build_action_checklist(report: CollateReport) -> List[str]:
    codes = {issue.code for issue in report.issues}
    actions: List[str] = []

    if {"runs_root_missing", "missing_input_path", "invalid_input_path", "duplicate_input", "no_runs_discovered"} & codes:
        actions.append("Validate run selection inputs (`--runs`, `--runs-root`, `--run-glob`) and fix path typos.")

    if {"missing_table_file", "empty_file", "csv_parse_error"} & codes:
        actions.append(
            "Re-check downstream collation outputs for affected runs and regenerate missing/empty/malformed CSV files."
        )

    if "schema_drift" in codes:
        actions.append(
            "Inspect schema drift entries and confirm all runs were produced by compatible downstream-analysis versions."
        )

    if "skipped_table_no_schema" in codes:
        actions.append("Provide known table names or ensure at least one readable file exists to infer schema.")

    if not actions:
        actions.append("No action required.")

    return actions


def _render_markdown(report: CollateReport) -> str:
    lines: List[str] = []
    counts = report.issue_counts

    lines.append("# Collation Report")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append(f"- Generated: `{report.generated_at}`")
    lines.append(f"- Output directory: `{report.output_dir}`")
    lines.append(f"- Runs requested (explicit): **{report.runs_requested}**")
    lines.append(f"- Runs discovered: **{report.runs_discovered}**")
    lines.append(f"- Runs total (after dedupe/intake): **{report.runs_total}**")
    lines.append(f"- Runs usable: **{report.runs_usable}**")
    lines.append(f"- Runs unusable/skipped: **{report.runs_unusable}**")
    lines.append(f"- Soft issues (warning+error): **{report.soft_issue_count}**")
    lines.append("")
    lines.append("| Severity | Count |")
    lines.append("| --- | ---: |")
    lines.append(f"| error | {counts.get('error', 0)} |")
    lines.append(f"| warning | {counts.get('warning', 0)} |")
    lines.append(f"| info | {counts.get('info', 0)} |")
    lines.append("")

    lines.append("## Run Intake Table")
    lines.append("")
    lines.append("| Source | Input Path | Normalized Path | Run Name | Status | Reason |")
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for rec in report.run_records:
        lines.append(
            f"| {rec.source} | `{rec.input_path}` | `{rec.normalized_path or ''}` | "
            f"`{rec.run_name or ''}` | {rec.status} | {rec.reason or ''} |"
        )
    lines.append("")

    lines.append("## Per-Table Merge Summary")
    lines.append("")
    lines.append(
        "| Table | Rows Out | Cols Out | Runs Used | Missing Files | Empty Files | Parse Errors | Schema Drift | Output |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |")
    for rec in report.table_records:
        lines.append(
            f"| `{rec.table}` | {rec.rows_out} | {rec.cols_out} | {rec.participating_runs}/{rec.requested_runs} | "
            f"{rec.missing_files} | {rec.empty_files} | {rec.parse_errors} | {str(rec.schema_drift)} | "
            f"`{rec.output_path or ''}` |"
        )
    lines.append("")

    lines.append("## Issues by Severity")
    lines.append("")
    for severity in ("error", "warning", "info"):
        items = [i for i in report.issues if i.severity == severity]
        if not items:
            continue
        lines.append(f"### {severity.upper()} ({len(items)})")
        lines.append("")
        for issue in items:
            details = []
            if issue.run:
                details.append(f"run={issue.run}")
            if issue.table:
                details.append(f"table={issue.table}")
            if issue.path:
                details.append(f"path={issue.path}")
            suffix = f" ({', '.join(details)})" if details else ""
            lines.append(f"- `{issue.code}`: {issue.message}{suffix}")
        lines.append("")

    lines.append("## Action Checklist")
    lines.append("")
    for action in _build_action_checklist(report):
        lines.append(f"- {action}")
    lines.append("")

    return "\n".join(lines)


def _write_reports(report: CollateReport, output_dir: Path, report_prefix: str) -> None:
    json_path = output_dir / f"{report_prefix}.json"
    md_path = output_dir / f"{report_prefix}.md"
    report.output_files["report_json"] = str(json_path)
    report.output_files["report_markdown"] = str(md_path)

    json_path.write_text(json.dumps(report.to_dict(), indent=2), encoding="utf-8")
    md_path.write_text(_render_markdown(report), encoding="utf-8")


def collate(
    run_dirs: List[Path],
    output_dir: Path,
    *,
    tables: Optional[List[str]] = None,
    runs_root: Optional[Path] = None,
    run_glob: str = DEFAULT_RUN_GLOB,
    report_prefix: str = DEFAULT_REPORT_PREFIX,
    quiet: bool = False,
    strict: bool = False,
) -> CollateReport:
    table_names = tables if tables else list(DEFAULT_TABLES)
    report = CollateReport(
        generated_at=_now_iso(),
        output_dir=str(output_dir),
        strict=strict,
        run_glob=run_glob,
        tables_requested=table_names,
    )

    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as exc:  # pragma: no cover - hard to simulate reliably
        raise HardCollateError(f"Unable to create output directory '{output_dir}': {exc}") from exc

    explicit_paths = [Path(p) for p in run_dirs]
    report.runs_requested = len(explicit_paths)
    discovered_paths: List[Path] = []

    if runs_root is not None:
        root = Path(runs_root)
        if not root.exists():
            _add_issue(
                report,
                quiet=quiet,
                severity="warning",
                code="runs_root_missing",
                message="Discovery root does not exist",
                path=root,
            )
        elif not root.is_dir():
            _add_issue(
                report,
                quiet=quiet,
                severity="warning",
                code="runs_root_not_dir",
                message="Discovery root is not a directory",
                path=root,
            )
        else:
            discovered_paths = _discover_runs(root, run_glob)
            report.runs_discovered = len(discovered_paths)
            if not discovered_paths:
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="warning",
                    code="no_runs_discovered",
                    message=f"No runs discovered under root using glob '{run_glob}'",
                    path=root,
                )

    source_paths: List[Tuple[str, Path]] = [("explicit", p) for p in explicit_paths]
    source_paths.extend(("discovered", p) for p in discovered_paths)

    usable_runs: List[Tuple[str, Path, Path]] = []  # (run_name, downstream_dir, collation_dir)
    seen_normalized: Dict[str, str] = {}

    for source, input_path in source_paths:
        downstream_dir, collation_dir, reason = _normalize_input_path(input_path)
        normalized_key = (
            str((downstream_dir or input_path).resolve(strict=False))
            if downstream_dir is not None or input_path is not None
            else str(input_path)
        )
        run_name = _run_name_from_downstream(downstream_dir) if downstream_dir else None

        if normalized_key in seen_normalized:
            rec = RunRecord(
                source=source,
                input_path=str(input_path),
                normalized_path=str(downstream_dir) if downstream_dir else None,
                run_name=run_name,
                status="duplicate",
                reason=f"Duplicate of {seen_normalized[normalized_key]}",
            )
            report.run_records.append(rec)
            _add_issue(
                report,
                quiet=quiet,
                severity="warning",
                code="duplicate_input",
                message="Duplicate run input after normalization; skipping duplicate",
                run=run_name,
                path=input_path,
            )
            continue

        seen_normalized[normalized_key] = str(input_path)

        if downstream_dir is None or collation_dir is None:
            status = "missing" if "does not exist" in reason else "invalid"
            rec = RunRecord(
                source=source,
                input_path=str(input_path),
                normalized_path=str(downstream_dir) if downstream_dir else None,
                run_name=run_name,
                status=status,
                reason=reason,
            )
            report.run_records.append(rec)
            _add_issue(
                report,
                quiet=quiet,
                severity="warning",
                code="missing_input_path" if status == "missing" else "invalid_input_path",
                message=reason,
                run=run_name,
                path=input_path,
            )
            continue

        rec = RunRecord(
            source=source,
            input_path=str(input_path),
            normalized_path=str(downstream_dir),
            run_name=run_name,
            status="usable",
            reason="",
        )
        report.run_records.append(rec)
        usable_runs.append((run_name, downstream_dir, collation_dir))

    report.runs_total = len(report.run_records)
    report.runs_usable = len(usable_runs)
    report.runs_unusable = len([r for r in report.run_records if r.status != "usable"])

    if not usable_runs:
        _add_issue(
            report,
            quiet=quiet,
            severity="error",
            code="no_usable_runs",
            message="No usable runs resolved from inputs",
        )
        _write_reports(report, output_dir=output_dir, report_prefix=report_prefix)
        raise HardCollateError("No usable runs resolved from inputs")

    for table_name in table_names:
        frames: List[pd.DataFrame] = []
        columns_by_run: Dict[str, List[str]] = {}
        missing_files = 0
        empty_files = 0
        parse_errors = 0

        for run_name, _downstream_dir, collation_dir in usable_runs:
            table_path = collation_dir / table_name

            if not table_path.exists() or not table_path.is_file():
                missing_files += 1
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="warning",
                    code="missing_table_file",
                    message="Table file not found for run",
                    run=run_name,
                    table=table_name,
                    path=table_path,
                )
                continue

            if table_path.stat().st_size == 0:
                empty_files += 1
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="warning",
                    code="empty_file",
                    message="Table file is empty (0 bytes)",
                    run=run_name,
                    table=table_name,
                    path=table_path,
                )
                continue

            try:
                df = pd.read_csv(table_path)
            except EmptyDataError:
                empty_files += 1
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="warning",
                    code="empty_file",
                    message="CSV file has no readable content",
                    run=run_name,
                    table=table_name,
                    path=table_path,
                )
                continue
            except Exception as exc:
                parse_errors += 1
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="warning",
                    code="csv_parse_error",
                    message=f"Failed to parse CSV: {exc}",
                    run=run_name,
                    table=table_name,
                    path=table_path,
                )
                continue

            columns_by_run[run_name] = list(df.columns)
            frames.append(df)

        union_cols = _ordered_union(columns_by_run)
        intersection_cols = _ordered_intersection(columns_by_run, union_cols)
        schema_variants = len({tuple(cols) for cols in columns_by_run.values()}) if columns_by_run else 0
        schema_drift = schema_variants > 1
        drift_details: Dict[str, Dict[str, List[str]]] = {}

        output_path: Optional[Path] = None
        rows_out = 0
        cols_out = 0

        if frames:
            baseline_cols = next(iter(columns_by_run.values()))
            baseline_set = set(baseline_cols)

            if schema_drift:
                for run_name, cols in columns_by_run.items():
                    cols_set = set(cols)
                    missing_vs_union = [c for c in union_cols if c not in cols_set]
                    extra_vs_baseline = [c for c in cols if c not in baseline_set]
                    if missing_vs_union or extra_vs_baseline:
                        drift_details[run_name] = {
                            "missing_vs_union": missing_vs_union,
                            "extra_vs_baseline": extra_vs_baseline,
                        }
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="warning",
                    code="schema_drift",
                    message=f"Schema drift detected across runs ({schema_variants} schema variants)",
                    table=table_name,
                )

            normalized_frames = [df.reindex(columns=union_cols) for df in frames]
            merged = pd.concat(normalized_frames, ignore_index=True, sort=False)
            output_path = output_dir / table_name
            merged.to_csv(output_path, index=False)
            rows_out = len(merged)
            cols_out = len(merged.columns)
            report.output_files[table_name] = str(output_path)
        else:
            known_schema = DEFAULT_TABLE_SCHEMAS.get(table_name, [])
            if known_schema:
                output_path = output_dir / table_name
                pd.DataFrame(columns=known_schema).to_csv(output_path, index=False)
                rows_out = 0
                cols_out = len(known_schema)
                union_cols = known_schema
                intersection_cols = known_schema
                report.output_files[table_name] = str(output_path)
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="info",
                    code="header_only_output_written",
                    message="No readable input frames; wrote header-only output using known schema",
                    table=table_name,
                    path=output_path,
                )
            else:
                _add_issue(
                    report,
                    quiet=quiet,
                    severity="warning",
                    code="skipped_table_no_schema",
                    message="No readable input frames and no known schema; skipped writing table",
                    table=table_name,
                )

        report.table_records.append(
            TableRecord(
                table=table_name,
                output_path=str(output_path) if output_path else None,
                rows_out=rows_out,
                cols_out=cols_out,
                requested_runs=len(usable_runs),
                participating_runs=len(frames),
                missing_files=missing_files,
                empty_files=empty_files,
                parse_errors=parse_errors,
                schema_drift=schema_drift,
                schema_variants=schema_variants,
                columns_union=union_cols,
                columns_intersection=intersection_cols,
                drift_details=drift_details,
            )
        )

    _write_reports(report, output_dir=output_dir, report_prefix=report_prefix)

    counts = report.issue_counts
    print(
        f"Collation complete: runs_usable={report.runs_usable}, "
        f"tables={len(report.table_records)}, warnings={counts.get('warning', 0)}, "
        f"errors={counts.get('error', 0)}, info={counts.get('info', 0)}"
    )
    print(f"Report JSON: {report.output_files.get('report_json', '')}")
    print(f"Report Markdown: {report.output_files.get('report_markdown', '')}")

    return report


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Collate downstream run tables with diagnostics")
    parser.add_argument("--runs", nargs="+", help="Run paths (results_*, downstream_analysis, or collation)")
    parser.add_argument("--runs-root", help="Root directory for run discovery")
    parser.add_argument(
        "--run-glob",
        default=DEFAULT_RUN_GLOB,
        help=f"Glob pattern used under --runs-root (default: {DEFAULT_RUN_GLOB})",
    )
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--tables", nargs="+", help="Optional table names to collate")
    parser.add_argument(
        "--report-prefix",
        default=DEFAULT_REPORT_PREFIX,
        help=f"Prefix for report files (default: {DEFAULT_REPORT_PREFIX})",
    )
    parser.add_argument("--strict", action="store_true", help="Return non-zero if any soft issues are found")
    parser.add_argument("--quiet", action="store_true", help="Suppress per-issue console output")
    args = parser.parse_args(argv)

    if not args.runs and not args.runs_root:
        parser.error("At least one of --runs or --runs-root must be provided.")

    try:
        report = collate(
            run_dirs=[Path(p) for p in (args.runs or [])],
            output_dir=Path(args.output),
            tables=args.tables,
            runs_root=Path(args.runs_root) if args.runs_root else None,
            run_glob=args.run_glob,
            report_prefix=args.report_prefix,
            quiet=args.quiet,
            strict=args.strict,
        )
    except HardCollateError as exc:
        print(f"HARD FAILURE: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover - safety net
        print(f"UNEXPECTED FAILURE: {exc}", file=sys.stderr)
        return 2

    if args.strict and report.soft_issue_count > 0:
        print(
            f"Strict mode enabled and {report.soft_issue_count} soft issues detected; returning non-zero.",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
