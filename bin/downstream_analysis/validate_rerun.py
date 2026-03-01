#!/usr/bin/env python3
"""Validate a full downstream rerun against an existing collated baseline."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLLATION_FILES = [
    "metrics_long.csv",
    "metrics_summary_long.csv",
    "window_metrics_long.csv",
    "window_metrics_summary_long.csv",
    "lag_metrics_long.csv",
    "lag_metrics_summary_long.csv",
]

PARSER_REASON_BY_TOOL = {
    "differr": "differr_score_semantics",
    "xpore": "xpore_bh_fdr_and_shift",
    "epinano": "epinano_delta_error_score",
    "yanocomp": "yanocomp_center_shift",
    "nanocompore": "nanocompore_center_shift",
    "drummer": "drummer_score_priority",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", required=True, help="Coverage runs root directory")
    parser.add_argument("--old-collated", required=True, help="Existing collated baseline directory")
    parser.add_argument("--new-collated", required=True, help="New collated rerun directory")
    parser.add_argument(
        "--run-glob",
        default="results_*/downstream_20260301_parserfix_rernaeval",
        help="Glob for per-run rerun directories under runs root",
    )
    parser.add_argument(
        "--spotcheck-run",
        default="results_1000x",
        help="Coverage run used for parser spotchecks",
    )
    return parser.parse_args()


def _coerce_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    values = series.astype(str).str.strip().str.lower()
    return values.map({"true": True, "false": False}).fillna(False)


def _line_count_csv(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        total = sum(1 for _ in handle)
    return max(total - 1, 0)


def _load_header(path: Path) -> list[str]:
    return list(pd.read_csv(path, nrows=0).columns)


def compare_schemas(old_dir: Path, new_dir: Path) -> pd.DataFrame:
    records = []
    for filename in REQUIRED_COLLATION_FILES:
        old_cols = _load_header(old_dir / filename)
        new_cols = _load_header(new_dir / filename)
        records.append(
            {
                "file": filename,
                "schema_match": old_cols == new_cols,
                "old_column_count": len(old_cols),
                "new_column_count": len(new_cols),
                "old_columns": "|".join(old_cols),
                "new_columns": "|".join(new_cols),
            }
        )
    return pd.DataFrame(records)


def build_metrics_valid_fraction(new_collated: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = pd.read_csv(new_collated / "metrics_long.csv")
    metrics["metrics_valid"] = _coerce_bool(metrics["metrics_valid"])
    metrics["no_call_formula_ok"] = (
        pd.to_numeric(metrics["n_no_call"], errors="coerce")
        == pd.to_numeric(metrics["n_universe"], errors="coerce")
        - pd.to_numeric(metrics["n_reported"], errors="coerce")
    )
    grouped = (
        metrics.groupby(["coverage_label", "reference", "tool"], dropna=False)
        .agg(
            metrics_valid_fraction=("metrics_valid", "mean"),
            n_rows=("metrics_valid", "size"),
            n_valid=("metrics_valid", "sum"),
            invalid_score_fraction_mean=("invalid_score_fraction", "mean"),
            n_universe_min=("n_universe", "min"),
            n_universe_max=("n_universe", "max"),
            n_reported_mean=("n_reported", "mean"),
            n_no_call_mean=("n_no_call", "mean"),
            no_call_rate_mean=("no_call_rate", "mean"),
            no_call_formula_ok_fraction=("no_call_formula_ok", "mean"),
        )
        .reset_index()
        .sort_values(["coverage_label", "reference", "tool"])
    )
    return metrics, grouped


def build_metric_deltas(old_collated: Path, new_collated: Path) -> pd.DataFrame:
    old_summary = pd.read_csv(old_collated / "metrics_summary_long.csv")
    new_summary = pd.read_csv(new_collated / "metrics_summary_long.csv")
    key_cols = [col for col in ["coverage_label", "quality_label", "reference", "tool"] if col in old_summary.columns and col in new_summary.columns]
    merged = old_summary.merge(
        new_summary,
        on=key_cols,
        how="outer",
        suffixes=("_old", "_new"),
        indicator=True,
    )

    value_cols = [
        col
        for col in new_summary.columns
        if col in old_summary.columns and col not in key_cols and col != "run_id"
    ]
    numeric_cols = [col for col in value_cols if pd.api.types.is_numeric_dtype(new_summary[col]) or pd.api.types.is_numeric_dtype(old_summary[col])]

    for col in numeric_cols:
        merged[f"{col}_delta"] = pd.to_numeric(merged.get(f"{col}_new"), errors="coerce") - pd.to_numeric(
            merged.get(f"{col}_old"), errors="coerce"
        )

    changed_columns = []
    max_abs_deltas = []
    material_changes = []
    for _, row in merged.iterrows():
        row_changed = []
        row_max_abs = 0.0
        for col in numeric_cols:
            old_value = row.get(f"{col}_old")
            new_value = row.get(f"{col}_new")
            if pd.isna(old_value) and pd.isna(new_value):
                continue
            if pd.isna(old_value) != pd.isna(new_value):
                row_changed.append(col)
                row_max_abs = math.inf
                continue
            if not np.isclose(float(old_value), float(new_value), equal_nan=True, atol=1e-9, rtol=0.0):
                row_changed.append(col)
                row_max_abs = max(row_max_abs, abs(float(new_value) - float(old_value)))
        changed_columns.append("|".join(row_changed))
        max_abs_deltas.append(row_max_abs)
        material_changes.append(bool(row_changed))

    merged["changed_columns"] = changed_columns
    merged["max_abs_delta"] = max_abs_deltas
    merged["material_change"] = material_changes
    merged["parser_reason"] = merged["tool"].map(PARSER_REASON_BY_TOOL).fillna("")
    merged["interval_reason"] = "rrna_interval_only"
    merged["change_reason"] = np.where(
        merged["material_change"],
        np.where(
            merged["parser_reason"] != "",
            merged["parser_reason"] + ";rrna_interval_only",
            "rrna_interval_only",
        ),
        "",
    )
    return merged.sort_values(key_cols).reset_index(drop=True)


def _find_first_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError(f"No candidate paths exist: {paths}")


def build_parser_spotchecks(runs_root: Path, run_glob: str, spotcheck_run: str) -> pd.DataFrame:
    rerun_dir = runs_root / spotcheck_run / Path(run_glob).name
    reference = "16s_88_rrsE"
    parsed_base = rerun_dir / "by_reference" / reference / "01_parsed_raw"
    raw_base = runs_root / spotcheck_run / "modifications"
    records: list[dict[str, object]] = []

    def add_check(tool: str, check_name: str, status: bool, observed: object, expected: object, detail: str) -> None:
        records.append(
            {
                "tool": tool,
                "coverage_label": spotcheck_run.removeprefix("results_"),
                "reference": reference,
                "check_name": check_name,
                "status": bool(status),
                "observed": observed,
                "expected": expected,
                "detail": detail,
            }
        )

    differr = pd.read_csv(parsed_base / "differr.csv")
    differr = differr[differr["replicate"] == "rep1"].copy()
    differr_row = differr[np.isfinite(pd.to_numeric(differr["score"], errors="coerce"))].iloc[0]
    add_check(
        "differr",
        "score_type",
        differr_row["score_type"] == "neglog10_fdr",
        differr_row["score_type"],
        "neglog10_fdr",
        "Parsed DiffErr rows must retain transformed significance semantics.",
    )
    add_check(
        "differr",
        "score_equals_neglog10_fdr",
        np.isclose(float(differr_row["score"]), float(differr_row["g_fdr_neglog10"]), atol=1e-12, rtol=0.0),
        differr_row["score"],
        differr_row["g_fdr_neglog10"],
        "Parsed score should equal the explicit -log10(FDR) column.",
    )

    xpore_raw_path = raw_base / "xpore" / "16s_rep1_diffmod" / "diffmod.table"
    xpore_raw = pd.read_csv(xpore_raw_path)
    pvalue_col = next(col for col in xpore_raw.columns if "pval" in col.lower())
    xpore_raw["bh_fdr_expected"] = _benjamini_hochberg(pd.to_numeric(xpore_raw[pvalue_col], errors="coerce"))
    xpore_raw["parsed_position"] = pd.to_numeric(xpore_raw["position"], errors="coerce") + 2
    xpore_raw = xpore_raw.rename(columns={pvalue_col: "raw_pvalue_expected"})
    xpore = pd.read_csv(parsed_base / "xpore.csv")
    xpore = xpore[xpore["replicate"] == "rep1"].copy()
    xpore_row = xpore.merge(
        xpore_raw[["parsed_position", "raw_pvalue_expected", "bh_fdr_expected"]],
        left_on="position",
        right_on="parsed_position",
        how="inner",
    ).iloc[0]
    add_check(
        "xpore",
        "score_type",
        xpore_row["score_type"] == "fdr",
        xpore_row["score_type"],
        "fdr",
        "xPore score type should be BH-adjusted FDR.",
    )
    add_check(
        "xpore",
        "position_shift",
        float(xpore_row["position"]) == float(xpore_row["parsed_position"]),
        xpore_row["position"],
        xpore_row["parsed_position"],
        "Parsed xPore positions should be raw position +2.",
    )
    add_check(
        "xpore",
        "raw_pvalue_preserved",
        np.isclose(float(xpore_row["pvalue"]), float(xpore_row["raw_pvalue_expected"]), atol=1e-12, rtol=0.0),
        xpore_row["pvalue"],
        xpore_row["raw_pvalue_expected"],
        "Parsed pvalue should retain the native raw p-value.",
    )
    add_check(
        "xpore",
        "bh_fdr_used_as_score",
        np.isclose(float(xpore_row["score"]), float(xpore_row["bh_fdr_expected"]), atol=1e-12, rtol=0.0),
        xpore_row["score"],
        xpore_row["bh_fdr_expected"],
        "Parsed score should equal the per-file BH-adjusted FDR.",
    )

    epinano_raw_dir = raw_base / "epinano" / "16s" / "16s_rep1_epinano_error"
    epinano_raw_path = _find_first_existing(
        [
            epinano_raw_dir / "16s_rep1_sumerr.delta-sum_err.prediction.csv",
            epinano_raw_dir / "16s_rep1_mismatch.delta-mis.prediction.csv",
        ]
    )
    epinano_raw = pd.read_csv(epinano_raw_path)
    score_col = next(col for col in ["delta_sum_err", "delta_mis", "sum_err", "z_scores"] if col in epinano_raw.columns)
    epinano_raw["position"] = epinano_raw["chr_pos"].astype(str).str.split().str[1].astype(int)
    epinano = pd.read_csv(parsed_base / "epinano.csv")
    epinano = epinano[epinano["replicate"] == "rep1"].copy()
    epinano_raw = epinano_raw.rename(columns={score_col: "expected_score"})
    epinano_row = epinano.merge(
        epinano_raw[["position", "expected_score"]],
        on="position",
        how="inner",
    ).iloc[0]
    add_check(
        "epinano",
        "score_type",
        epinano_row["score_type"] == ("error" if score_col != "z_scores" else "zscore"),
        epinano_row["score_type"],
        "error" if score_col != "z_scores" else "zscore",
        "Parsed EpiNano score type should follow the selected score column semantics.",
    )
    add_check(
        "epinano",
        "delta_score_used",
        np.isclose(float(epinano_row["score"]), float(epinano_row["expected_score"]), atol=1e-12, rtol=0.0),
        epinano_row["score"],
        epinano_row["expected_score"],
        f"Parsed EpiNano score should use {score_col}.",
    )

    yanocomp_raw = pd.read_csv(raw_base / "yanocomp" / "16s" / "16s_rep1.bed", sep="\t", header=None)
    yanocomp_raw.columns = [f"c{i}" for i in range(len(yanocomp_raw.columns))]
    yanocomp_raw["position_expected"] = pd.to_numeric(yanocomp_raw["c1"], errors="coerce") + 2
    yanocomp = pd.read_csv(parsed_base / "yanocomp.csv")
    yanocomp = yanocomp[yanocomp["replicate"] == "rep1"].copy()
    yanocomp_row = yanocomp.merge(
        yanocomp_raw[["position_expected"]],
        left_on="position",
        right_on="position_expected",
        how="inner",
    ).iloc[0]
    add_check(
        "yanocomp",
        "position_shift",
        float(yanocomp_row["position"]) == float(yanocomp_row["position_expected"]),
        yanocomp_row["position"],
        yanocomp_row["position_expected"],
        "Parsed Yanocomp positions should be raw start +2.",
    )

    nanocompore_raw = pd.read_csv(
        raw_base / "nanocompore" / "16s_rep1_sampcomp" / "outnanocompore_results.tsv",
        sep="\t",
    )
    genomic_pos = pd.to_numeric(nanocompore_raw.get("genomicPos"), errors="coerce")
    raw_pos = pd.to_numeric(nanocompore_raw.get("pos"), errors="coerce")
    nanocompore_raw["position_expected"] = genomic_pos.where(genomic_pos.notna(), raw_pos) + 2
    nanocompore = pd.read_csv(parsed_base / "nanocompore.csv")
    nanocompore = nanocompore[nanocompore["replicate"] == "rep1"].copy()
    nanocompore_row = nanocompore.merge(
        nanocompore_raw[["position_expected", "GMM_logit_pvalue"]],
        left_on="position",
        right_on="position_expected",
        how="inner",
    ).iloc[0]
    add_check(
        "nanocompore",
        "position_shift",
        float(nanocompore_row["position"]) == float(nanocompore_row["position_expected"]),
        nanocompore_row["position"],
        nanocompore_row["position_expected"],
        "Parsed Nanocompore positions should be centered with +2.",
    )
    add_check(
        "nanocompore",
        "score_source",
        np.isclose(float(nanocompore_row["score"]), float(nanocompore_row["GMM_logit_pvalue"]), atol=1e-12, rtol=0.0),
        nanocompore_row["score"],
        nanocompore_row["GMM_logit_pvalue"],
        "Parsed Nanocompore scores should use GMM_logit_pvalue.",
    )

    drummer = pd.read_csv(parsed_base / "drummer.csv")
    drummer = drummer[drummer["replicate"] == "rep1"].copy()
    drummer_row = drummer.iloc[0]
    drummer_score_col = next(
        col for col in ["max.G_padj", "G_padj", "OR_padj", "pval", "pvalue", "p_value"] if col in drummer.columns and pd.notna(drummer_row.get(col))
    )
    add_check(
        "drummer",
        "score_priority",
        np.isclose(float(drummer_row["score"]), float(drummer_row[drummer_score_col]), atol=1e-12, rtol=0.0),
        drummer_row["score"],
        drummer_row[drummer_score_col],
        f"Parsed DRUMMER score should use {drummer_score_col}.",
    )

    return pd.DataFrame(records)


def _benjamini_hochberg(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=numeric.index, dtype=float)
    valid = numeric.dropna()
    if valid.empty:
        return result
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for idx in range(n - 1, -1, -1):
        rank = idx + 1
        value = ranked[idx] * n / rank
        running = min(running, value)
        adjusted[idx] = min(running, 1.0)
    adjusted_series = pd.Series(adjusted, index=valid.index[order])
    result.loc[adjusted_series.index] = adjusted_series.sort_index()
    return result


def load_evaluation_intervals(runs_root: Path, run_glob: str) -> pd.DataFrame:
    records = []
    for rerun_dir in sorted(runs_root.glob(run_glob)):
        interval_path = rerun_dir / "metadata" / "evaluation_intervals.csv"
        if not interval_path.exists():
            continue
        coverage_label = rerun_dir.parent.name.removeprefix("results_")
        intervals = pd.read_csv(interval_path)
        intervals["coverage_label"] = coverage_label
        records.append(intervals)
    if not records:
        return pd.DataFrame(columns=["reference", "eval_start", "eval_end", "eval_length", "coverage_label"])
    return pd.concat(records, ignore_index=True)


def scan_metric_ready_counts(runs_root: Path, run_glob: str) -> pd.DataFrame:
    records = []
    for rerun_dir in sorted(runs_root.glob(run_glob)):
        coverage_label = rerun_dir.parent.name.removeprefix("results_")
        interval_path = rerun_dir / "metadata" / "evaluation_intervals.csv"
        intervals = pd.read_csv(interval_path).set_index("reference") if interval_path.exists() else pd.DataFrame()
        for metric_ready in sorted(rerun_dir.glob("by_reference/*/03_metric_ready/*.csv")):
            reference = metric_ready.parents[1].name
            tool = metric_ready.stem
            expected = int(intervals.loc[reference, "eval_length"]) if not intervals.empty and reference in intervals.index else np.nan
            replicate_counts = pd.read_csv(metric_ready, usecols=["replicate"]).value_counts(sort=False).reset_index(name="observed_rows")
            for row in replicate_counts.itertuples(index=False):
                records.append(
                    {
                        "coverage_label": coverage_label,
                        "reference": reference,
                        "tool": tool,
                        "replicate": row.replicate,
                        "observed_rows": int(row.observed_rows),
                        "expected_rows": expected,
                        "rows_match_universe": int(row.observed_rows) == expected,
                    }
                )
    return pd.DataFrame(records)


def write_validation_report(
    output_path: Path,
    schema_comparison: pd.DataFrame,
    metrics: pd.DataFrame,
    metrics_valid_fraction: pd.DataFrame,
    metric_deltas: pd.DataFrame,
    parser_spotchecks: pd.DataFrame,
    evaluation_intervals: pd.DataFrame,
    metric_ready_counts: pd.DataFrame,
) -> None:
    tools_present = sorted(metrics["tool"].dropna().astype(str).unique())
    differr_fraction = metrics_valid_fraction.loc[
        metrics_valid_fraction["tool"] == "differr", "metrics_valid_fraction"
    ]
    universe_16s = evaluation_intervals.loc[evaluation_intervals["reference"] == "16s_88_rrsE", "eval_length"]
    universe_23s = evaluation_intervals.loc[evaluation_intervals["reference"] == "23s_78_rrlB", "eval_length"]
    no_call_ok = bool(metrics["no_call_formula_ok"].all())
    metric_ready_ok = bool(metric_ready_counts["rows_match_universe"].all()) if not metric_ready_counts.empty else False
    changed = metric_deltas[metric_deltas["material_change"]]

    lines = [
        "# Validation Report",
        "",
        "## Summary",
        f"- tools_present ({len(tools_present)}): {', '.join(tools_present)}",
        f"- schema_unchanged: {bool(schema_comparison['schema_match'].all())}",
        f"- differr_metrics_valid_fraction_min: {differr_fraction.min() if not differr_fraction.empty else 'NA'}",
        f"- rrna_eval_length_16s_unique: {sorted(universe_16s.dropna().astype(int).unique().tolist())}",
        f"- rrna_eval_length_23s_unique: {sorted(universe_23s.dropna().astype(int).unique().tolist())}",
        f"- no_call_formula_ok_all_rows: {no_call_ok}",
        f"- metric_ready_rows_match_universe: {metric_ready_ok}",
        f"- parser_spotchecks_all_pass: {bool(parser_spotchecks['status'].all())}",
        "",
        "## Schema Comparison",
    ]
    for row in schema_comparison.itertuples(index=False):
        lines.append(f"- {row.file}: schema_match={row.schema_match}")

    lines.extend(
        [
            "",
            "## Parser Spotchecks",
        ]
    )
    for row in parser_spotchecks.itertuples(index=False):
        lines.append(
            f"- {row.tool} {row.check_name}: status={row.status} observed={row.observed} expected={row.expected}"
        )

    lines.extend(
        [
            "",
            "## Material Metric Changes",
        ]
    )
    if changed.empty:
        lines.append("- No material summary-metric deltas detected.")
    else:
        for row in changed.sort_values(["tool", "coverage_label", "reference"]).itertuples(index=False):
            lines.append(
                f"- {row.tool} {row.coverage_label} {row.reference}: reasons={row.change_reason} changed_columns={row.changed_columns}"
            )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    runs_root = Path(args.runs_root)
    old_collated = Path(args.old_collated)
    new_collated = Path(args.new_collated)
    new_collated.mkdir(parents=True, exist_ok=True)

    schema_comparison = compare_schemas(old_collated, new_collated)
    metrics, metrics_valid_fraction = build_metrics_valid_fraction(new_collated)
    metric_deltas = build_metric_deltas(old_collated, new_collated)
    parser_spotchecks = build_parser_spotchecks(runs_root, args.run_glob, args.spotcheck_run)
    evaluation_intervals = load_evaluation_intervals(runs_root, args.run_glob)
    metric_ready_counts = scan_metric_ready_counts(runs_root, args.run_glob)

    metrics_valid_fraction.to_csv(new_collated / "validation_metrics_valid_fraction.csv", index=False)
    metric_deltas.to_csv(new_collated / "validation_metric_deltas.csv", index=False)
    parser_spotchecks.to_csv(new_collated / "validation_parser_spotchecks.csv", index=False)
    write_validation_report(
        new_collated / "validation_report.md",
        schema_comparison,
        metrics,
        metrics_valid_fraction,
        metric_deltas,
        parser_spotchecks,
        evaluation_intervals,
        metric_ready_counts,
    )

    summary = {
        "tools_present": sorted(metrics["tool"].dropna().astype(str).unique().tolist()),
        "schema_unchanged": bool(schema_comparison["schema_match"].all()),
        "parser_spotchecks_all_pass": bool(parser_spotchecks["status"].all()),
        "no_call_formula_ok_all_rows": bool(metrics["no_call_formula_ok"].all()),
        "metric_ready_rows_match_universe": bool(metric_ready_counts["rows_match_universe"].all()) if not metric_ready_counts.empty else False,
    }
    (new_collated / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
