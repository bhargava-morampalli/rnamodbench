#!/usr/bin/env python3
"""
Downstream analysis orchestration.

Key guarantees:
- never mixes references (e.g. 16S and 23S)
- preserves raw tool outputs (writes only to new downstream folder)
- standardizes each tool/replicate to full reference universe (1..L)
- emits per-reference metrics + collation-ready long tables
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve

# allow script-mode local imports
sys.path.insert(0, str(Path(__file__).parent))

from parse_outputs import (
    build_reference_aliases,
    load_all_tool_outputs,
    parse_tool_output,
)
from benchmark_metrics import _metrics_from_arrays, calculate_metrics, load_ground_truth
from data_quality import generate_quality_summary, validate_all_outputs
from position_standardization import (
    build_metric_ready_table,
    canonicalize_tool_references,
    compute_no_call_diagnostics,
    get_reference_eval_regions,
    get_reference_lengths,
    load_reference_catalog,
    standardize_to_reference_universe,
)
from replicate_analysis import calculate_concordance, consensus_calling
from tool_comparison import (
    calculate_pairwise_agreement,
    compare_tools,
    generate_upset_data,
    get_consensus_positions,
)
try:
    from visualization import (
        plot_lag_scan_panels,
        plot_pr_tool_grid,
        plot_roc_tool_grid,
        plot_window_tolerance_panels,
    )
except ImportError:
    plot_lag_scan_panels = None
    plot_pr_tool_grid = None
    plot_roc_tool_grid = None
    plot_window_tolerance_panels = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

SUPPORTED_TOOLS = [
    "tombo",
    "yanocomp",
    "nanocompore",
    "xpore",
    "eligos",
    "epinano",
    "differr",
    "drummer",
    "jacusa2",
    "nanodoc",
]

WINDOW_GRID = list(range(0, 11))
LAG_GRID = list(range(-10, 11))
DIFFERR_PRIMARY_SCORE_FIELD = "g_fdr_neglog10"
DIFFERR_SECONDARY_SCORE_FIELD = "g_stat"
DIFFERR_BOTH_SCORE_FIELD = "both"
DIFFERR_GSTAT_SUBDIR = "differr_gstat"


def load_tool_outputs_from_args(
    args: argparse.Namespace,
    differr_score_field: str = DIFFERR_PRIMARY_SCORE_FIELD,
    reference_aliases: Optional[Dict[str, str]] = None,
) -> Dict[str, pd.DataFrame]:
    tool_outputs: Dict[str, pd.DataFrame] = {}

    for tool in SUPPORTED_TOOLS:
        arg_name = tool.replace("-", "_")
        p = getattr(args, arg_name, None)
        if not p:
            continue

        path = Path(p)
        if not path.exists():
            logger.warning("Tool path not found: %s", path)
            continue

        try:
            df = parse_tool_output(
                tool,
                path,
                differr_score_field=differr_score_field,
                reference_aliases=reference_aliases,
            )
            if not df.empty:
                tool_outputs[tool] = df
                logger.info("Loaded %d rows for %s", len(df), tool)
        except Exception as exc:
            logger.error("Failed loading %s from %s: %s", tool, path, exc)

    return tool_outputs


def _canonicalize_ground_truth(
    gt_df: Optional[pd.DataFrame], catalog: Optional[pd.DataFrame]
) -> Optional[pd.DataFrame]:
    if gt_df is None or gt_df.empty or catalog is None or catalog.empty:
        return gt_df

    alias_map = {}
    for _, row in catalog.iterrows():
        target = str(row["target"]).lower()
        ref_id = str(row["reference_id"])
        alias_map[target] = ref_id
        alias_map[ref_id.lower()] = ref_id

    out = gt_df.copy()

    def canon(ref: str) -> str:
        r = str(ref).strip()
        rl = r.lower()
        if rl in alias_map:
            return alias_map[rl]
        for target, ref_id in alias_map.items():
            if target in {"16s", "23s", "5s"} and rl.startswith(target):
                return ref_id
        return r

    out["reference"] = out["reference"].astype(str).map(canon)
    out["position"] = pd.to_numeric(out["position"], errors="coerce")
    out = out[out["position"].notna()].copy()
    out["position"] = out["position"].astype(int)
    return out


def _infer_reference_lengths(
    tool_outputs: Dict[str, pd.DataFrame], ground_truth: Optional[pd.DataFrame]
) -> Dict[str, int]:
    max_by_ref: Dict[str, int] = {}

    for df in tool_outputs.values():
        if df.empty:
            continue
        tmp = df.copy()
        tmp["position"] = pd.to_numeric(tmp["position"], errors="coerce")
        tmp = tmp[tmp["position"].notna()]
        for ref, grp in tmp.groupby("reference"):
            max_by_ref[str(ref)] = max(max_by_ref.get(str(ref), 0), int(grp["position"].max()))

    if ground_truth is not None and not ground_truth.empty:
        gt = ground_truth.copy()
        gt["position"] = pd.to_numeric(gt["position"], errors="coerce")
        gt = gt[gt["position"].notna()]
        for ref, grp in gt.groupby("reference"):
            max_by_ref[str(ref)] = max(max_by_ref.get(str(ref), 0), int(grp["position"].max()))

    return {k: int(v) for k, v in max_by_ref.items() if v > 0}


def _make_dirs(base: Path) -> Dict[str, Path]:
    dirs = {
        "parsed": base / "01_parsed_raw",
        "imputed": base / "02_imputed_nan",
        "metric_ready": base / "03_metric_ready",
        "metrics": base / "04_metrics",
        "curves": base / "05_curves",
        "comparison": base / "06_tool_comparison",
        "replicate": base / "07_replicate_analysis",
        "visualizations": base / "08_visualizations",
        "tolerance": base / "09_tolerance_analysis",
        "lag": base / "10_lag_diagnostics",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def _extract_curve_points(
    metric_ready_df: pd.DataFrame,
    tool: str,
    reference: str,
    replicate: str,
) -> Dict[str, pd.DataFrame]:
    rows_pr = []
    rows_roc = []

    sub = metric_ready_df[
        (metric_ready_df["tool"] == tool)
        & (metric_ready_df["reference"] == reference)
        & (metric_ready_df["replicate"].astype(str) == str(replicate))
    ].copy()

    if sub.empty:
        return {
            "pr": pd.DataFrame(rows_pr),
            "roc": pd.DataFrame(rows_roc),
        }

    def _append_scope(scope_df: pd.DataFrame, scope: str) -> None:
        if scope_df.empty:
            return

        y_true = pd.to_numeric(scope_df["label"], errors="coerce").fillna(0).astype(int)
        y_score = pd.to_numeric(scope_df["score_metric"], errors="coerce")
        valid = y_score.notna() & np.isfinite(y_score.to_numpy())
        if int(valid.sum()) == 0:
            return

        y_true = y_true.loc[valid].to_numpy()
        y_score = y_score.loc[valid].astype(float).to_numpy()

        if len(np.unique(y_true)) < 2:
            return

        precision, recall, pr_thresholds = precision_recall_curve(y_true, y_score)
        for i, (p, r) in enumerate(zip(precision, recall)):
            rows_pr.append(
                {
                    "reference": reference,
                    "tool": tool,
                    "replicate": replicate,
                    "scope": scope,
                    "index": i,
                    "precision": float(p),
                    "recall": float(r),
                    "threshold": float(pr_thresholds[i]) if i < len(pr_thresholds) else np.nan,
                }
            )

        fpr, tpr, roc_thresholds = roc_curve(y_true, y_score)
        for i, (f, t, th) in enumerate(zip(fpr, tpr, roc_thresholds)):
            rows_roc.append(
                {
                    "reference": reference,
                    "tool": tool,
                    "replicate": replicate,
                    "scope": scope,
                    "index": i,
                    "fpr": float(f),
                    "tpr": float(t),
                    "threshold": float(th),
                }
            )

    _append_scope(sub, "universe")
    if "is_reported" in sub.columns:
        _append_scope(sub[sub["is_reported"] == True].copy(), "reported")

    return {"pr": pd.DataFrame(rows_pr), "roc": pd.DataFrame(rows_roc)}


def _compute_subset_metrics(sub_df: pd.DataFrame) -> Optional[Dict[str, float]]:
    if sub_df.empty:
        return None

    y_true = pd.to_numeric(sub_df.get("label", np.nan), errors="coerce")
    y_score = pd.to_numeric(sub_df.get("score_metric", np.nan), errors="coerce")
    valid = y_true.notna() & y_score.notna() & np.isfinite(y_score.to_numpy())
    if int(valid.sum()) == 0:
        return None

    y_true_arr = y_true[valid].astype(int).to_numpy()
    y_score_arr = y_score[valid].astype(float).to_numpy()

    if len(np.unique(y_true_arr)) < 2:
        return None

    (
        auprc,
        auroc,
        f1_optimal,
        _precision_optimal,
        _recall_optimal,
        _optimal_threshold,
        _tp,
        _fp,
        _fn,
        _tn,
    ) = _metrics_from_arrays(y_true_arr, y_score_arr)

    return {
        "auprc_reported": float(auprc),
        "auroc_reported": float(auroc),
        "f1_reported": float(f1_optimal),
    }


def _invalid_score_stats(rep_df_reported: pd.DataFrame) -> tuple[int, int, float]:
    if rep_df_reported.empty:
        return 0, 0, 0.0

    score_type = "unknown"
    if "score_type" in rep_df_reported.columns and rep_df_reported["score_type"].notna().any():
        score_type = str(rep_df_reported["score_type"].dropna().iloc[0]).strip().lower()

    if score_type not in {"pvalue", "fdr"}:
        return 0, 0, 0.0

    if "score_raw" not in rep_df_reported.columns:
        return 0, 0, 0.0

    raw = pd.to_numeric(rep_df_reported["score_raw"], errors="coerce")
    invalid_mask = (~np.isfinite(raw.to_numpy())) | (raw.to_numpy() < 0) | (raw.to_numpy() > 1)

    invalid_n = int(invalid_mask.sum())
    invalid_total = int(len(raw))
    invalid_fraction = float(invalid_n / invalid_total) if invalid_total else 0.0
    return invalid_n, invalid_total, invalid_fraction


def _expand_positions_by_window(
    ground_truth_positions: Set[int],
    window_size: int,
    eval_start: int,
    eval_end: Optional[int] = None,
) -> Set[int]:
    if not ground_truth_positions:
        return set()

    if eval_end is None:
        eval_end = int(eval_start)
        eval_start = 1

    expanded: Set[int] = set()
    for pos in ground_truth_positions:
        p = int(pos)
        for offset in range(-int(window_size), int(window_size) + 1):
            q = p + offset
            if int(eval_start) <= q <= int(eval_end):
                expanded.add(q)
    return expanded


def _shift_positions(
    ground_truth_positions: Set[int],
    delta: int,
    eval_start: int,
    eval_end: Optional[int] = None,
) -> Set[int]:
    if not ground_truth_positions:
        return set()

    if eval_end is None:
        eval_end = int(eval_start)
        eval_start = 1

    shifted: Set[int] = set()
    for pos in ground_truth_positions:
        q = int(pos) + int(delta)
        if int(eval_start) <= q <= int(eval_end):
            shifted.add(q)
    return shifted


def _compute_labeled_scope_metrics(
    scope_df: pd.DataFrame,
    positive_positions: Set[int],
) -> Dict[str, float]:
    n_universe = int(len(scope_df))
    if n_universe == 0:
        return {
            "auprc": np.nan,
            "auroc": np.nan,
            "f1_optimal": np.nan,
            "precision_optimal": np.nan,
            "recall_optimal": np.nan,
            "n_true_positives": np.nan,
            "n_false_positives": np.nan,
            "n_false_negatives": np.nan,
            "n_true_negatives": np.nan,
            "optimal_threshold": np.nan,
            "n_positive": 0,
            "n_universe": 0,
            "prevalence": np.nan,
            "metric_scope_note": "empty_scope",
        }

    scores = pd.to_numeric(scope_df.get("score_metric", np.nan), errors="coerce")
    valid = scores.notna() & np.isfinite(scores.to_numpy())
    if int(valid.sum()) == 0:
        return {
            "auprc": np.nan,
            "auroc": np.nan,
            "f1_optimal": np.nan,
            "precision_optimal": np.nan,
            "recall_optimal": np.nan,
            "n_true_positives": np.nan,
            "n_false_positives": np.nan,
            "n_false_negatives": np.nan,
            "n_true_negatives": np.nan,
            "optimal_threshold": np.nan,
            "n_positive": 0,
            "n_universe": n_universe,
            "prevalence": 0.0,
            "metric_scope_note": "no_finite_scores",
        }

    sub = scope_df.loc[valid].copy()
    y_scores = scores.loc[valid].astype(float).to_numpy()
    y_true = sub["position"].astype(int).isin(positive_positions).astype(int).to_numpy()
    n_positive = int(y_true.sum())
    n_total = int(len(y_true))
    prevalence = float(n_positive / n_total) if n_total else np.nan

    if n_total == 0:
        return {
            "auprc": np.nan,
            "auroc": np.nan,
            "f1_optimal": np.nan,
            "precision_optimal": np.nan,
            "recall_optimal": np.nan,
            "n_true_positives": np.nan,
            "n_false_positives": np.nan,
            "n_false_negatives": np.nan,
            "n_true_negatives": np.nan,
            "optimal_threshold": np.nan,
            "n_positive": 0,
            "n_universe": 0,
            "prevalence": np.nan,
            "metric_scope_note": "empty_scope",
        }

    if len(np.unique(y_true)) < 2:
        return {
            "auprc": np.nan,
            "auroc": np.nan,
            "f1_optimal": np.nan,
            "precision_optimal": np.nan,
            "recall_optimal": np.nan,
            "n_true_positives": np.nan,
            "n_false_positives": np.nan,
            "n_false_negatives": np.nan,
            "n_true_negatives": np.nan,
            "optimal_threshold": np.nan,
            "n_positive": n_positive,
            "n_universe": n_total,
            "prevalence": prevalence,
            "metric_scope_note": "single_class_labels",
        }

    (
        auprc,
        auroc,
        f1_optimal,
        precision_optimal,
        recall_optimal,
        optimal_threshold,
        tp,
        fp,
        fn,
        tn,
    ) = _metrics_from_arrays(y_true, y_scores)

    return {
        "auprc": float(auprc),
        "auroc": float(auroc),
        "f1_optimal": float(f1_optimal),
        "precision_optimal": float(precision_optimal),
        "recall_optimal": float(recall_optimal),
        "n_true_positives": int(tp),
        "n_false_positives": int(fp),
        "n_false_negatives": int(fn),
        "n_true_negatives": int(tn),
        "optimal_threshold": float(optimal_threshold),
        "n_positive": n_positive,
        "n_universe": n_total,
        "prevalence": prevalence,
        "metric_scope_note": "ok",
    }


def _summarize_sweep_metrics(df: pd.DataFrame, sweep_col: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()

    return (
        df.groupby(["reference", "tool", "scope", sweep_col], as_index=False)
        .agg(
            auprc_mean=("auprc", "mean"),
            auprc_std=("auprc", "std"),
            auprc_median=("auprc", "median"),
            auroc_mean=("auroc", "mean"),
            auroc_std=("auroc", "std"),
            auroc_median=("auroc", "median"),
            f1_optimal_mean=("f1_optimal", "mean"),
            f1_optimal_std=("f1_optimal", "std"),
            f1_optimal_median=("f1_optimal", "median"),
            precision_optimal_mean=("precision_optimal", "mean"),
            recall_optimal_mean=("recall_optimal", "mean"),
            prevalence_mean=("prevalence", "mean"),
            n_positive_mean=("n_positive", "mean"),
            n_universe_mean=("n_universe", "mean"),
            metrics_valid_fraction=("metrics_valid", "mean"),
        )
    )


def run_analysis(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: Optional[pd.DataFrame],
    output_dir: Path,
    references_csv: Optional[Path] = None,
    references_filter: Optional[List[str]] = None,
    min_replicates: int = 2,
    score_threshold: float = None,
    expected_replicates: int = 3,
    run_id: Optional[str] = None,
    coverage_label: Optional[str] = None,
    quality_label: Optional[str] = None,
    differr_score_field: str = "g_fdr_neglog10",
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_dir = output_dir / "metadata"
    by_ref_dir = output_dir / "by_reference"
    collation_dir = output_dir / "collation"
    metadata_dir.mkdir(exist_ok=True)
    by_ref_dir.mkdir(exist_ok=True)
    collation_dir.mkdir(exist_ok=True)

    # Data quality checks
    quality_report = validate_all_outputs(tool_outputs, expected_replicates=expected_replicates)
    quality_report.save(metadata_dir / "data_quality_report.csv")
    quality_report.save_warnings(metadata_dir / "analysis_warnings.txt")
    generate_quality_summary(quality_report).to_csv(metadata_dir / "quality_summary.csv", index=False)
    with open(metadata_dir / "quality_summary.txt", "w", encoding="utf-8") as handle:
        handle.write(quality_report.summary())

    # Drop hard-failed tools only
    valid_tools: Dict[str, pd.DataFrame] = {}
    for tool_name, tool_df in tool_outputs.items():
        status_rows = [s for s in quality_report.tool_statuses if s.tool == tool_name]
        if status_rows and all(s.status.value == "FAILED" for s in status_rows):
            logger.warning("Excluding %s due to FAILED status", tool_name)
            continue
        valid_tools[tool_name] = tool_df

    tool_outputs = valid_tools
    if not tool_outputs:
        logger.error("No valid tool outputs after quality filtering")
        return

    # Reference catalog and canonicalization
    reference_catalog = None
    reference_lengths: Dict[str, int] = {}
    evaluation_regions: Dict[str, Tuple[int, int]] = {}

    if references_csv is not None:
        reference_catalog = load_reference_catalog(Path(references_csv))
        reference_catalog.to_csv(metadata_dir / "reference_catalog.csv", index=False)
        reference_lengths = get_reference_lengths(reference_catalog)
        evaluation_regions = get_reference_eval_regions(reference_catalog)

        for tool, df in tool_outputs.items():
            if not df.empty:
                tool_outputs[tool] = canonicalize_tool_references(df, reference_catalog)

        ground_truth = _canonicalize_ground_truth(ground_truth, reference_catalog)
    else:
        logger.warning("No --references-csv provided. Reference lengths inferred from observed max positions.")

    if not reference_lengths:
        reference_lengths = _infer_reference_lengths(tool_outputs, ground_truth)
    if not evaluation_regions:
        evaluation_regions = {ref: (1, int(length)) for ref, length in reference_lengths.items()}

    # Decide references to process
    refs_from_data: Set[str] = set()
    for df in tool_outputs.values():
        if not df.empty:
            refs_from_data.update(df["reference"].dropna().astype(str).tolist())

    if ground_truth is not None and not ground_truth.empty:
        refs_from_data.update(ground_truth["reference"].dropna().astype(str).tolist())

    refs = sorted(set(reference_lengths.keys()) | refs_from_data)

    if references_filter:
        wanted = {r.lower() for r in references_filter}
        if reference_catalog is not None and not reference_catalog.empty:
            for _, row in reference_catalog.iterrows():
                target = str(row["target"]).lower()
                ref_id = str(row["reference_id"]).lower()
                if target in wanted:
                    wanted.add(ref_id)
        refs = [r for r in refs if r.lower() in wanted]

    if not refs:
        logger.error("No references detected for downstream analysis")
        return

    all_metrics_rows: List[pd.DataFrame] = []
    all_summary_rows: List[pd.DataFrame] = []
    all_window_rows: List[pd.DataFrame] = []
    all_window_summary_rows: List[pd.DataFrame] = []
    all_lag_rows: List[pd.DataFrame] = []
    all_lag_summary_rows: List[pd.DataFrame] = []

    for ref in refs:
        if ref not in reference_lengths:
            logger.warning("Skipping reference %s: no known length", ref)
            continue

        ref_len = int(reference_lengths[ref])
        eval_start, eval_end = evaluation_regions.get(ref, (1, ref_len))
        eval_positions = list(range(int(eval_start), int(eval_end) + 1))
        logger.info(
            "Processing reference %s (length=%d, eval_start=%d, eval_end=%d)",
            ref,
            ref_len,
            eval_start,
            eval_end,
        )

        ref_root = by_ref_dir / ref
        dirs = _make_dirs(ref_root)

        gt_ref = (
            ground_truth[ground_truth["reference"].astype(str) == ref].copy()
            if ground_truth is not None and not ground_truth.empty
            else pd.DataFrame(columns=["reference", "position"])
        )
        gt_positions = (
            {
                int(pos)
                for pos in gt_ref["position"].astype(int).tolist()
                if int(eval_start) <= int(pos) <= int(eval_end)
            }
            if not gt_ref.empty
            else set()
        )

        # Keep per-reference tool slices
        parsed_ref_tools: Dict[str, pd.DataFrame] = {}
        metric_ready_by_tool: Dict[str, pd.DataFrame] = {}
        replicate_gate: Dict[Tuple[str, str], Dict[str, float]] = {}

        # Standardization report for this reference
        from position_standardization import StandardizationReport

        std_report = StandardizationReport()

        metric_rows = []
        no_call_rows = []
        pr_rows = []
        roc_rows = []

        for tool, df in tool_outputs.items():
            if df.empty:
                continue

            tool_ref = df[df["reference"].astype(str) == ref].copy()
            if tool_ref.empty:
                continue

            parsed_ref_tools[tool] = tool_ref
            tool_ref.to_csv(dirs["parsed"] / f"{tool}.csv", index=False)

            reps = sorted(tool_ref["replicate"].astype(str).unique().tolist()) if "replicate" in tool_ref.columns else ["rep1"]

            imputed = standardize_to_reference_universe(
                tool_ref,
                reference_id=ref,
                reference_length=ref_len,
                universe_positions=eval_positions,
                replicates=reps,
                report=std_report,
            )
            imputed.to_csv(dirs["imputed"] / f"{tool}.csv", index=False)

            metric_ready = build_metric_ready_table(imputed, gt_positions)
            metric_ready.to_csv(dirs["metric_ready"] / f"{tool}.csv", index=False)
            metric_ready_by_tool[tool] = metric_ready

            # No-call diagnostics
            no_call = compute_no_call_diagnostics(metric_ready, gt_positions)
            no_call_rows.append(no_call)

            # Metrics per replicate
            for rep in sorted(metric_ready["replicate"].astype(str).unique().tolist()):
                rep_df = metric_ready[metric_ready["replicate"].astype(str) == rep].copy()
                if rep_df.empty:
                    continue

                rep_df_reported = (
                    rep_df[rep_df["is_reported"] == True].copy()
                    if "is_reported" in rep_df.columns
                    else rep_df.copy()
                )

                invalid_n, invalid_total, invalid_fraction = _invalid_score_stats(rep_df_reported)
                metrics_valid = not (invalid_total > 0 and invalid_fraction > 0.5)
                metric_scope_note = "ok"
                replicate_gate[(tool, str(rep))] = {
                    "metrics_valid": bool(metrics_valid),
                    "invalid_score_fraction": float(invalid_fraction),
                    "invalid_score_n": int(invalid_n),
                    "invalid_score_total": int(invalid_total),
                }

                m = calculate_metrics(rep_df, ground_truth=None, reference=ref)
                metric_row = m.to_dict()

                # reported-only metrics complement full-universe metrics
                reported_metrics = _compute_subset_metrics(rep_df_reported)
                if reported_metrics is None:
                    metric_row["auprc_reported"] = np.nan
                    metric_row["auroc_reported"] = np.nan
                    metric_row["f1_reported"] = np.nan
                    if metrics_valid:
                        metric_scope_note = "reported_single_class"
                else:
                    metric_row.update(reported_metrics)

                if not metrics_valid:
                    metric_scope_note = "invalid_scores"

                metric_row["metrics_valid"] = bool(metrics_valid)
                metric_row["invalid_score_fraction"] = invalid_fraction
                metric_row["invalid_score_n"] = invalid_n
                metric_row["invalid_score_total"] = invalid_total
                metric_row["metric_scope_note"] = metric_scope_note
                metric_rows.append(pd.DataFrame([metric_row]))

                curves = _extract_curve_points(metric_ready, tool=tool, reference=ref, replicate=rep)
                if not curves["pr"].empty:
                    pr_rows.append(curves["pr"])
                if not curves["roc"].empty:
                    roc_rows.append(curves["roc"])

        # Save standardization reports
        std_report.save(ref_root)

        # Metrics output files
        metrics_per_rep = pd.concat(metric_rows, ignore_index=True) if metric_rows else pd.DataFrame()
        no_call_df = pd.concat(no_call_rows, ignore_index=True) if no_call_rows else pd.DataFrame()

        if not metrics_per_rep.empty and not no_call_df.empty:
            metrics_per_rep = metrics_per_rep.merge(
                no_call_df[
                    [
                        "reference",
                        "tool",
                        "replicate",
                        "n_universe",
                        "n_reported",
                        "n_no_call",
                        "no_call_rate",
                        "n_gt_total",
                        "n_gt_reported",
                        "gt_recall_raw",
                    ]
                ],
                on=["reference", "tool", "replicate"],
                how="left",
            )

        if not metrics_per_rep.empty:
            summary = (
                metrics_per_rep.groupby(["reference", "tool"], as_index=False)
                .agg(
                    auprc_mean=("auprc", "mean"),
                    auprc_std=("auprc", "std"),
                    auprc_median=("auprc", "median"),
                    auroc_mean=("auroc", "mean"),
                    auroc_std=("auroc", "std"),
                    auroc_median=("auroc", "median"),
                    auprc_reported_mean=("auprc_reported", "mean"),
                    auprc_reported_std=("auprc_reported", "std"),
                    auprc_reported_median=("auprc_reported", "median"),
                    auroc_reported_mean=("auroc_reported", "mean"),
                    auroc_reported_std=("auroc_reported", "std"),
                    auroc_reported_median=("auroc_reported", "median"),
                    f1_optimal_mean=("f1_optimal", "mean"),
                    f1_optimal_std=("f1_optimal", "std"),
                    f1_optimal_median=("f1_optimal", "median"),
                    f1_reported_mean=("f1_reported", "mean"),
                    f1_reported_std=("f1_reported", "std"),
                    f1_reported_median=("f1_reported", "median"),
                    precision_optimal_mean=("precision_optimal", "mean"),
                    recall_optimal_mean=("recall_optimal", "mean"),
                    metrics_valid_fraction=("metrics_valid", "mean"),
                    no_call_rate_mean=("no_call_rate", "mean"),
                    gt_recall_raw_mean=("gt_recall_raw", "mean"),
                )
            )
            all_metrics_rows.append(metrics_per_rep)
            all_summary_rows.append(summary)
        else:
            summary = pd.DataFrame()

        metrics_per_rep.to_csv(dirs["metrics"] / "metrics_per_replicate.csv", index=False)
        summary.to_csv(dirs["metrics"] / "metrics_summary.csv", index=False)

        no_call_df.to_csv(dirs["metrics"] / "no_call_diagnostics.csv", index=False)

        # Always-on tolerance and lag analyses (additive; exact outputs unchanged)
        window_rows_ref: List[pd.DataFrame] = []
        lag_rows_ref: List[pd.DataFrame] = []
        for tool, metric_ready in metric_ready_by_tool.items():
            reps = sorted(metric_ready["replicate"].astype(str).unique().tolist())
            for rep in reps:
                rep_df = metric_ready[metric_ready["replicate"].astype(str) == str(rep)].copy()
                rep_reported = rep_df[rep_df["is_reported"] == True].copy() if "is_reported" in rep_df.columns else rep_df.copy()
                gate = replicate_gate.get(
                    (tool, str(rep)),
                    {
                        "metrics_valid": True,
                        "invalid_score_fraction": 0.0,
                        "invalid_score_n": 0,
                        "invalid_score_total": 0,
                    },
                )
                for scope_name, scope_df in [("universe", rep_df), ("reported", rep_reported)]:
                    for window_size in WINDOW_GRID:
                        positives = _expand_positions_by_window(
                            gt_positions, window_size, eval_start, eval_end
                        )
                        row = _compute_labeled_scope_metrics(scope_df, positives)
                        row.update(
                            {
                                "reference": ref,
                                "tool": tool,
                                "replicate": str(rep),
                                "scope": scope_name,
                                "window_size": int(window_size),
                                "metrics_valid": bool(gate["metrics_valid"]),
                                "invalid_score_fraction": float(gate["invalid_score_fraction"]),
                                "invalid_score_n": int(gate["invalid_score_n"]),
                                "invalid_score_total": int(gate["invalid_score_total"]),
                            }
                        )
                        if not bool(gate["metrics_valid"]):
                            row["metric_scope_note"] = "invalid_scores"
                        window_rows_ref.append(pd.DataFrame([row]))

                    for delta in LAG_GRID:
                        positives = _shift_positions(gt_positions, delta, eval_start, eval_end)
                        row = _compute_labeled_scope_metrics(scope_df, positives)
                        row.update(
                            {
                                "reference": ref,
                                "tool": tool,
                                "replicate": str(rep),
                                "scope": scope_name,
                                "delta": int(delta),
                                "metrics_valid": bool(gate["metrics_valid"]),
                                "invalid_score_fraction": float(gate["invalid_score_fraction"]),
                                "invalid_score_n": int(gate["invalid_score_n"]),
                                "invalid_score_total": int(gate["invalid_score_total"]),
                            }
                        )
                        if not bool(gate["metrics_valid"]):
                            row["metric_scope_note"] = "invalid_scores"
                        lag_rows_ref.append(pd.DataFrame([row]))

        window_metrics = pd.concat(window_rows_ref, ignore_index=True) if window_rows_ref else pd.DataFrame()
        lag_metrics = pd.concat(lag_rows_ref, ignore_index=True) if lag_rows_ref else pd.DataFrame()

        window_summary = _summarize_sweep_metrics(window_metrics, "window_size")
        lag_summary = _summarize_sweep_metrics(lag_metrics, "delta")

        best_offset = pd.DataFrame()
        if not lag_summary.empty:
            tmp = lag_summary.copy()
            auroc_best = (
                tmp.sort_values("auroc_mean", ascending=False)
                .groupby(["reference", "tool", "scope"], as_index=False)
                .first()[["reference", "tool", "scope", "delta", "auroc_mean"]]
                .rename(columns={"delta": "best_delta_auroc", "auroc_mean": "best_auroc_mean"})
            )
            auprc_best = (
                tmp.sort_values("auprc_mean", ascending=False)
                .groupby(["reference", "tool", "scope"], as_index=False)
                .first()[["reference", "tool", "scope", "delta", "auprc_mean"]]
                .rename(columns={"delta": "best_delta_auprc", "auprc_mean": "best_auprc_mean"})
            )
            best_offset = auroc_best.merge(
                auprc_best, on=["reference", "tool", "scope"], how="outer"
            )

        window_metrics.to_csv(dirs["tolerance"] / "window_metrics_per_replicate.csv", index=False)
        window_summary.to_csv(dirs["tolerance"] / "window_metrics_summary.csv", index=False)
        lag_metrics.to_csv(dirs["lag"] / "lag_metrics_per_replicate.csv", index=False)
        lag_summary.to_csv(dirs["lag"] / "lag_metrics_summary.csv", index=False)
        best_offset.to_csv(dirs["lag"] / "lag_best_offset_by_tool.csv", index=False)

        if not window_metrics.empty:
            all_window_rows.append(window_metrics)
        if not window_summary.empty:
            all_window_summary_rows.append(window_summary)
        if not lag_metrics.empty:
            all_lag_rows.append(lag_metrics)
        if not lag_summary.empty:
            all_lag_summary_rows.append(lag_summary)

        pr_df = pd.concat(pr_rows, ignore_index=True) if pr_rows else pd.DataFrame()
        roc_df = pd.concat(roc_rows, ignore_index=True) if roc_rows else pd.DataFrame()

        pr_df.to_csv(dirs["curves"] / "pr_curve_points.csv", index=False)
        roc_df.to_csv(dirs["curves"] / "roc_curve_points.csv", index=False)

        # Publication-style per-reference visualizations (non-fatal)
        try:
            if not metrics_per_rep.empty:
                pr_universe = pr_df[pr_df["scope"] == "universe"].copy() if "scope" in pr_df.columns else pr_df
                roc_universe = roc_df[roc_df["scope"] == "universe"].copy() if "scope" in roc_df.columns else roc_df
                if plot_roc_tool_grid is not None:
                    plot_roc_tool_grid(
                        curve_df=roc_universe,
                        metrics_df=metrics_per_rep,
                        reference=ref,
                        output_dir=dirs["visualizations"],
                        ncols=3,
                        output_stem="roc_tool_grid",
                    )
                if plot_pr_tool_grid is not None:
                    plot_pr_tool_grid(
                        curve_df=pr_universe,
                        metrics_df=metrics_per_rep,
                        reference=ref,
                        output_dir=dirs["visualizations"],
                        ncols=3,
                        output_stem="pr_tool_grid",
                    )

                # reported-only visualization using reported-only metrics and prevalence context
                metrics_reported = metrics_per_rep.copy()
                metrics_reported["auroc"] = metrics_reported["auroc_reported"]
                metrics_reported["auprc"] = metrics_reported["auprc_reported"]
                metrics_reported["n_gt_total"] = metrics_reported["n_gt_reported"]
                metrics_reported["n_universe"] = metrics_reported["n_reported"]

                pr_reported = pr_df[pr_df["scope"] == "reported"].copy() if "scope" in pr_df.columns else pd.DataFrame()
                roc_reported = roc_df[roc_df["scope"] == "reported"].copy() if "scope" in roc_df.columns else pd.DataFrame()

                if plot_roc_tool_grid is not None:
                    plot_roc_tool_grid(
                        curve_df=roc_reported,
                        metrics_df=metrics_reported,
                        reference=ref,
                        output_dir=dirs["visualizations"],
                        ncols=3,
                        output_stem="roc_tool_grid_reported",
                    )
                if plot_pr_tool_grid is not None:
                    plot_pr_tool_grid(
                        curve_df=pr_reported,
                        metrics_df=metrics_reported,
                        reference=ref,
                        output_dir=dirs["visualizations"],
                        ncols=3,
                        output_stem="pr_tool_grid_reported",
                    )

                # always-on tolerance and lag panel figures (AUROC + AUPRC)
                for scope in ("universe", "reported"):
                    scope_window = window_metrics[window_metrics["scope"] == scope].copy() if not window_metrics.empty else pd.DataFrame()
                    scope_lag = lag_metrics[lag_metrics["scope"] == scope].copy() if not lag_metrics.empty else pd.DataFrame()

                    if not scope_window.empty and plot_window_tolerance_panels is not None:
                        plot_window_tolerance_panels(
                            sweep_df=scope_window,
                            reference=ref,
                            output_dir=dirs["visualizations"],
                            output_stem=f"window_tolerance_{scope}",
                        )
                    if not scope_lag.empty and plot_lag_scan_panels is not None:
                        plot_lag_scan_panels(
                            sweep_df=scope_lag,
                            reference=ref,
                            output_dir=dirs["visualizations"],
                            output_stem=f"lag_scan_{scope}",
                        )
        except Exception as exc:
            logger.warning("Visualization generation failed for reference %s: %s", ref, exc)

        # Tool comparison (reference-specific only)
        if parsed_ref_tools:
            cmp = compare_tools(parsed_ref_tools, threshold=score_threshold, reference=ref)
            pd.DataFrame([cmp.to_dict()]).to_csv(dirs["comparison"] / "summary.csv", index=False)

            upset = generate_upset_data(parsed_ref_tools, threshold=score_threshold, reference=ref)
            if not upset.empty:
                upset.to_csv(dirs["comparison"] / "upset_data.csv", index=False)

            pairwise = calculate_pairwise_agreement(parsed_ref_tools, threshold=score_threshold, reference=ref)
            if not pairwise.empty:
                pairwise.to_csv(dirs["comparison"] / "pairwise_agreement.csv", index=False)

            for min_tools in [2, 3, len(parsed_ref_tools)]:
                if min_tools <= len(parsed_ref_tools):
                    cons = get_consensus_positions(
                        parsed_ref_tools,
                        min_tools=min_tools,
                        threshold=score_threshold,
                        reference=ref,
                    )
                    if not cons.empty:
                        cons.to_csv(dirs["comparison"] / f"consensus_{min_tools}tools.csv", index=False)

        # Replicate analysis (reference-specific only)
        rep_stats_rows = []
        for tool, df in parsed_ref_tools.items():
            reps = sorted(df["replicate"].astype(str).unique()) if "replicate" in df.columns else ["rep1"]
            if len(reps) < 2:
                continue

            stats_obj = calculate_concordance(df, threshold=score_threshold, reference=ref)
            rep_stats_rows.append(pd.DataFrame([stats_obj.to_dict()]))

            consensus = consensus_calling(
                df,
                min_replicates=min_replicates,
                threshold=score_threshold,
                reference=ref,
            )
            if not consensus.empty:
                consensus.to_csv(dirs["replicate"] / f"{tool}_consensus.csv", index=False)

        if rep_stats_rows:
            pd.concat(rep_stats_rows, ignore_index=True).to_csv(
                dirs["replicate"] / "concordance.csv", index=False
            )

    # Collation-ready long tables
    if all_metrics_rows:
        metrics_long = pd.concat(all_metrics_rows, ignore_index=True)
        metrics_long.insert(0, "run_id", run_id or output_dir.name)
        metrics_long.insert(1, "coverage_label", coverage_label if coverage_label is not None else "NA")
        metrics_long.insert(2, "quality_label", quality_label if quality_label is not None else "NA")
        metrics_long.to_csv(collation_dir / "metrics_long.csv", index=False)

    if all_summary_rows:
        summary_long = pd.concat(all_summary_rows, ignore_index=True)
        summary_long.insert(0, "run_id", run_id or output_dir.name)
        summary_long.insert(1, "coverage_label", coverage_label if coverage_label is not None else "NA")
        summary_long.insert(2, "quality_label", quality_label if quality_label is not None else "NA")
        summary_long.to_csv(collation_dir / "metrics_summary_long.csv", index=False)

    if all_window_rows:
        window_long = pd.concat(all_window_rows, ignore_index=True)
        window_long.insert(0, "run_id", run_id or output_dir.name)
        window_long.insert(1, "coverage_label", coverage_label if coverage_label is not None else "NA")
        window_long.insert(2, "quality_label", quality_label if quality_label is not None else "NA")
        window_long.to_csv(collation_dir / "window_metrics_long.csv", index=False)

    if all_window_summary_rows:
        window_summary_long = pd.concat(all_window_summary_rows, ignore_index=True)
        window_summary_long.insert(0, "run_id", run_id or output_dir.name)
        window_summary_long.insert(1, "coverage_label", coverage_label if coverage_label is not None else "NA")
        window_summary_long.insert(2, "quality_label", quality_label if quality_label is not None else "NA")
        window_summary_long.to_csv(collation_dir / "window_metrics_summary_long.csv", index=False)

    if all_lag_rows:
        lag_long = pd.concat(all_lag_rows, ignore_index=True)
        lag_long.insert(0, "run_id", run_id or output_dir.name)
        lag_long.insert(1, "coverage_label", coverage_label if coverage_label is not None else "NA")
        lag_long.insert(2, "quality_label", quality_label if quality_label is not None else "NA")
        lag_long.to_csv(collation_dir / "lag_metrics_long.csv", index=False)

    if all_lag_summary_rows:
        lag_summary_long = pd.concat(all_lag_summary_rows, ignore_index=True)
        lag_summary_long.insert(0, "run_id", run_id or output_dir.name)
        lag_summary_long.insert(1, "coverage_label", coverage_label if coverage_label is not None else "NA")
        lag_summary_long.insert(2, "quality_label", quality_label if quality_label is not None else "NA")
        lag_summary_long.to_csv(collation_dir / "lag_metrics_summary_long.csv", index=False)

    evaluation_intervals_df = pd.DataFrame(
        [
            {
                "reference": ref,
                "full_reference_length": int(reference_lengths[ref]),
                "eval_start": int(start),
                "eval_end": int(end),
                "eval_length": int(end - start + 1),
            }
            for ref, (start, end) in sorted(evaluation_regions.items())
            if ref in reference_lengths
        ]
    )
    evaluation_intervals_df.to_csv(metadata_dir / "evaluation_intervals.csv", index=False)

    run_meta = {
        "run_id": run_id or output_dir.name,
        "coverage_label": coverage_label,
        "quality_label": quality_label,
        "references_processed": refs,
        "reference_lengths": reference_lengths,
        "evaluation_intervals": {
            ref: {"eval_start": int(a), "eval_end": int(b)}
            for ref, (a, b) in evaluation_regions.items()
        },
        "evaluation_lengths": {
            ref: int(b - a + 1) for ref, (a, b) in evaluation_regions.items()
        },
        "tools_loaded": sorted([t for t, df in tool_outputs.items() if not df.empty]),
        "tool_score_overrides": {
            "differr": {
                "score_field": str(differr_score_field),
                "score_type": "g_stat" if str(differr_score_field) == "g_stat" else "neglog10_fdr",
                "higher_is_better": True,
            }
        },
    }
    with open(metadata_dir / "run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(run_meta, handle, indent=2)

    logger.info("Downstream analysis complete: %s", output_dir)


def _load_tool_outputs_for_mode(
    args: argparse.Namespace,
    differr_score_field: str,
    reference_aliases: Optional[Dict[str, str]] = None,
) -> Dict[str, pd.DataFrame]:
    if args.input_dir:
        in_dir = Path(args.input_dir)
        if not in_dir.exists():
            raise FileNotFoundError(f"Input directory not found: {in_dir}")
        return load_all_tool_outputs(
            in_dir,
            SUPPORTED_TOOLS,
            coverage_dirs=args.coverage_dirs,
            differr_score_field=differr_score_field,
            reference_aliases=reference_aliases,
        )

    return load_tool_outputs_from_args(
        args,
        differr_score_field=differr_score_field,
        reference_aliases=reference_aliases,
    )


def _resolve_differr_runs(
    output_dir: Path,
    differr_score_field: str,
) -> List[Tuple[str, Path]]:
    requested = str(differr_score_field).strip()
    if requested == DIFFERR_BOTH_SCORE_FIELD:
        return [
            (DIFFERR_PRIMARY_SCORE_FIELD, output_dir),
            (DIFFERR_SECONDARY_SCORE_FIELD, output_dir / DIFFERR_GSTAT_SUBDIR),
        ]
    if requested in {DIFFERR_PRIMARY_SCORE_FIELD, DIFFERR_SECONDARY_SCORE_FIELD}:
        return [(requested, output_dir)]
    raise ValueError(
        f"Unsupported DiffErr score field mode: {requested}. "
        f"Expected one of {DIFFERR_PRIMARY_SCORE_FIELD}, "
        f"{DIFFERR_SECONDARY_SCORE_FIELD}, {DIFFERR_BOTH_SCORE_FIELD}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run downstream analysis")

    input_group = parser.add_argument_group("Input")
    input_group.add_argument("--input-dir", "-i", help="modifications directory")
    for tool in SUPPORTED_TOOLS:
        input_group.add_argument(f"--{tool}", help=f"{tool} output path")

    parser.add_argument("--ground-truth", "-g", help="ground truth CSV/TSV")
    parser.add_argument("--references-csv", help="references CSV (target,reference) for canonical IDs and lengths")

    parser.add_argument("--output-dir", "-o", required=True)
    parser.add_argument("--references", "-r", nargs="+", help="optional subset of canonical references")
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--expected-replicates", type=int, default=3)
    parser.add_argument("--threshold", type=float, help="score threshold for significance")
    parser.add_argument("--coverage-dirs", action="store_true")
    parser.add_argument(
        "--differr-score-field",
        choices=[DIFFERR_PRIMARY_SCORE_FIELD, DIFFERR_SECONDARY_SCORE_FIELD, DIFFERR_BOTH_SCORE_FIELD],
        default=DIFFERR_BOTH_SCORE_FIELD,
        help=(
            "DiffErr score mode: "
            "'g_fdr_neglog10' (single run), "
            "'g_stat' (single run), or "
            "'both' (write canonical g_fdr_neglog10 and differr_gstat sub-run)"
        ),
    )

    parser.add_argument("--run-id", help="run identifier for collation tables")
    parser.add_argument("--coverage-label", help="coverage label for collation (e.g., 100x)")
    parser.add_argument("--quality-label", help="quality label for collation")

    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    ground_truth = None
    if args.ground_truth:
        gt_path = Path(args.ground_truth)
        if gt_path.exists():
            ground_truth = load_ground_truth(gt_path)
        else:
            logger.warning("Ground truth not found: %s", gt_path)

    references_csv = Path(args.references_csv) if args.references_csv else None

    # Build reference aliases up-front so filename-based reference inference
    # (e.g. Tombo) routes rows to canonical FASTA-header IDs that match the
    # reference catalog. Failure here is non-fatal — without aliases the
    # parser falls back to the legacy 16s/23s heuristic.
    reference_aliases: Optional[Dict[str, str]] = None
    if references_csv is not None:
        try:
            reference_aliases = build_reference_aliases(references_csv)
            logger.info(
                "Built %d reference alias entries from %s",
                len(reference_aliases),
                references_csv,
            )
        except Exception as exc:
            logger.warning("Failed to build reference aliases from %s: %s", references_csv, exc)

    try:
        planned_runs = _resolve_differr_runs(Path(args.output_dir), args.differr_score_field)
    except ValueError as exc:
        logger.error(str(exc))
        sys.exit(1)

    for score_field, mode_output_dir in planned_runs:
        logger.info(
            "Running downstream analysis with DiffErr score field '%s' -> %s",
            score_field,
            mode_output_dir,
        )
        try:
            tool_outputs = _load_tool_outputs_for_mode(
                args,
                differr_score_field=score_field,
                reference_aliases=reference_aliases,
            )
        except FileNotFoundError as exc:
            logger.error(str(exc))
            sys.exit(1)

        if not tool_outputs:
            logger.error("No tool outputs loaded for DiffErr score field '%s'", score_field)
            sys.exit(1)

        run_analysis(
            tool_outputs=tool_outputs,
            ground_truth=ground_truth,
            output_dir=mode_output_dir,
            references_csv=references_csv,
            references_filter=args.references,
            min_replicates=args.min_replicates,
            score_threshold=args.threshold,
            expected_replicates=args.expected_replicates,
            run_id=args.run_id,
            coverage_label=args.coverage_label,
            quality_label=args.quality_label,
            differr_score_field=score_field,
        )


if __name__ == "__main__":
    main()
