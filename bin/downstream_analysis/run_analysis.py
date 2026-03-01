#!/usr/bin/env python3
"""
Main Analysis Script for RNA Modification Detection Pipeline

Orchestrates all downstream analyses including:
- Output parsing from 9 modification detection tools
- Benchmark metrics calculation (AUPRC, AUROC, F1)
- Cross-replicate analysis and consensus calling
- Cross-tool comparison
- Visualization and report generation

Usage:
    python run_analysis.py \\
        --input-dir /path/to/modifications \\
        --ground-truth /path/to/known_sites.csv \\
        --output-dir /path/to/downstream_analysis

    Or provide individual tool outputs:
    python run_analysis.py \\
        --tombo /path/to/tombo.csv \\
        --drummer /path/to/drummer_dir \\
        --ground-truth /path/to/known_sites.csv \\
        --output-dir /path/to/downstream_analysis
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from parse_outputs import load_all_tool_outputs, parse_tool_output
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
from tool_comparison import (
    compare_tools,
    generate_upset_data,
    calculate_pairwise_agreement,
    rank_tools,
    get_consensus_positions,
)
from visualization import (
    plot_pr_curves,
    plot_roc_curves,
    plot_score_distributions,
    plot_upset,
    plot_venn,
    plot_metrics_comparison,
    plot_replicate_concordance,
    plot_correlation_heatmap,
    plot_coverage_performance,
    plot_coverage_heatmap,
    plot_coverage_stability,
    plot_saturation_curves,
    plot_multi_metric_coverage,
    generate_html_report,
)
from data_quality import (
    DataQualityReport,
    validate_all_outputs,
    generate_quality_summary,
)
from coverage_analysis import (
    analyze_coverage,
    compare_tools_by_coverage,
    generate_coverage_summary,
    get_coverage_levels,
)
from position_standardization import (
    standardize_all_tool_outputs,
    StandardizationReport,
    generate_availability_report,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
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

# Tool names
SUPPORTED_TOOLS = [
    'tombo', 'yanocomp', 'nanocompore', 'xpore', 'eligos',
    'epinano', 'differr', 'drummer', 'jacusa2'
]


def load_tool_outputs_from_args(args: argparse.Namespace) -> Dict[str, pd.DataFrame]:
    """Load tool outputs from command line arguments."""
    tool_outputs = {}

    for tool in SUPPORTED_TOOLS:
        arg_name = tool.replace('-', '_')
        filepath = getattr(args, arg_name, None)

        if filepath:
            filepath = Path(filepath)
            if filepath.exists():
                try:
                    df = parse_tool_output(tool, filepath)
                    if not df.empty:
                        tool_outputs[tool] = df
                        logger.info(f"Loaded {len(df)} positions from {tool}")
                except Exception as e:
                    logger.error(f"Failed to load {tool}: {e}")

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
    eval_end: int,
) -> Set[int]:
    if not ground_truth_positions:
        return set()

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
    eval_end: int,
) -> Set[int]:
    if not ground_truth_positions:
        return set()

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
    references: list = None,
    min_replicates: int = 2,
    score_threshold: float = None,
    expected_replicates: int = 3
) -> None:
    """
    Run complete downstream analysis.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Ground truth DataFrame (optional)
        output_dir: Output directory
        references: List of references to analyze
        min_replicates: Minimum replicates for consensus
        score_threshold: Score threshold for significance
        expected_replicates: Expected number of replicates per tool (for quality checks)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create subdirectories
    metrics_dir = output_dir / 'metrics'
    comparison_dir = output_dir / 'tool_comparison'
    replicate_dir = output_dir / 'replicate_analysis'
    viz_dir = output_dir / 'visualizations'

    for d in [metrics_dir, comparison_dir, replicate_dir, viz_dir]:
        d.mkdir(exist_ok=True)

    logger.info(f"Running analysis with {len(tool_outputs)} tools")
    logger.info(f"Output directory: {output_dir}")

    # =========================================================================
    # 0. Data Quality Validation
    # =========================================================================
    logger.info("Validating data quality...")

    quality_report = validate_all_outputs(tool_outputs, expected_replicates=expected_replicates)

    # Save quality report and warnings
    quality_report.save(output_dir / 'data_quality_report.csv')
    quality_report.save_warnings(output_dir / 'analysis_warnings.txt')

    # Generate and save quality summary
    quality_summary = generate_quality_summary(quality_report)
    with open(output_dir / 'quality_summary.txt', 'w') as f:
        f.write(quality_summary)

    # Log quality summary
    logger.info(f"Data Quality Summary:\n{quality_report.summary()}")

    # Log any critical warnings
    critical_warnings = [w for w in quality_report.warnings if w.level.value == 'ERROR']
    if critical_warnings:
        for warning in critical_warnings:
            logger.error(f"[{warning.tool}] {warning.message}")

    # Filter out tools with FAILED status (keep PARTIAL and OK)
    valid_tool_outputs = {}
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
                    for col in [
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
                        "auprc_reported",
                        "auroc_reported",
                        "f1_reported",
                    ]:
                        metric_row[col] = np.nan

                metric_row["metrics_valid"] = bool(metrics_valid)
                metric_row["invalid_score_fraction"] = invalid_fraction
                metric_row["invalid_score_n"] = invalid_n
                metric_row["invalid_score_total"] = invalid_total
                metric_row["metric_scope_note"] = metric_scope_note
                metric_rows.append(pd.DataFrame([metric_row]))

                if metrics_valid:
                    curves = _extract_curve_points(metric_ready, tool=tool, reference=ref, replicate=rep)
                    if not curves["pr"].empty:
                        pr_rows.append(curves["pr"])
                    if not curves["roc"].empty:
                        roc_rows.append(curves["roc"])
                else:
                    logger.warning(
                        "Skipping curve generation for %s/%s/%s due to invalid score fraction %.3f (%d/%d)",
                        tool,
                        ref,
                        rep,
                        invalid_fraction,
                        invalid_n,
                        invalid_total,
                    )

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
            valid_tool_outputs[tool_name] = tool_df

    if not valid_tool_outputs:
        logger.error("No valid tool outputs to analyze after quality filtering")
        return

    tool_outputs = valid_tool_outputs
    logger.info(f"Proceeding with {len(tool_outputs)} valid tools")

    # =========================================================================
    # 0.5. Position Standardization
    # =========================================================================
    logger.info("Standardizing positions across replicates and tools...")

    standardized_outputs, std_report = standardize_all_tool_outputs(
        tool_outputs, ground_truth
    )

    # Save standardization reports
    std_report.save(output_dir)

    # Log summary
    logger.info(f"Position Standardization Summary:\n{std_report.summary()}")

    # Use standardized outputs for all downstream analysis
    tool_outputs = standardized_outputs

    # =========================================================================
    # 1. Benchmark Metrics (if ground truth provided)
    # =========================================================================
    if ground_truth is not None and not ground_truth.empty:
        logger.info("Calculating benchmark metrics...")

        # Calculate metrics for all tools
        metrics_df = calculate_metrics_all_tools(tool_outputs, ground_truth, references)

        if not metrics_df.empty:
            metrics_df.to_csv(metrics_dir / 'all_metrics.csv', index=False)
            logger.info(f"Saved metrics to {metrics_dir / 'all_metrics.csv'}")

            # Generate metrics comparison plot
            plot_metrics_comparison(metrics_df, viz_dir / 'metrics_comparison.png')

            # PR and ROC curves
            for ref in (references or [None]):
                ref_suffix = f"_{ref}" if ref else ""
                plot_pr_curves(tool_outputs, ground_truth,
                              viz_dir / f'pr_curves{ref_suffix}.png', ref)
                plot_roc_curves(tool_outputs, ground_truth,
                               viz_dir / f'roc_curves{ref_suffix}.png', ref)

            # Tool ranking
            ranking = rank_tools(tool_outputs, ground_truth)
            ranking.to_csv(metrics_dir / 'tool_ranking.csv', index=False)

    # =========================================================================
    # 2. Replicate Analysis
    # =========================================================================
    logger.info("Analyzing replicates...")

    concordance_results = []
    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            continue

        # Check if there are multiple replicates
        replicates = tool_df['replicate'].unique()
        if len(replicates) > 1:
            # Calculate concordance
            stats = calculate_concordance(tool_df, score_threshold)
            concordance_results.append(stats.to_dict())

            # Consensus calling
            consensus = consensus_calling(tool_df, min_replicates, score_threshold)
            if not consensus.empty:
                consensus.to_csv(
                    replicate_dir / f'{tool_name}_consensus.csv',
                    index=False
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
                            for c in metric_cols:
                                row[c] = np.nan
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
                            for c in metric_cols:
                                row[c] = np.nan
                            row["metric_scope_note"] = "invalid_scores"
                        lag_rows_ref.append(pd.DataFrame([row]))

    # =========================================================================
    # 3. Tool Comparison
    # =========================================================================
    logger.info("Comparing tools...")

    # Generate UpSet data
    upset_df = generate_upset_data(tool_outputs, score_threshold)
    if not upset_df.empty:
        upset_df.to_csv(comparison_dir / 'upset_data.csv', index=False)
        try:
            plot_upset(tool_outputs, viz_dir / 'upset_plot.png',
                      threshold=score_threshold)
        except Exception as e:
            logger.warning(f"Could not generate UpSet plot: {e}")

    # Pairwise agreement
    pairwise = calculate_pairwise_agreement(tool_outputs, score_threshold)
    if not pairwise.empty:
        pairwise.to_csv(comparison_dir / 'pairwise_agreement.csv', index=False)
        plot_correlation_heatmap(tool_outputs, viz_dir / 'tool_correlation.png',
                                threshold=score_threshold)

    # Consensus positions across tools
    for min_tools in [2, 3, len(tool_outputs)]:
        if min_tools <= len(tool_outputs):
            consensus = get_consensus_positions(
                tool_outputs, min_tools=min_tools, threshold=score_threshold
            )
            if not consensus.empty:
                consensus.to_csv(
                    comparison_dir / f'consensus_{min_tools}tools.csv',
                    index=False
                )

    # Tool comparison summary
    comparison = compare_tools(tool_outputs, score_threshold)
    summary_data = [{
        'metric': 'n_tools',
        'value': comparison.n_tools
    }, {
        'metric': 'n_union',
        'value': comparison.n_union
    }, {
        'metric': 'n_intersection',
        'value': comparison.n_intersection
    }]
    for tool, n in comparison.n_positions_per_tool.items():
        summary_data.append({'metric': f'n_{tool}', 'value': n})
    pd.DataFrame(summary_data).to_csv(comparison_dir / 'summary.csv', index=False)

    # =========================================================================
    # 4. Score Distributions
    # =========================================================================
    logger.info("Generating score distributions...")
    plot_score_distributions(tool_outputs, viz_dir / 'score_distributions.png')

    # =========================================================================
    # 5. Venn Diagrams (for top tools)
    # =========================================================================
    if len(tool_outputs) >= 2:
        tools_to_compare = list(tool_outputs.keys())[:3]
        try:
            plot_venn(tool_outputs, viz_dir / 'venn_diagram.png',
                     tools=tools_to_compare, threshold=score_threshold)
        except Exception as e:
            logger.warning(f"Could not generate Venn diagram: {e}")

    # =========================================================================
    # 6. Coverage Analysis (if coverage data available)
    # =========================================================================
    coverage_dir = output_dir / 'coverage_analysis'
    has_coverage_data = False

    # Check if any tool has coverage data
    for tool_df in tool_outputs.values():
        if 'coverage' in tool_df.columns and tool_df['coverage'].notna().any():
            has_coverage_data = True
            break

    if has_coverage_data and ground_truth is not None:
        logger.info("Running coverage analysis...")
        coverage_dir.mkdir(exist_ok=True)

        try:
            # Analyze coverage
            coverage_result = analyze_coverage(tool_outputs, ground_truth, references)

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
    }
    with open(metadata_dir / "run_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(run_meta, handle, indent=2)

            # Save comparison table
            comparison_table = compare_tools_by_coverage(coverage_result)
            if not comparison_table.empty:
                comparison_table.to_csv(coverage_dir / 'coverage_comparison.csv')

            # Save summary
            cov_summary = generate_coverage_summary(coverage_result)
            with open(coverage_dir / 'coverage_summary.txt', 'w') as f:
                f.write(cov_summary)

            # Generate visualizations
            cov_metrics_df = coverage_result.to_dataframe()

            if not cov_metrics_df.empty:
                # Performance vs coverage
                plot_coverage_performance(
                    cov_metrics_df,
                    viz_dir / 'coverage_performance.png'
                )

                # Coverage heatmap
                plot_coverage_heatmap(
                    cov_metrics_df,
                    viz_dir / 'coverage_heatmap.png'
                )

                # Saturation curves with marked saturation points
                plot_saturation_curves(
                    cov_metrics_df,
                    viz_dir / 'saturation_curves.png',
                    saturation_points=coverage_result.saturation_points
                )

                # Stability comparison
                plot_coverage_stability(
                    coverage_result.coverage_stability,
                    viz_dir / 'coverage_stability.png'
                )

                # Multi-metric coverage plot
                plot_multi_metric_coverage(
                    cov_metrics_df,
                    viz_dir / 'multi_metric_coverage.png'
                )

            logger.info(f"Coverage analysis complete: {coverage_dir}")

        except Exception as e:
            logger.warning(f"Coverage analysis failed: {e}")
    elif has_coverage_data:
        logger.info("Coverage data found but no ground truth - skipping coverage analysis")

    # =========================================================================
    # 7. HTML Report
    # =========================================================================
    logger.info("Generating HTML report...")
    if ground_truth is not None:
        report_path = generate_html_report(
            tool_outputs, ground_truth, output_dir / 'report'
        )
        logger.info(f"Generated HTML report: {report_path}")

    # =========================================================================
    # Summary
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Analysis Complete!")
    logger.info("=" * 60)
    logger.info(f"Tools analyzed: {len(tool_outputs)}")
    logger.info(f"Output directory: {output_dir}")

    # Quality summary
    n_warnings = len(quality_report.warnings)
    n_errors = len([w for w in quality_report.warnings if w.level.value == 'ERROR'])
    if n_warnings > 0:
        logger.info(f"Data quality: {n_warnings} warnings ({n_errors} errors)")
        logger.info(f"See {output_dir / 'analysis_warnings.txt'} for details")

    if ground_truth is not None and 'metrics_df' in dir() and not metrics_df.empty:
        best_tool = metrics_df.loc[metrics_df['auprc'].idxmax()]
        logger.info(f"Best performing tool: {best_tool['tool']} (AUPRC={best_tool['auprc']:.4f})")

    logger.info(f"Total unique positions: {comparison.n_union}")
    logger.info(f"Positions detected by all tools: {comparison.n_intersection}")


def main():
    parser = argparse.ArgumentParser(
        description='Run downstream analysis for RNA modification detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input options
    input_group = parser.add_argument_group('Input Options')
    input_group.add_argument(
        '--input-dir', '-i',
        help='Directory containing all tool outputs (standard pipeline output structure)'
    )

    # Individual tool inputs
    for tool in SUPPORTED_TOOLS:
        input_group.add_argument(
            f'--{tool}',
            help=f'Path to {tool} output file or directory'
        )

    # Ground truth
    parser.add_argument(
        '--ground-truth', '-g',
        help='Path to ground truth CSV/TSV file with known modification sites'
    )

    # Output options
    parser.add_argument(
        '--output-dir', '-o',
        required=True,
        help='Output directory for analysis results'
    )

    # Analysis options
    parser.add_argument(
        '--references', '-r',
        nargs='+',
        help='Filter to specific references (e.g., 16s 23s)'
    )
    parser.add_argument(
        '--min-replicates',
        type=int,
        default=2,
        help='Minimum replicates for consensus calling (default: 2)'
    )
    parser.add_argument(
        '--expected-replicates',
        type=int,
        default=3,
        help='Expected number of replicates per tool for quality checks (default: 3)'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        help='Score threshold for significance (p-value or FDR cutoff)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '--coverage-dirs',
        action='store_true',
        help='Expect {coverage}/tool/ directory structure for coverage analysis'
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Load tool outputs
    tool_outputs = {}

    if args.input_dir:
        # Load from standard directory structure
        input_dir = Path(args.input_dir)
        if not input_dir.exists():
            logger.error(f"Input directory not found: {input_dir}")
            sys.exit(1)

        tool_outputs = load_all_tool_outputs(
            input_dir, SUPPORTED_TOOLS, coverage_dirs=args.coverage_dirs
        )
    else:
        # Load from individual arguments
        tool_outputs = load_tool_outputs_from_args(args)

    if not tool_outputs:
        logger.error("No tool outputs loaded. Provide --input-dir or individual tool paths.")
        sys.exit(1)

    # Load ground truth
    ground_truth = None
    if args.ground_truth:
        ground_truth_path = Path(args.ground_truth)
        if ground_truth_path.exists():
            ground_truth = load_ground_truth(ground_truth_path)
        else:
            logger.warning(f"Ground truth file not found: {ground_truth_path}")

    # Run analysis
    run_analysis(
        tool_outputs=tool_outputs,
        ground_truth=ground_truth,
        output_dir=Path(args.output_dir),
        references=args.references,
        min_replicates=args.min_replicates,
        score_threshold=args.threshold,
        expected_replicates=args.expected_replicates
    )


if __name__ == '__main__':
    main()
