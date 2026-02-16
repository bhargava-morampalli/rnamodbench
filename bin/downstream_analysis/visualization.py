#!/usr/bin/env python3
"""
Visualization Module for RNA Modification Detection

Generates publication-quality plots for:
- Precision-Recall curves
- ROC curves
- Score distributions
- UpSet plots for tool comparison
- Venn diagrams
- Coverage vs performance plots
- Correlation heatmaps
- Summary reports

Usage:
    from visualization import plot_pr_curves, plot_upset, generate_report

    plot_pr_curves(tool_outputs, ground_truth, output_path)
    plot_upset(tool_outputs, output_path)
"""

import logging
import math
from pathlib import Path
from typing import Union, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Set default style
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    plt.style.use('seaborn-v0_8-whitegrid')
except ImportError:
    logger.warning("matplotlib not available - plotting functions will not work")

# Color palette for tools
TOOL_COLORS = {
    'tombo': '#1f77b4',      # Blue
    'yanocomp': '#ff7f0e',   # Orange
    'nanocompore': '#2ca02c', # Green
    'xpore': '#d62728',      # Red
    'eligos': '#9467bd',     # Purple
    'epinano': '#8c564b',    # Brown
    'differr': '#e377c2',    # Pink
    'drummer': '#7f7f7f',    # Gray
    'jacusa2': '#bcbd22',    # Yellow-green
}


def get_tool_color(tool_name: str) -> str:
    """Get color for a tool."""
    return TOOL_COLORS.get(tool_name.lower(), '#17becf')


REPLICATE_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
REPLICATE_LINESTYLES = ['-', '--', '-.', ':']
DISPLAY_NAMES = {
    "tombo": "Tombo",
    "yanocomp": "Yanocomp",
    "nanocompore": "Nanocompore",
    "xpore": "xPore",
    "eligos": "ELIGOS",
    "epinano": "EpiNano",
    "differr": "DiffErr",
    "drummer": "DRUMMER",
    "jacusa2": "JACUSA2",
}


def _apply_nature_minimal_style() -> bool:
    """Apply deterministic publication-style rcParams."""
    if "plt" not in globals() or "matplotlib" not in globals():
        logger.warning("matplotlib not available - cannot render grid figures")
        return False

    matplotlib.rcParams.update(
        {
            "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "axes.titlesize": 9,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "figure.titlesize": 11,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "grid.color": "#9a9a9a",
            "grid.alpha": 0.2,
            "grid.linewidth": 0.6,
            "axes.grid": False,
        }
    )
    return True


def _fmt_metric(value: Union[float, int, None]) -> str:
    """Format metric with NA handling."""
    try:
        if value is None or pd.isna(value):
            return "NA"
        return f"{float(value):.3f}"
    except Exception:
        return "NA"


def _replicate_sort_key(rep: str) -> Tuple[int, str]:
    rep_s = str(rep).strip().lower()
    if rep_s.startswith("rep"):
        suffix = rep_s[3:]
        if suffix.isdigit():
            return (int(suffix), rep_s)
    digits = "".join(ch for ch in rep_s if ch.isdigit())
    if digits:
        return (int(digits), rep_s)
    return (999, rep_s)


def _replicate_style(rep: str, idx: int) -> Tuple[str, str]:
    # fixed colors for rep1/rep2/rep3, fall back by index
    rep_s = str(rep).strip().lower()
    rep_to_idx = {"rep1": 0, "rep2": 1, "rep3": 2}
    use_idx = rep_to_idx.get(rep_s, idx)
    color = REPLICATE_COLORS[use_idx % len(REPLICATE_COLORS)]
    linestyle = REPLICATE_LINESTYLES[use_idx % len(REPLICATE_LINESTYLES)]
    return color, linestyle


def _display_name(tool: str) -> str:
    return DISPLAY_NAMES.get(str(tool).lower(), str(tool))


def _metric_summary(tool_metrics: pd.DataFrame, metric_col: str) -> Tuple[str, str]:
    vals = pd.to_numeric(tool_metrics.get(metric_col, np.nan), errors="coerce").dropna()
    if vals.empty:
        return ("NA", "NA")
    mean_v = float(vals.mean())
    # pandas std defaults to ddof=1, keep this for replicate variation
    sd_v = float(vals.std()) if len(vals) > 1 else 0.0
    return (_fmt_metric(mean_v), _fmt_metric(sd_v))


def _figure_context_text(metrics_df: pd.DataFrame, reference: str) -> str:
    ref_df = metrics_df.copy()
    if "reference" in ref_df.columns:
        ref_df = ref_df[ref_df["reference"].astype(str) == str(reference)]

    gt_total = pd.to_numeric(ref_df.get("n_gt_total", np.nan), errors="coerce").dropna()
    n_universe = pd.to_numeric(ref_df.get("n_universe", np.nan), errors="coerce").dropna()

    gt_txt = str(int(gt_total.iloc[0])) if not gt_total.empty else "NA"
    universe_txt = str(int(n_universe.iloc[0])) if not n_universe.empty else "NA"
    return f"Reference: {reference} | Ground-truth positives: {gt_txt} | Universe positions: {universe_txt}"


def _save_figure_outputs(fig, output_prefix: Path) -> List[Path]:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for ext in ("pdf", "svg", "png"):
        out = output_prefix.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.05)
        else:
            fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
        written.append(out)
    return written


def plot_roc_tool_grid(
    curve_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    reference: str,
    output_dir: Union[str, Path],
    ncols: int = 3,
    output_stem: str = "roc_tool_grid",
) -> List[Path]:
    """
    Plot per-tool ROC subplots with replicate overlays for a single reference.
    Returns list of generated file paths.
    """
    if not _apply_nature_minimal_style():
        return []

    if metrics_df is None or metrics_df.empty:
        logger.warning("No metrics data for ROC tool grid (%s)", reference)
        return []

    output_dir = Path(output_dir)
    ref_metrics = metrics_df.copy()
    if "reference" in ref_metrics.columns:
        ref_metrics = ref_metrics[ref_metrics["reference"].astype(str) == str(reference)]
    if ref_metrics.empty:
        logger.warning("No reference-specific metrics for ROC grid (%s)", reference)
        return []

    if curve_df is None:
        curve_df = pd.DataFrame()
    ref_curves = curve_df.copy()
    if not ref_curves.empty and "reference" in ref_curves.columns:
        ref_curves = ref_curves[ref_curves["reference"].astype(str) == str(reference)]

    tools = sorted(ref_metrics["tool"].dropna().astype(str).unique().tolist(), key=str.lower)
    if not tools:
        logger.warning("No tools available for ROC grid (%s)", reference)
        return []

    n_tools = len(tools)
    nrows = int(math.ceil(n_tools / max(1, ncols)))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 3.5 * nrows),
        constrained_layout=False,
        squeeze=False,
    )
    ax_list = axes.flatten()

    for idx, tool in enumerate(tools):
        ax = ax_list[idx]
        tool_metrics = ref_metrics[ref_metrics["tool"].astype(str) == tool].copy()
        tool_curves = (
            ref_curves[ref_curves["tool"].astype(str) == tool].copy() if not ref_curves.empty else pd.DataFrame()
        )

        reps = sorted(
            tool_metrics["replicate"].astype(str).unique().tolist() if "replicate" in tool_metrics.columns else [],
            key=_replicate_sort_key,
        )
        plotted = 0

        for rep_idx, rep in enumerate(reps):
            rep_curve = (
                tool_curves[tool_curves["replicate"].astype(str) == str(rep)].copy()
                if not tool_curves.empty and "replicate" in tool_curves.columns
                else pd.DataFrame()
            )
            if rep_curve.empty:
                continue

            sort_cols = [c for c in ("index", "fpr", "tpr") if c in rep_curve.columns]
            if sort_cols:
                rep_curve = rep_curve.sort_values(sort_cols)

            fpr = pd.to_numeric(rep_curve.get("fpr", np.nan), errors="coerce")
            tpr = pd.to_numeric(rep_curve.get("tpr", np.nan), errors="coerce")
            valid = (~fpr.isna()) & (~tpr.isna())
            if valid.sum() < 2:
                continue

            rep_metrics = tool_metrics[tool_metrics["replicate"].astype(str) == str(rep)]
            rep_auroc = _fmt_metric(
                pd.to_numeric(rep_metrics.get("auroc", np.nan), errors="coerce").dropna().iloc[0]
                if not rep_metrics.empty
                else np.nan
            )

            color, linestyle = _replicate_style(rep, rep_idx)
            ax.plot(
                fpr[valid].to_numpy(),
                tpr[valid].to_numpy(),
                color=color,
                linestyle=linestyle,
                linewidth=1.3,
                label=f"{rep} (AUROC={rep_auroc})",
            )
            plotted += 1

        ax.plot([0, 1], [0, 1], color="#666666", linestyle="--", linewidth=0.9, alpha=0.8)

        if plotted == 0:
            ax.text(0.5, 0.5, "No curve data", ha="center", va="center", fontsize=8, color="#444444")

        auroc_mean, auroc_sd = _metric_summary(tool_metrics, "auroc")
        ax.set_title(f"{_display_name(tool)} | mean AUROC={auroc_mean} \u00b1 {auroc_sd}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ticks = np.linspace(0, 1, 6)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.grid(True, which="major", alpha=0.2)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            inside = len(labels) <= 3
            if inside:
                ax.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.9)
            else:
                ax.legend(
                    loc="center left",
                    bbox_to_anchor=(1.02, 0.5),
                    borderaxespad=0.0,
                    frameon=True,
                    facecolor="white",
                    framealpha=0.9,
                )

    for idx in range(n_tools, len(ax_list)):
        ax_list[idx].set_visible(False)

    fig.suptitle("ROC Curves by Tool (Replicate Overlay)", fontsize=11, y=0.988)
    fig.text(0.5, 0.963, _figure_context_text(ref_metrics, reference), ha="center", va="top", fontsize=8)
    fig.text(
        0.5,
        0.944,
        "Interpretation: higher AUROC indicates better rank-order separation.",
        ha="center",
        va="top",
        fontsize=8,
    )
    fig.tight_layout(rect=[0.01, 0.02, 0.99, 0.92])

    written = _save_figure_outputs(fig, Path(output_dir) / output_stem)
    plt.close(fig)
    logger.info("Saved ROC tool grid for %s: %s", reference, ", ".join(str(p) for p in written))
    return written


def plot_pr_tool_grid(
    curve_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    reference: str,
    output_dir: Union[str, Path],
    ncols: int = 3,
    output_stem: str = "pr_tool_grid",
) -> List[Path]:
    """
    Plot per-tool PR subplots with replicate overlays for a single reference.
    Returns list of generated file paths.
    """
    if not _apply_nature_minimal_style():
        return []

    if metrics_df is None or metrics_df.empty:
        logger.warning("No metrics data for PR tool grid (%s)", reference)
        return []

    output_dir = Path(output_dir)
    ref_metrics = metrics_df.copy()
    if "reference" in ref_metrics.columns:
        ref_metrics = ref_metrics[ref_metrics["reference"].astype(str) == str(reference)]
    if ref_metrics.empty:
        logger.warning("No reference-specific metrics for PR grid (%s)", reference)
        return []

    if curve_df is None:
        curve_df = pd.DataFrame()
    ref_curves = curve_df.copy()
    if not ref_curves.empty and "reference" in ref_curves.columns:
        ref_curves = ref_curves[ref_curves["reference"].astype(str) == str(reference)]

    tools = sorted(ref_metrics["tool"].dropna().astype(str).unique().tolist(), key=str.lower)
    if not tools:
        logger.warning("No tools available for PR grid (%s)", reference)
        return []

    n_tools = len(tools)
    nrows = int(math.ceil(n_tools / max(1, ncols)))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 3.5 * nrows),
        constrained_layout=False,
        squeeze=False,
    )
    ax_list = axes.flatten()

    gt_total_series = pd.to_numeric(ref_metrics.get("n_gt_total", np.nan), errors="coerce").dropna()
    n_universe_series = pd.to_numeric(ref_metrics.get("n_universe", np.nan), errors="coerce").dropna()
    prevalence = np.nan
    if not gt_total_series.empty and not n_universe_series.empty:
        denom = float(n_universe_series.iloc[0])
        if denom > 0:
            prevalence = float(gt_total_series.iloc[0] / denom)

    for idx, tool in enumerate(tools):
        ax = ax_list[idx]
        tool_metrics = ref_metrics[ref_metrics["tool"].astype(str) == tool].copy()
        tool_curves = (
            ref_curves[ref_curves["tool"].astype(str) == tool].copy() if not ref_curves.empty else pd.DataFrame()
        )

        reps = sorted(
            tool_metrics["replicate"].astype(str).unique().tolist() if "replicate" in tool_metrics.columns else [],
            key=_replicate_sort_key,
        )
        plotted = 0

        for rep_idx, rep in enumerate(reps):
            rep_curve = (
                tool_curves[tool_curves["replicate"].astype(str) == str(rep)].copy()
                if not tool_curves.empty and "replicate" in tool_curves.columns
                else pd.DataFrame()
            )
            if rep_curve.empty:
                continue

            sort_cols = [c for c in ("index", "recall", "precision") if c in rep_curve.columns]
            if sort_cols:
                rep_curve = rep_curve.sort_values(sort_cols)

            recall = pd.to_numeric(rep_curve.get("recall", np.nan), errors="coerce")
            precision = pd.to_numeric(rep_curve.get("precision", np.nan), errors="coerce")
            valid = (~recall.isna()) & (~precision.isna())
            if valid.sum() < 2:
                continue

            rep_metrics = tool_metrics[tool_metrics["replicate"].astype(str) == str(rep)]
            rep_auprc = _fmt_metric(
                pd.to_numeric(rep_metrics.get("auprc", np.nan), errors="coerce").dropna().iloc[0]
                if not rep_metrics.empty
                else np.nan
            )

            color, linestyle = _replicate_style(rep, rep_idx)
            ax.plot(
                recall[valid].to_numpy(),
                precision[valid].to_numpy(),
                color=color,
                linestyle=linestyle,
                linewidth=1.3,
                label=f"{rep} (AUPRC={rep_auprc})",
            )
            plotted += 1

        if not pd.isna(prevalence):
            ax.axhline(prevalence, color="#666666", linestyle="--", linewidth=0.9, alpha=0.8)

        if plotted == 0:
            ax.text(0.5, 0.5, "No curve data", ha="center", va="center", fontsize=8, color="#444444")

        auprc_mean, auprc_sd = _metric_summary(tool_metrics, "auprc")
        ax.set_title(f"{_display_name(tool)} | mean AUPRC={auprc_mean} \u00b1 {auprc_sd}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ticks = np.linspace(0, 1, 6)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.grid(True, which="major", alpha=0.2)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            inside = len(labels) <= 3
            if inside:
                ax.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.9)
            else:
                ax.legend(
                    loc="center left",
                    bbox_to_anchor=(1.02, 0.5),
                    borderaxespad=0.0,
                    frameon=True,
                    facecolor="white",
                    framealpha=0.9,
                )

    for idx in range(n_tools, len(ax_list)):
        ax_list[idx].set_visible(False)

    fig.suptitle("Precision-Recall Curves by Tool (Replicate Overlay)", fontsize=11, y=0.988)
    fig.text(0.5, 0.963, _figure_context_text(ref_metrics, reference), ha="center", va="top", fontsize=8)
    if not pd.isna(prevalence):
        pr_msg = f"Interpretation: higher AUPRC is better; baseline (dashed) equals prevalence={prevalence:.3f}."
    else:
        pr_msg = "Interpretation: higher AUPRC is better; prevalence baseline unavailable."
    fig.text(0.5, 0.944, pr_msg, ha="center", va="top", fontsize=8)
    fig.tight_layout(rect=[0.01, 0.02, 0.99, 0.92])

    written = _save_figure_outputs(fig, Path(output_dir) / output_stem)
    plt.close(fig)
    logger.info("Saved PR tool grid for %s: %s", reference, ", ".join(str(p) for p in written))
    return written


def plot_pr_curves(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    output_path: Union[str, Path],
    reference: str = None,
    figsize: Tuple[int, int] = (10, 8),
    title: str = None
) -> None:
    """
    Plot Precision-Recall curves for all tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Ground truth DataFrame
        output_path: Path to save the plot
        reference: Filter to specific reference
        figsize: Figure size
        title: Plot title
    """
    try:
        from .benchmark_metrics import get_precision_recall_curve, calculate_metrics
    except Exception:
        from benchmark_metrics import get_precision_recall_curve, calculate_metrics

    fig, ax = plt.subplots(figsize=figsize)

    metrics_data = []

    for tool_name, tool_df in sorted(tool_outputs.items()):
        if tool_df.empty:
            continue

        try:
            precision, recall, _ = get_precision_recall_curve(
                tool_df, ground_truth, reference
            )
            metrics = calculate_metrics(tool_df, ground_truth, reference)

            color = get_tool_color(tool_name)
            ax.plot(recall, precision, color=color, lw=2,
                    label=f'{tool_name} (AUPRC={metrics.auprc:.3f})')

            metrics_data.append({
                'tool': tool_name,
                'auprc': metrics.auprc
            })
        except Exception as e:
            logger.warning(f"Could not plot PR curve for {tool_name}: {e}")

    # Add baseline
    if metrics_data:
        n_positive = len(ground_truth)
        baseline = n_positive / (n_positive + 1000)  # Approximate
        ax.axhline(y=baseline, color='gray', linestyle='--', lw=1, alpha=0.5)

    ax.set_xlabel('Recall', fontsize=12)
    ax.set_ylabel('Precision', fontsize=12)
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.legend(loc='lower left', fontsize=10)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ref_str = f" ({reference})" if reference else ""
        ax.set_title(f'Precision-Recall Curves{ref_str}', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved PR curves to {output_path}")


def plot_roc_curves(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    output_path: Union[str, Path],
    reference: str = None,
    figsize: Tuple[int, int] = (10, 8),
    title: str = None
) -> None:
    """
    Plot ROC curves for all tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Ground truth DataFrame
        output_path: Path to save the plot
        reference: Filter to specific reference
        figsize: Figure size
        title: Plot title
    """
    try:
        from .benchmark_metrics import get_roc_curve, calculate_metrics
    except Exception:
        from benchmark_metrics import get_roc_curve, calculate_metrics

    fig, ax = plt.subplots(figsize=figsize)

    for tool_name, tool_df in sorted(tool_outputs.items()):
        if tool_df.empty:
            continue

        try:
            fpr, tpr, _ = get_roc_curve(tool_df, ground_truth, reference)
            metrics = calculate_metrics(tool_df, ground_truth, reference)

            color = get_tool_color(tool_name)
            ax.plot(fpr, tpr, color=color, lw=2,
                    label=f'{tool_name} (AUROC={metrics.auroc:.3f})')
        except Exception as e:
            logger.warning(f"Could not plot ROC curve for {tool_name}: {e}")

    # Add diagonal baseline
    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)

    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.legend(loc='lower right', fontsize=10)

    if title:
        ax.set_title(title, fontsize=14)
    else:
        ref_str = f" ({reference})" if reference else ""
        ax.set_title(f'ROC Curves{ref_str}', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved ROC curves to {output_path}")


def plot_score_distributions(
    tool_outputs: Dict[str, pd.DataFrame],
    output_path: Union[str, Path],
    reference: str = None,
    figsize: Tuple[int, int] = (14, 10),
    log_scale: bool = True
) -> None:
    """
    Plot score distributions for all tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        output_path: Path to save the plot
        reference: Filter to specific reference
        figsize: Figure size
        log_scale: Use log scale for x-axis (p-values)
    """
    import seaborn as sns

    n_tools = len([t for t, df in tool_outputs.items() if not df.empty])
    if n_tools == 0:
        logger.warning("No tool outputs to plot")
        return

    # Calculate grid size
    n_cols = min(3, n_tools)
    n_rows = (n_tools + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_tools == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    plot_idx = 0
    for tool_name, tool_df in sorted(tool_outputs.items()):
        if tool_df.empty:
            continue

        ax = axes[plot_idx]
        df = tool_df.copy()

        if reference:
            df = df[df['reference'].str.contains(reference, case=False, na=False)]

        if df.empty:
            plot_idx += 1
            continue

        scores = df['score'].dropna()
        score_type = df['score_type'].iloc[0] if 'score_type' in df.columns else 'score'

        color = get_tool_color(tool_name)

        # Transform scores based on type
        if score_type in ['pvalue', 'fdr'] and log_scale:
            # Use -log10 scale
            scores = -np.log10(scores.clip(lower=1e-300))
            xlabel = f'-log10({score_type})'
        else:
            xlabel = score_type

        sns.histplot(scores, ax=ax, color=color, kde=True, alpha=0.7)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(f'{tool_name}', fontsize=12)

        plot_idx += 1

    # Hide unused axes
    for idx in range(plot_idx, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('Score Distributions by Tool', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved score distributions to {output_path}")


def plot_upset(
    tool_outputs: Dict[str, pd.DataFrame],
    output_path: Union[str, Path],
    threshold: float = None,
    reference: str = None,
    figsize: Tuple[int, int] = (12, 8),
    min_subset_size: int = 1
) -> None:
    """
    Generate UpSet plot showing tool overlaps.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        output_path: Path to save the plot
        threshold: Score threshold
        reference: Filter to specific reference
        figsize: Figure size
        min_subset_size: Minimum subset size to show
    """
    try:
        from upsetplot import UpSet, from_memberships
    except ImportError:
        logger.error("upsetplot not installed. Run: pip install upsetplot")
        return

    try:
        from .tool_comparison import generate_upset_data
    except Exception:
        from tool_comparison import generate_upset_data

    # Generate upset data
    upset_df = generate_upset_data(tool_outputs, threshold, reference)

    if upset_df.empty:
        logger.warning("No data for UpSet plot")
        return

    # Get tool columns
    tool_cols = [c for c in upset_df.columns if c != 'position']

    # Create membership lists
    memberships = []
    for _, row in upset_df.iterrows():
        members = [tool for tool in tool_cols if row[tool]]
        if members:
            memberships.append(members)

    if not memberships:
        logger.warning("No memberships for UpSet plot")
        return

    # Create UpSet data
    upset_data = from_memberships(memberships)

    # Create plot
    fig = plt.figure(figsize=figsize)
    upset = UpSet(upset_data, min_subset_size=min_subset_size,
                  show_counts=True, sort_by='cardinality')
    upset.plot(fig=fig)

    plt.suptitle('Tool Agreement (UpSet Plot)', fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved UpSet plot to {output_path}")


def plot_venn(
    tool_outputs: Dict[str, pd.DataFrame],
    output_path: Union[str, Path],
    tools: List[str] = None,
    threshold: float = None,
    reference: str = None,
    figsize: Tuple[int, int] = (10, 8)
) -> None:
    """
    Generate Venn diagram for 2-3 tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        output_path: Path to save the plot
        tools: List of 2-3 tools to compare
        threshold: Score threshold
        reference: Filter to specific reference
        figsize: Figure size
    """
    try:
        from matplotlib_venn import venn2, venn3
    except ImportError:
        logger.error("matplotlib-venn not installed. Run: pip install matplotlib-venn")
        return

    try:
        from .tool_comparison import generate_venn_data
    except Exception:
        from tool_comparison import generate_venn_data

    # Generate Venn data
    venn_data = generate_venn_data(tool_outputs, tools, threshold, reference)

    if not venn_data:
        logger.warning("No data for Venn diagram")
        return

    labels = venn_data.get('labels', [])
    n_tools = len(labels)

    fig, ax = plt.subplots(figsize=figsize)

    if n_tools == 2:
        venn2(subsets=(
            venn_data['10'],
            venn_data['01'],
            venn_data['11']
        ), set_labels=labels, ax=ax)
    elif n_tools == 3:
        venn3(subsets=(
            venn_data['100'],
            venn_data['010'],
            venn_data['110'],
            venn_data['001'],
            venn_data['101'],
            venn_data['011'],
            venn_data['111']
        ), set_labels=labels, ax=ax)
    else:
        logger.error(f"Venn diagram supports 2-3 tools, got {n_tools}")
        plt.close()
        return

    plt.title('Tool Overlap (Venn Diagram)', fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved Venn diagram to {output_path}")


def plot_metrics_comparison(
    metrics_df: pd.DataFrame,
    output_path: Union[str, Path],
    metric: str = 'auprc',
    figsize: Tuple[int, int] = (12, 6)
) -> None:
    """
    Plot bar chart comparing metrics across tools.

    Args:
        metrics_df: DataFrame with metrics per tool
        output_path: Path to save the plot
        metric: Metric to plot
        figsize: Figure size
    """
    if metrics_df.empty:
        logger.warning("No metrics to plot")
        return

    fig, ax = plt.subplots(figsize=figsize)

    # Sort by metric
    df = metrics_df.sort_values(metric, ascending=False)

    # Get colors
    colors = [get_tool_color(t) for t in df['tool']]

    bars = ax.bar(df['tool'], df[metric], color=colors, alpha=0.8, edgecolor='black')

    # Add value labels
    for bar, val in zip(bars, df[metric]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Tool', fontsize=12)
    ax.set_ylabel(metric.upper(), fontsize=12)
    ax.set_title(f'{metric.upper()} Comparison Across Tools', fontsize=14)
    ax.set_ylim(0, 1.1 * df[metric].max())

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved metrics comparison to {output_path}")


def plot_replicate_concordance(
    concordance_df: pd.DataFrame,
    output_path: Union[str, Path],
    figsize: Tuple[int, int] = (12, 6)
) -> None:
    """
    Plot replicate concordance (Jaccard similarity) across tools.

    Args:
        concordance_df: DataFrame with concordance metrics
        output_path: Path to save the plot
        figsize: Figure size
    """
    if concordance_df.empty:
        logger.warning("No concordance data to plot")
        return

    fig, ax = plt.subplots(figsize=figsize)

    df = concordance_df.sort_values('mean_jaccard', ascending=False)
    colors = [get_tool_color(t) for t in df['tool']]

    x = range(len(df))
    bars = ax.bar(x, df['mean_jaccard'], color=colors, alpha=0.8, edgecolor='black')

    # Add error bars if std available
    if 'std_jaccard' in df.columns:
        ax.errorbar(x, df['mean_jaccard'], yerr=df['std_jaccard'],
                    fmt='none', color='black', capsize=3)

    ax.set_xticks(x)
    ax.set_xticklabels(df['tool'], rotation=45, ha='right')
    ax.set_xlabel('Tool', fontsize=12)
    ax.set_ylabel('Mean Jaccard Similarity', fontsize=12)
    ax.set_title('Replicate Concordance by Tool', fontsize=14)
    ax.set_ylim(0, 1.05)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved concordance plot to {output_path}")


def plot_correlation_heatmap(
    tool_outputs: Dict[str, pd.DataFrame],
    output_path: Union[str, Path],
    threshold: float = None,
    reference: str = None,
    figsize: Tuple[int, int] = (10, 8)
) -> None:
    """
    Plot heatmap of pairwise Jaccard similarities between tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        output_path: Path to save the plot
        threshold: Score threshold
        reference: Filter to specific reference
        figsize: Figure size
    """
    import seaborn as sns
    try:
        from .tool_comparison import calculate_pairwise_agreement
    except Exception:
        from tool_comparison import calculate_pairwise_agreement

    pairwise = calculate_pairwise_agreement(tool_outputs, threshold, reference)

    if pairwise.empty:
        logger.warning("No pairwise data for heatmap")
        return

    # Create correlation matrix
    tools = sorted(set(pairwise['tool_1'].tolist() + pairwise['tool_2'].tolist()))

    corr_matrix = pd.DataFrame(1.0, index=tools, columns=tools)

    for _, row in pairwise.iterrows():
        corr_matrix.loc[row['tool_1'], row['tool_2']] = row['jaccard']
        corr_matrix.loc[row['tool_2'], row['tool_1']] = row['jaccard']

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(corr_matrix, annot=True, fmt='.3f', cmap='RdYlGn',
                vmin=0, vmax=1, ax=ax, square=True,
                cbar_kws={'label': 'Jaccard Similarity'})

    ax.set_title('Tool Agreement (Jaccard Similarity)', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved correlation heatmap to {output_path}")


def plot_coverage_performance(
    coverage_metrics: pd.DataFrame,
    output_path: Union[str, Path],
    metric: str = 'auprc',
    figsize: Tuple[int, int] = (12, 8)
) -> None:
    """
    Plot performance vs coverage for each tool.

    Args:
        coverage_metrics: DataFrame with metrics at different coverage levels
        output_path: Path to save the plot
        metric: Metric to plot
        figsize: Figure size
    """
    if coverage_metrics.empty:
        logger.warning("No coverage metrics to plot")
        return

    fig, ax = plt.subplots(figsize=figsize)

    for tool in coverage_metrics['tool'].unique():
        tool_data = coverage_metrics[coverage_metrics['tool'] == tool]
        tool_data = tool_data.sort_values('coverage')

        color = get_tool_color(tool)
        ax.plot(tool_data['coverage'], tool_data[metric],
                marker='o', color=color, lw=2, label=tool)

    ax.set_xlabel('Coverage', fontsize=12)
    ax.set_ylabel(metric.upper(), fontsize=12)
    ax.set_title(f'{metric.upper()} vs Coverage', fontsize=14)
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xscale('log')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved coverage performance plot to {output_path}")


def plot_coverage_heatmap(
    coverage_metrics: pd.DataFrame,
    output_path: Union[str, Path],
    metric: str = 'auprc',
    figsize: Tuple[int, int] = (14, 8)
) -> None:
    """
    Plot heatmap of metric values across tools and coverage levels.

    Args:
        coverage_metrics: DataFrame with 'tool', 'coverage', and metric columns
        output_path: Path to save the plot
        metric: Metric to plot (auprc, f1, auroc, etc.)
        figsize: Figure size
    """
    import seaborn as sns

    if coverage_metrics.empty:
        logger.warning("No coverage metrics for heatmap")
        return

    # Create pivot table
    pivot = coverage_metrics.pivot_table(
        index='tool',
        columns='coverage',
        values=metric,
        aggfunc='mean'
    )

    if pivot.empty:
        logger.warning("Could not create pivot table for coverage heatmap")
        return

    # Sort columns by coverage value (numeric)
    def parse_cov(c):
        try:
            return int(str(c).lower().rstrip('x'))
        except:
            return 0

    sorted_cols = sorted(pivot.columns, key=parse_cov)
    pivot = pivot[sorted_cols]

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn',
                vmin=0, vmax=1, ax=ax,
                cbar_kws={'label': metric.upper()})

    ax.set_title(f'{metric.upper()} by Tool and Coverage Level', fontsize=14)
    ax.set_xlabel('Coverage', fontsize=12)
    ax.set_ylabel('Tool', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved coverage heatmap to {output_path}")


def plot_coverage_stability(
    stability_scores: Dict[str, float],
    output_path: Union[str, Path],
    figsize: Tuple[int, int] = (12, 6)
) -> None:
    """
    Plot bar chart of coverage stability scores for each tool.

    Args:
        stability_scores: Dictionary mapping tool names to stability scores
        output_path: Path to save the plot
        figsize: Figure size
    """
    if not stability_scores:
        logger.warning("No stability scores to plot")
        return

    # Filter out NaN values
    valid_scores = {k: v for k, v in stability_scores.items() if not pd.isna(v)}

    if not valid_scores:
        logger.warning("No valid stability scores to plot")
        return

    # Sort by stability
    sorted_items = sorted(valid_scores.items(), key=lambda x: x[1], reverse=True)
    tools = [item[0] for item in sorted_items]
    scores = [item[1] for item in sorted_items]
    colors = [get_tool_color(t) for t in tools]

    fig, ax = plt.subplots(figsize=figsize)

    bars = ax.bar(tools, scores, color=colors, alpha=0.8, edgecolor='black')

    # Add value labels
    for bar, score in zip(bars, scores):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{score:.3f}', ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Tool', fontsize=12)
    ax.set_ylabel('Stability Score', fontsize=12)
    ax.set_title('Coverage Stability by Tool\n(Higher = More Consistent Across Coverages)', fontsize=14)
    ax.set_ylim(0, 1.1)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved stability plot to {output_path}")


def plot_saturation_curves(
    coverage_metrics: pd.DataFrame,
    output_path: Union[str, Path],
    metric: str = 'auprc',
    saturation_points: Dict[str, Optional[str]] = None,
    figsize: Tuple[int, int] = (12, 8)
) -> None:
    """
    Plot performance saturation curves with saturation points marked.

    Args:
        coverage_metrics: DataFrame with coverage metrics
        output_path: Path to save the plot
        metric: Metric to plot
        saturation_points: Dictionary of tool -> saturation coverage
        figsize: Figure size
    """
    if coverage_metrics.empty:
        logger.warning("No coverage metrics for saturation plot")
        return

    fig, ax = plt.subplots(figsize=figsize)

    for tool in coverage_metrics['tool'].unique():
        tool_data = coverage_metrics[coverage_metrics['tool'] == tool].copy()

        if tool_data.empty:
            continue

        # Sort by numeric coverage
        def parse_cov(c):
            try:
                return int(str(c).lower().rstrip('x'))
            except:
                return 0

        tool_data['cov_num'] = tool_data['coverage'].apply(parse_cov)
        tool_data = tool_data.sort_values('cov_num')

        color = get_tool_color(tool)

        # Plot curve
        ax.plot(tool_data['cov_num'], tool_data[metric],
                marker='o', color=color, lw=2, label=tool, markersize=8)

        # Mark saturation point if provided
        if saturation_points and tool in saturation_points:
            sat_cov = saturation_points[tool]
            if sat_cov:
                sat_num = parse_cov(sat_cov)
                sat_val = tool_data[tool_data['cov_num'] == sat_num][metric].values
                if len(sat_val) > 0:
                    ax.axvline(x=sat_num, color=color, linestyle='--', alpha=0.5)
                    ax.scatter([sat_num], sat_val, marker='*', s=200, color=color,
                              edgecolors='black', zorder=5)

    ax.set_xlabel('Coverage (x)', fontsize=12)
    ax.set_ylabel(metric.upper(), fontsize=12)
    ax.set_title(f'{metric.upper()} vs Coverage\n(Stars indicate saturation points)', fontsize=14)
    ax.legend(loc='lower right', fontsize=10)
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved saturation curves to {output_path}")


def plot_multi_metric_coverage(
    coverage_metrics: pd.DataFrame,
    output_path: Union[str, Path],
    tool: str = None,
    metrics: List[str] = None,
    figsize: Tuple[int, int] = (14, 10)
) -> None:
    """
    Plot multiple metrics vs coverage for a specific tool or all tools.

    Args:
        coverage_metrics: DataFrame with coverage metrics
        output_path: Path to save the plot
        tool: Specific tool to plot (None for all tools averaged)
        metrics: List of metrics to plot (default: auprc, f1, precision, recall)
        figsize: Figure size
    """
    if coverage_metrics.empty:
        logger.warning("No coverage metrics for multi-metric plot")
        return

    if metrics is None:
        metrics = ['auprc', 'f1', 'precision', 'recall']

    # Filter to tool if specified
    if tool:
        df = coverage_metrics[coverage_metrics['tool'] == tool].copy()
        title_suffix = f"for {tool}"
    else:
        # Average across tools
        df = coverage_metrics.copy()
        title_suffix = "(All Tools Averaged)"

    if df.empty:
        logger.warning(f"No data for tool: {tool}")
        return

    # Parse coverage
    def parse_cov(c):
        try:
            return int(str(c).lower().rstrip('x'))
        except:
            return 0

    df['cov_num'] = df['coverage'].apply(parse_cov)

    # Group by coverage and average
    grouped = df.groupby('cov_num')[metrics].mean().reset_index()
    grouped = grouped.sort_values('cov_num')

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()

    colors = plt.cm.Set1.colors

    for idx, metric in enumerate(metrics[:4]):
        ax = axes[idx]
        color = colors[idx % len(colors)]

        ax.plot(grouped['cov_num'], grouped[metric],
                marker='o', color=color, lw=2, markersize=8)
        ax.fill_between(grouped['cov_num'], 0, grouped[metric],
                        alpha=0.2, color=color)

        ax.set_xlabel('Coverage (x)', fontsize=10)
        ax.set_ylabel(metric.upper(), fontsize=10)
        ax.set_title(f'{metric.upper()} vs Coverage', fontsize=12)
        ax.set_xscale('log')
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)

    plt.suptitle(f'Performance Metrics vs Coverage {title_suffix}', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved multi-metric coverage plot to {output_path}")


def generate_html_report(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    output_dir: Union[str, Path],
    reference: str = None
) -> str:
    """
    Generate comprehensive HTML report with all analyses.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Ground truth DataFrame
        output_dir: Directory to save report and plots
        reference: Filter to specific reference

    Returns:
        Path to HTML report
    """
    try:
        from .benchmark_metrics import calculate_metrics_all_tools
    except Exception:
        from benchmark_metrics import calculate_metrics_all_tools
    try:
        from .tool_comparison import summarize_tool_comparison
    except Exception:
        from tool_comparison import summarize_tool_comparison

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    plots = {}

    # PR curves
    pr_path = output_dir / 'pr_curves.png'
    plot_pr_curves(tool_outputs, ground_truth, pr_path, reference)
    plots['pr_curves'] = 'pr_curves.png'

    # ROC curves
    roc_path = output_dir / 'roc_curves.png'
    plot_roc_curves(tool_outputs, ground_truth, roc_path, reference)
    plots['roc_curves'] = 'roc_curves.png'

    # Score distributions
    dist_path = output_dir / 'score_distributions.png'
    plot_score_distributions(tool_outputs, dist_path, reference)
    plots['score_dist'] = 'score_distributions.png'

    # UpSet plot
    try:
        upset_path = output_dir / 'upset_plot.png'
        plot_upset(tool_outputs, upset_path, reference=reference)
        plots['upset'] = 'upset_plot.png'
    except Exception as e:
        logger.warning(f"Could not generate UpSet plot: {e}")

    # Correlation heatmap
    heatmap_path = output_dir / 'correlation_heatmap.png'
    plot_correlation_heatmap(tool_outputs, heatmap_path, reference=reference)
    plots['heatmap'] = 'correlation_heatmap.png'

    # Calculate metrics
    metrics_df = calculate_metrics_all_tools(tool_outputs, ground_truth, [reference] if reference else None)

    # Metrics comparison
    if not metrics_df.empty:
        metrics_path = output_dir / 'metrics_comparison.png'
        plot_metrics_comparison(metrics_df, metrics_path)
        plots['metrics'] = 'metrics_comparison.png'

    # Generate summary
    summary = summarize_tool_comparison(tool_outputs, ground_truth, reference=reference)

    # Build HTML
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>RNA Modification Detection - Downstream Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #333; border-bottom: 2px solid #4CAF50; padding-bottom: 10px; }}
        h2 {{ color: #555; margin-top: 30px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        img {{ max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #ddd; border-radius: 4px; }}
        .metric {{ font-weight: bold; color: #4CAF50; }}
        .summary-box {{ background: #e8f5e9; padding: 20px; border-radius: 4px; margin: 20px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>RNA Modification Detection - Downstream Analysis Report</h1>

        <div class="summary-box">
            <h2>Summary</h2>
            <p><strong>Tools analyzed:</strong> {', '.join(summary.get('tools', []))}</p>
            <p><strong>Total unique positions:</strong> {summary.get('n_union', 0)}</p>
            <p><strong>Positions detected by all tools:</strong> {summary.get('n_intersection', 0)}</p>
            <p><strong>Mean pairwise Jaccard:</strong> {summary.get('mean_pairwise_jaccard', 0):.3f}</p>
        </div>

        <h2>Benchmark Metrics</h2>
        {metrics_df.to_html(index=False, classes='metrics-table') if not metrics_df.empty else '<p>No metrics available</p>'}

        <h2>Precision-Recall Curves</h2>
        <img src="{plots.get('pr_curves', '')}" alt="PR Curves">

        <h2>ROC Curves</h2>
        <img src="{plots.get('roc_curves', '')}" alt="ROC Curves">

        <h2>Score Distributions</h2>
        <img src="{plots.get('score_dist', '')}" alt="Score Distributions">

        <h2>Tool Agreement</h2>
        {'<img src="' + plots.get('upset', '') + '" alt="UpSet Plot">' if 'upset' in plots else ''}

        <h2>Tool Correlation</h2>
        <img src="{plots.get('heatmap', '')}" alt="Correlation Heatmap">

        <footer style="margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; color: #666;">
            Generated by RNA Modification Detection Pipeline - Downstream Analysis Module
        </footer>
    </div>
</body>
</html>
    """

    report_path = output_dir / 'report.html'
    with open(report_path, 'w') as f:
        f.write(html_content)

    # Save metrics CSV
    if not metrics_df.empty:
        metrics_df.to_csv(output_dir / 'metrics.csv', index=False)

    logger.info(f"Generated HTML report at {report_path}")
    return str(report_path)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate visualizations for RNA modification detection'
    )
    parser.add_argument('--tools', nargs='+', required=True,
                        help='Tool names and their output files (format: tool:file)')
    parser.add_argument('--ground-truth', required=True,
                        help='Path to ground truth file')
    parser.add_argument('--output-dir', '-o', required=True,
                        help='Output directory for plots and report')
    parser.add_argument('--reference', help='Filter to specific reference')

    args = parser.parse_args()

    from parse_outputs import parse_tool_output
    from benchmark_metrics import load_ground_truth

    # Parse tool:file pairs
    tool_outputs = {}
    for item in args.tools:
        tool, filepath = item.split(':')
        tool_outputs[tool] = parse_tool_output(tool, filepath)

    # Load ground truth
    ground_truth = load_ground_truth(args.ground_truth)

    # Generate report
    report_path = generate_html_report(
        tool_outputs, ground_truth, args.output_dir, args.reference
    )
    print(f"Generated report: {report_path}")
