#!/usr/bin/env python3
"""
Benchmark Metrics Module for RNA Modification Detection

Calculates performance metrics for modification detection tools using
known modification sites as ground truth.

Metrics calculated:
- AUPRC (Area Under Precision-Recall Curve) - preferred for imbalanced data
- AUROC (Area Under ROC Curve)
- F1 Score at optimal threshold
- Precision and Recall at various thresholds
- Optimal threshold determination (maximizes F1 or Youden's J)

Usage:
    from benchmark_metrics import calculate_metrics, load_ground_truth

    ground_truth = load_ground_truth("known_sites.csv")
    results = calculate_metrics(tool_output, ground_truth)
"""

import logging
from pathlib import Path
from typing import Union, Dict, List, Tuple, Optional
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.metrics import (
    precision_recall_curve,
    average_precision_score,
    roc_curve,
    roc_auc_score,
    precision_recall_fscore_support,
    confusion_matrix
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class MetricsResult:
    """Container for benchmark metrics results."""
    tool: str
    reference: str
    replicate: str
    auprc: float
    auroc: float
    f1_optimal: float
    precision_optimal: float
    recall_optimal: float
    optimal_threshold: float
    n_true_positives: int
    n_false_positives: int
    n_false_negatives: int
    n_true_negatives: int
    n_ground_truth: int
    n_predictions: int

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'tool': self.tool,
            'reference': self.reference,
            'replicate': self.replicate,
            'auprc': self.auprc,
            'auroc': self.auroc,
            'f1_optimal': self.f1_optimal,
            'precision_optimal': self.precision_optimal,
            'recall_optimal': self.recall_optimal,
            'optimal_threshold': self.optimal_threshold,
            'n_true_positives': self.n_true_positives,
            'n_false_positives': self.n_false_positives,
            'n_false_negatives': self.n_false_negatives,
            'n_true_negatives': self.n_true_negatives,
            'n_ground_truth': self.n_ground_truth,
            'n_predictions': self.n_predictions,
        }


def load_ground_truth(
    filepath: Union[str, Path],
    reference_col: str = 'reference',
    position_col: str = 'position',
    mod_type_col: str = 'modification_type'
) -> pd.DataFrame:
    """
    Load ground truth modification sites from CSV/TSV file.

    Expected format:
        reference,position,modification_type,strand
        16S_rRNA,1402,m6A,+
        23S_rRNA,2445,Nm,+

    Args:
        filepath: Path to ground truth file
        reference_col: Column name for reference/chromosome
        position_col: Column name for position
        mod_type_col: Column name for modification type (optional)

    Returns:
        DataFrame with ground truth sites
    """
    filepath = Path(filepath)
    logger.info(f"Loading ground truth from: {filepath}")

    if not filepath.exists():
        raise FileNotFoundError(f"Ground truth file not found: {filepath}")

    # Detect separator
    sep = '\t' if filepath.suffix in ['.tsv', '.bed'] else ','

    df = pd.read_csv(filepath, sep=sep)

    # Standardize column names
    col_mapping = {}
    for col in df.columns:
        col_lower = col.lower()
        if 'ref' in col_lower or 'chr' in col_lower or 'contig' in col_lower:
            col_mapping[col] = 'reference'
        elif 'pos' in col_lower or 'start' in col_lower:
            col_mapping[col] = 'position'
        elif 'mod' in col_lower or 'type' in col_lower:
            col_mapping[col] = 'modification_type'
        elif 'strand' in col_lower:
            col_mapping[col] = 'strand'

    df = df.rename(columns=col_mapping)

    # Ensure required columns exist
    if 'reference' not in df.columns:
        df['reference'] = 'unknown'
    if 'position' not in df.columns:
        raise ValueError("Ground truth file must have a position column")

    df['position'] = df['position'].astype(int)

    logger.info(f"Loaded {len(df)} ground truth sites")
    return df


def prepare_labels(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame,
    reference: str = None,
    window: int = 0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Prepare binary labels and scores for metrics calculation.

    For each position in the tool output:
    - Label = 1 if position matches a ground truth site (within window)
    - Label = 0 otherwise

    Args:
        tool_output: Tool output DataFrame (must have 'position', 'score', 'reference')
        ground_truth: Ground truth DataFrame (must have 'position', 'reference')
        reference: Filter to specific reference (optional)
        window: Position matching window (0 = exact match)

    Returns:
        Tuple of (y_true, y_scores, positions)
    """
    # Filter by reference if specified
    if reference:
        tool_output = tool_output[
            tool_output['reference'].str.contains(reference, case=False, na=False)
        ].copy()
        ground_truth = ground_truth[
            ground_truth['reference'].str.contains(reference, case=False, na=False)
        ].copy()

    if tool_output.empty:
        logger.warning(f"No tool output for reference: {reference}")
        return np.array([]), np.array([]), np.array([])

    # Get ground truth positions as a set for fast lookup
    gt_positions = set(ground_truth['position'].values)

    # Create expanded set if using window matching
    if window > 0:
        gt_positions_expanded = set()
        for pos in gt_positions:
            for offset in range(-window, window + 1):
                gt_positions_expanded.add(pos + offset)
        gt_positions = gt_positions_expanded

    # Prepare labels
    positions = tool_output['position'].values
    scores = tool_output['score'].values

    # Handle different score types (lower p-value = more significant)
    score_type = tool_output['score_type'].iloc[0] if 'score_type' in tool_output.columns else 'pvalue'

    # For p-values/FDR, we need to invert (lower = better = higher score)
    if score_type in ['pvalue', 'fdr']:
        # Convert to -log10 scale, handling 0 values
        scores = np.where(scores > 0, -np.log10(scores + 1e-300), 300)
    elif score_type == 'zscore':
        # Absolute z-score (higher = more significant)
        scores = np.abs(scores)

    # Create binary labels
    y_true = np.array([1 if pos in gt_positions else 0 for pos in positions])

    # Handle NaN scores
    valid_mask = ~np.isnan(scores)
    y_true = y_true[valid_mask]
    scores = scores[valid_mask]
    positions = positions[valid_mask]

    return y_true, scores, positions


def _metrics_from_arrays(
    y_true: np.ndarray,
    y_scores: np.ndarray,
) -> Tuple[float, float, float, float, float, float, int, int, int, int]:
    """Compute AUPRC/AUROC/optimal F1 and confusion values."""
    if len(y_true) == 0:
        return 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0

    # AUPRC
    try:
        auprc = float(average_precision_score(y_true, y_scores))
    except Exception:
        auprc = 0.0

    # AUROC
    try:
        auroc = float(roc_auc_score(y_true, y_scores))
    except Exception:
        auroc = 0.5

    # optimal threshold from PR curve (max F1)
    try:
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)

        if len(thresholds) == 0:
            optimal_threshold = 0.0
            f1_optimal = 0.0
            precision_optimal = 0.0
            recall_optimal = 0.0
        else:
            optimal_idx = int(np.argmax(f1_scores[:-1]))
            optimal_threshold = float(thresholds[optimal_idx])
            f1_optimal = float(f1_scores[optimal_idx])
            precision_optimal = float(precision[optimal_idx])
            recall_optimal = float(recall[optimal_idx])
    except Exception:
        optimal_threshold = 0.0
        f1_optimal = 0.0
        precision_optimal = 0.0
        recall_optimal = 0.0

    y_pred = (y_scores >= optimal_threshold).astype(int)
    try:
        cm_flat = np.asarray(confusion_matrix(y_true, y_pred, labels=[0, 1]), dtype=int).reshape(-1)
        if cm_flat.size >= 4:
            tn, fp, fn, tp = (int(cm_flat[0]), int(cm_flat[1]), int(cm_flat[2]), int(cm_flat[3]))
        else:
            tn = fp = fn = tp = 0
    except Exception:
        tn = fp = fn = tp = 0

    return (
        auprc,
        auroc,
        f1_optimal,
        precision_optimal,
        recall_optimal,
        optimal_threshold,
        int(tp),
        int(fp),
        int(fn),
        int(tn),
    )


def calculate_metrics(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame,
    reference: str = None,
    window: int = 0
) -> MetricsResult:
    """
    Calculate benchmark metrics for a tool's output against ground truth.

    Args:
        tool_output: Tool output DataFrame
        ground_truth: Ground truth DataFrame
        reference: Filter to specific reference (optional)
        window: Position matching window (0 = exact match)

    Returns:
        MetricsResult with all calculated metrics
    """
    tool_name = tool_output['tool'].iloc[0] if 'tool' in tool_output.columns else 'unknown'
    replicate = tool_output['replicate'].iloc[0] if 'replicate' in tool_output.columns else 'unknown'

    logger.info(f"Calculating metrics for {tool_name} (ref={reference}, rep={replicate})")

    # Prepare labels
    y_true, y_scores, positions = prepare_labels(
        tool_output, ground_truth, reference, window
    )

    if len(y_true) == 0 or y_true.sum() == 0:
        logger.warning(f"No positive samples found for {tool_name}")
        return MetricsResult(
            tool=tool_name,
            reference=reference or 'all',
            replicate=replicate,
            auprc=0.0,
            auroc=0.5,
            f1_optimal=0.0,
            precision_optimal=0.0,
            recall_optimal=0.0,
            optimal_threshold=0.0,
            n_true_positives=0,
            n_false_positives=0,
            n_false_negatives=int(y_true.sum()),
            n_true_negatives=int((1 - y_true).sum()),
            n_ground_truth=len(ground_truth) if reference is None else len(
                ground_truth[ground_truth['reference'].str.contains(reference, case=False, na=False)]
            ),
            n_predictions=len(tool_output)
        )

    # Calculate AUPRC (Average Precision)
    try:
        auprc = average_precision_score(y_true, y_scores)
    except Exception as e:
        logger.warning(f"Could not calculate AUPRC: {e}")
        auprc = 0.0

    # Calculate AUROC
    try:
        auroc = roc_auc_score(y_true, y_scores)
    except Exception as e:
        logger.warning(f"Could not calculate AUROC: {e}")
        auroc = 0.5

    # Calculate precision-recall curve and find optimal threshold
    try:
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)

        # Calculate F1 at each threshold
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)

        # Find optimal threshold (maximizes F1)
        optimal_idx = np.argmax(f1_scores[:-1])  # Last element has no threshold
        optimal_threshold = thresholds[optimal_idx]
        f1_optimal = f1_scores[optimal_idx]
        precision_optimal = precision[optimal_idx]
        recall_optimal = recall[optimal_idx]

    except Exception as e:
        logger.warning(f"Could not calculate PR curve: {e}")
        optimal_threshold = 0.0
        f1_optimal = 0.0
        precision_optimal = 0.0
        recall_optimal = 0.0

    # Calculate confusion matrix at optimal threshold
    y_pred = (y_scores >= optimal_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    # Get ground truth count
    n_gt = len(ground_truth)
    if reference:
        gt_filtered = ground_truth[
            ground_truth['reference'].str.contains(reference, case=False, na=False)
        ]
        n_gt = len(gt_filtered)

    return MetricsResult(
        tool=tool_name,
        reference=reference or 'all',
        replicate=replicate,
        auprc=float(auprc),
        auroc=float(auroc),
        f1_optimal=float(f1_optimal),
        precision_optimal=float(precision_optimal),
        recall_optimal=float(recall_optimal),
        optimal_threshold=float(optimal_threshold),
        n_true_positives=int(tp),
        n_false_positives=int(fp),
        n_false_negatives=int(fn),
        n_true_negatives=int(tn),
        n_ground_truth=n_gt,
        n_predictions=len(tool_output)
    )


def calculate_metrics_all_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    references: List[str] = None,
    window: int = 0
) -> pd.DataFrame:
    """
    Calculate metrics for all tools and references.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Ground truth DataFrame
        references: List of references to analyze (default: all unique in ground truth)
        window: Position matching window

    Returns:
        DataFrame with metrics for all tool/reference combinations
    """
    if references is None:
        references = ground_truth['reference'].unique().tolist()
        references.append(None)  # Also calculate across all references

    results = []

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            logger.warning(f"Empty output for {tool_name}, skipping")
            continue

        for ref in references:
            try:
                metrics = calculate_metrics(tool_df, ground_truth, ref, window)
                results.append(metrics.to_dict())
            except Exception as e:
                logger.error(f"Error calculating metrics for {tool_name}/{ref}: {e}")

    return pd.DataFrame(results)


def get_precision_recall_curve(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame,
    reference: str = None,
    window: int = 0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get precision-recall curve data for plotting.

    Args:
        tool_output: Tool output DataFrame
        ground_truth: Ground truth DataFrame
        reference: Filter to specific reference (optional)
        window: Position matching window

    Returns:
        Tuple of (precision, recall, thresholds)
    """
    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0 or y_true.sum() == 0:
        return np.array([0, 1]), np.array([1, 0]), np.array([0])

    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    return precision, recall, thresholds


def get_roc_curve(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame,
    reference: str = None,
    window: int = 0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get ROC curve data for plotting.

    Args:
        tool_output: Tool output DataFrame
        ground_truth: Ground truth DataFrame
        reference: Filter to specific reference (optional)
        window: Position matching window

    Returns:
        Tuple of (fpr, tpr, thresholds)
    """
    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0 or y_true.sum() == 0:
        return np.array([0, 1]), np.array([0, 1]), np.array([0])

    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    return fpr, tpr, thresholds


def find_optimal_threshold(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame,
    reference: str = None,
    method: str = 'f1',
    window: int = 0
) -> Tuple[float, float]:
    """
    Find optimal score threshold for calling modifications.

    Args:
        tool_output: Tool output DataFrame
        ground_truth: Ground truth DataFrame
        reference: Filter to specific reference (optional)
        method: Optimization method ('f1' or 'youden')
        window: Position matching window

    Returns:
        Tuple of (optimal_threshold, metric_value)
    """
    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0 or y_true.sum() == 0:
        return 0.0, 0.0

    if method == 'f1':
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
        optimal_idx = np.argmax(f1_scores[:-1])
        return float(thresholds[optimal_idx]), float(f1_scores[optimal_idx])

    elif method == 'youden':
        # Youden's J statistic: J = sensitivity + specificity - 1 = TPR - FPR
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        j_scores = tpr - fpr
        optimal_idx = np.argmax(j_scores)
        return float(thresholds[optimal_idx]), float(j_scores[optimal_idx])

    else:
        raise ValueError(f"Unknown optimization method: {method}")


def calculate_metrics_at_threshold(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame,
    threshold: float,
    reference: str = None,
    window: int = 0
) -> Dict:
    """
    Calculate metrics at a specific threshold.

    Args:
        tool_output: Tool output DataFrame
        ground_truth: Ground truth DataFrame
        threshold: Score threshold for calling positives
        reference: Filter to specific reference (optional)
        window: Position matching window

    Returns:
        Dictionary with precision, recall, F1, etc. at threshold
    """
    y_true, y_scores, positions = prepare_labels(
        tool_output, ground_truth, reference, window
    )

    if len(y_true) == 0:
        return {
            'threshold': threshold,
            'precision': 0.0,
            'recall': 0.0,
            'f1': 0.0,
            'tp': 0,
            'fp': 0,
            'fn': 0,
            'tn': 0
        }

    y_pred = (y_scores >= threshold).astype(int)

    precision, recall, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, average='binary', zero_division=0
    )

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    return {
        'threshold': threshold,
        'precision': float(precision),
        'recall': float(recall),
        'f1': float(f1),
        'tp': int(tp),
        'fp': int(fp),
        'fn': int(fn),
        'tn': int(tn)
    }


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Calculate benchmark metrics for RNA modification detection'
    )
    parser.add_argument('tool_output', help='Path to tool output file')
    parser.add_argument('ground_truth', help='Path to ground truth file')
    parser.add_argument('--tool', required=True, help='Tool name')
    parser.add_argument('--reference', help='Filter to specific reference')
    parser.add_argument('--window', type=int, default=0, help='Position matching window')
    parser.add_argument('--output', '-o', help='Output CSV file for metrics')

    args = parser.parse_args()

    # Import parse function
    from parse_outputs import parse_tool_output

    # Load data
    tool_df = parse_tool_output(args.tool, args.tool_output)
    gt_df = load_ground_truth(args.ground_truth)

    # Calculate metrics
    metrics = calculate_metrics(tool_df, gt_df, args.reference, args.window)

    # Output
    if args.output:
        pd.DataFrame([metrics.to_dict()]).to_csv(args.output, index=False)
        print(f"Saved metrics to {args.output}")
    else:
        print("\nBenchmark Metrics:")
        print(f"  Tool: {metrics.tool}")
        print(f"  Reference: {metrics.reference}")
        print(f"  AUPRC: {metrics.auprc:.4f}")
        print(f"  AUROC: {metrics.auroc:.4f}")
        print(f"  F1 (optimal): {metrics.f1_optimal:.4f}")
        print(f"  Precision: {metrics.precision_optimal:.4f}")
        print(f"  Recall: {metrics.recall_optimal:.4f}")
        print(f"  Optimal threshold: {metrics.optimal_threshold:.4f}")
        print(f"  TP/FP/FN/TN: {metrics.n_true_positives}/{metrics.n_false_positives}/"
              f"{metrics.n_false_negatives}/{metrics.n_true_negatives}")
