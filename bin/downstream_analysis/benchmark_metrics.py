#!/usr/bin/env python3
"""
Benchmark metrics for RNA modification detection outputs.

Supports both:
- classic mode: tool outputs + ground truth matching
- labeled mode: precomputed label/score_metric columns from full-universe tables
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_recall_curve,
    precision_recall_fscore_support,
    roc_auc_score,
    roc_curve,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


@dataclass
class MetricsResult:
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
        return {
            "tool": self.tool,
            "reference": self.reference,
            "replicate": self.replicate,
            "auprc": self.auprc,
            "auroc": self.auroc,
            "f1_optimal": self.f1_optimal,
            "precision_optimal": self.precision_optimal,
            "recall_optimal": self.recall_optimal,
            "optimal_threshold": self.optimal_threshold,
            "n_true_positives": self.n_true_positives,
            "n_false_positives": self.n_false_positives,
            "n_false_negatives": self.n_false_negatives,
            "n_true_negatives": self.n_true_negatives,
            "n_ground_truth": self.n_ground_truth,
            "n_predictions": self.n_predictions,
        }


def _normalize_reference_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower()


def _filter_reference(df: pd.DataFrame, reference: Optional[str]) -> pd.DataFrame:
    if reference is None or df.empty or "reference" not in df.columns:
        return df

    ref_l = reference.strip().lower()
    return df[_normalize_reference_series(df["reference"]) == ref_l].copy()


def load_ground_truth(
    filepath: Union[str, Path],
    reference_col: str = "reference",
    position_col: str = "position",
    mod_type_col: str = "modification_type",
) -> pd.DataFrame:
    filepath = Path(filepath)
    logger.info("Loading ground truth from: %s", filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Ground truth file not found: {filepath}")

    sep = "\t" if filepath.suffix in {".tsv", ".bed"} else ","
    df = pd.read_csv(filepath, sep=sep)

    col_mapping = {}
    for col in df.columns:
        col_l = col.lower()
        if "ref" in col_l or "chr" in col_l or "contig" in col_l:
            col_mapping[col] = "reference"
        elif "pos" in col_l or "start" in col_l:
            col_mapping[col] = "position"
        elif "mod" in col_l or "type" in col_l:
            col_mapping[col] = "modification_type"
        elif "strand" in col_l:
            col_mapping[col] = "strand"

    df = df.rename(columns=col_mapping)

    if "reference" not in df.columns:
        df["reference"] = "unknown"
    if "position" not in df.columns:
        raise ValueError("Ground truth file must have a position column")

    df["reference"] = df["reference"].astype(str)
    df["position"] = pd.to_numeric(df["position"], errors="coerce").astype("Int64")
    df = df[df["position"].notna()].copy()
    df["position"] = df["position"].astype(int)

    if mod_type_col not in df.columns and "modification_type" not in df.columns:
        df["modification_type"] = "mod"

    logger.info("Loaded %d ground truth rows", len(df))
    return df


def _positions_with_window(gt_df: pd.DataFrame, window: int) -> set:
    if gt_df.empty:
        return set()

    if window <= 0:
        return set(zip(gt_df["reference"], gt_df["position"]))

    expanded = set()
    for ref, pos in zip(gt_df["reference"], gt_df["position"]):
        for offset in range(-window, window + 1):
            expanded.add((ref, int(pos) + offset))
    return expanded


def _transform_scores(df: pd.DataFrame) -> np.ndarray:
    """Return ranking scores where higher means more likely modified."""
    if "score_metric" in df.columns:
        return pd.to_numeric(df["score_metric"], errors="coerce").values

    scores = pd.to_numeric(df.get("score", np.nan), errors="coerce").values

    score_type = "unknown"
    if "score_type" in df.columns and not df["score_type"].empty:
        score_type = str(df["score_type"].iloc[0]).lower()

    if score_type in {"pvalue", "fdr"}:
        return np.where(scores > 0, -np.log10(np.clip(scores, 1e-300, 1.0)), np.nan)
    if score_type == "zscore":
        return np.abs(scores)
    return scores


def prepare_labels(
    tool_output: pd.DataFrame,
    ground_truth: Optional[pd.DataFrame],
    reference: str = None,
    window: int = 0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Prepare y_true/y_score with reference-aware matching.

    Returns positions as object array of (reference, position) tuples.
    """
    tool_df = _filter_reference(tool_output.copy(), reference)
    gt_df = _filter_reference(ground_truth.copy(), reference) if ground_truth is not None else pd.DataFrame()

    if tool_df.empty:
        return np.array([]), np.array([]), np.array([])

    tool_df = tool_df.copy()
    tool_df["reference"] = tool_df["reference"].astype(str)
    tool_df["position"] = pd.to_numeric(tool_df["position"], errors="coerce").astype("Int64")
    tool_df = tool_df[tool_df["position"].notna()].copy()
    tool_df["position"] = tool_df["position"].astype(int)

    scores = _transform_scores(tool_df)

    if "label" in tool_df.columns:
        y_true = pd.to_numeric(tool_df["label"], errors="coerce").fillna(0).astype(int).values
    else:
        if gt_df.empty:
            logger.warning("No ground truth provided for unlabeled tool outputs")
            return np.array([]), np.array([]), np.array([])

        gt_df = gt_df.copy()
        gt_df["reference"] = gt_df["reference"].astype(str)
        gt_df["position"] = pd.to_numeric(gt_df["position"], errors="coerce").astype("Int64")
        gt_df = gt_df[gt_df["position"].notna()].copy()
        gt_df["position"] = gt_df["position"].astype(int)

        gt_pos = _positions_with_window(gt_df, window)
        y_true = np.array(
            [
                1 if (ref, pos) in gt_pos else 0
                for ref, pos in zip(tool_df["reference"].values, tool_df["position"].values)
            ]
        )

    valid = ~np.isnan(scores)
    y_true = y_true[valid]
    scores = scores[valid]

    positions = np.array(
        list(
            zip(
                tool_df.loc[valid, "reference"].values,
                tool_df.loc[valid, "position"].astype(int).values,
            )
        ),
        dtype=object,
    )

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
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

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
    ground_truth: Optional[pd.DataFrame],
    reference: str = None,
    window: int = 0,
) -> MetricsResult:
    tool_name = tool_output["tool"].iloc[0] if "tool" in tool_output.columns and not tool_output.empty else "unknown"
    replicate = (
        tool_output["replicate"].iloc[0]
        if "replicate" in tool_output.columns and not tool_output.empty
        else "unknown"
    )

    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0:
        n_gt = 0
        if ground_truth is not None and not ground_truth.empty:
            gt_ref = _filter_reference(ground_truth, reference)
            n_gt = len(gt_ref)

        return MetricsResult(
            tool=tool_name,
            reference=reference or "all",
            replicate=replicate,
            auprc=0.0,
            auroc=0.5,
            f1_optimal=0.0,
            precision_optimal=0.0,
            recall_optimal=0.0,
            optimal_threshold=0.0,
            n_true_positives=0,
            n_false_positives=0,
            n_false_negatives=0,
            n_true_negatives=0,
            n_ground_truth=n_gt,
            n_predictions=len(tool_output),
        )

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

    n_gt = int(y_true.sum())
    if ground_truth is not None and not ground_truth.empty:
        n_gt = len(_filter_reference(ground_truth, reference)) if reference else len(ground_truth)

    return MetricsResult(
        tool=tool_name,
        reference=reference or "all",
        replicate=replicate,
        auprc=auprc,
        auroc=auroc,
        f1_optimal=f1_optimal,
        precision_optimal=precision_optimal,
        recall_optimal=recall_optimal,
        optimal_threshold=optimal_threshold,
        n_true_positives=tp,
        n_false_positives=fp,
        n_false_negatives=fn,
        n_true_negatives=tn,
        n_ground_truth=n_gt,
        n_predictions=len(tool_output),
    )


def calculate_metrics_all_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: Optional[pd.DataFrame],
    references: List[str] = None,
    window: int = 0,
) -> pd.DataFrame:
    if references is None:
        if ground_truth is not None and not ground_truth.empty and "reference" in ground_truth.columns:
            references = sorted(ground_truth["reference"].dropna().astype(str).unique().tolist())
        else:
            refs = []
            for df in tool_outputs.values():
                if not df.empty and "reference" in df.columns:
                    refs.extend(df["reference"].dropna().astype(str).unique().tolist())
            references = sorted(set(refs))

    rows = []
    for tool_name, df in tool_outputs.items():
        if df.empty:
            continue
        for ref in references:
            try:
                rows.append(calculate_metrics(df, ground_truth, reference=ref, window=window).to_dict())
            except Exception as exc:
                logger.error("Failed metrics for %s/%s: %s", tool_name, ref, exc)

    return pd.DataFrame(rows)


def get_precision_recall_curve(
    tool_output: pd.DataFrame,
    ground_truth: Optional[pd.DataFrame],
    reference: str = None,
    window: int = 0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0 or y_true.sum() == 0:
        return np.array([0.0, 1.0]), np.array([1.0, 0.0]), np.array([0.0])

    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    return precision, recall, thresholds


def get_roc_curve(
    tool_output: pd.DataFrame,
    ground_truth: Optional[pd.DataFrame],
    reference: str = None,
    window: int = 0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0 or len(np.unique(y_true)) < 2:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([0.0])

    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    return fpr, tpr, thresholds


def find_optimal_threshold(
    tool_output: pd.DataFrame,
    ground_truth: Optional[pd.DataFrame],
    reference: str = None,
    method: str = "f1",
    window: int = 0,
) -> Tuple[float, float]:
    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0:
        return 0.0, 0.0

    if method == "f1":
        precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
        if len(thresholds) == 0:
            return 0.0, 0.0
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
        idx = int(np.argmax(f1_scores[:-1]))
        return float(thresholds[idx]), float(f1_scores[idx])

    if method == "youden":
        fpr, tpr, thresholds = roc_curve(y_true, y_scores)
        idx = int(np.argmax(tpr - fpr))
        return float(thresholds[idx]), float((tpr - fpr)[idx])

    raise ValueError(f"Unknown optimization method: {method}")


def calculate_metrics_at_threshold(
    tool_output: pd.DataFrame,
    ground_truth: Optional[pd.DataFrame],
    threshold: float,
    reference: str = None,
    window: int = 0,
) -> Dict:
    y_true, y_scores, _ = prepare_labels(tool_output, ground_truth, reference, window)

    if len(y_true) == 0:
        return {
            "threshold": threshold,
            "precision": 0.0,
            "recall": 0.0,
            "f1": 0.0,
            "tp": 0,
            "fp": 0,
            "fn": 0,
            "tn": 0,
        }

    y_pred = (y_scores >= threshold).astype(int)
    precision, recall, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, average="binary", zero_division=0
    )
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    return {
        "threshold": threshold,
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "tn": int(tn),
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark metrics")
    parser.add_argument("tool_output")
    parser.add_argument("ground_truth")
    parser.add_argument("--tool", required=True)
    parser.add_argument("--reference")
    parser.add_argument("--window", type=int, default=0)
    parser.add_argument("--output", "-o")
    args = parser.parse_args()

    from parse_outputs import parse_tool_output

    tool_df = parse_tool_output(args.tool, args.tool_output)
    gt_df = load_ground_truth(args.ground_truth)
    metrics = calculate_metrics(tool_df, gt_df, args.reference, args.window)

    if args.output:
        pd.DataFrame([metrics.to_dict()]).to_csv(args.output, index=False)
    else:
        print(pd.DataFrame([metrics.to_dict()]).to_string(index=False))
