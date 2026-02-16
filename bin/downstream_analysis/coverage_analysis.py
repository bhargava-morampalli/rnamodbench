#!/usr/bin/env python3
"""
Coverage Analysis Module for RNA Modification Detection Pipeline

Analyzes how tool performance varies with sequencing coverage level.
Supports subsampled coverage experiments where the same samples are
analyzed at different coverage depths.

Key analyses:
- Per-coverage metrics calculation
- Performance vs coverage curves
- Saturation point detection
- Tool ranking by coverage stability

Expected data format:
    DataFrame with 'coverage' column containing values like "100x", "500x"
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from pathlib import Path

import numpy as np
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class CoverageMetrics:
    """Metrics for a single tool at a single coverage level."""
    tool: str
    coverage: str
    coverage_numeric: int
    n_positions: int
    auprc: float = np.nan
    auroc: float = np.nan
    f1: float = np.nan
    precision: float = np.nan
    recall: float = np.nan
    n_true_positives: int = 0
    n_false_positives: int = 0
    n_false_negatives: int = 0

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            'tool': self.tool,
            'coverage': self.coverage,
            'coverage_numeric': self.coverage_numeric,
            'n_positions': self.n_positions,
            'auprc': self.auprc,
            'auroc': self.auroc,
            'f1': self.f1,
            'precision': self.precision,
            'recall': self.recall,
            'n_true_positives': self.n_true_positives,
            'n_false_positives': self.n_false_positives,
            'n_false_negatives': self.n_false_negatives,
        }


@dataclass
class CoverageAnalysisResult:
    """Complete coverage analysis results for all tools."""
    metrics_by_coverage: List[CoverageMetrics]
    saturation_points: Dict[str, Optional[str]]  # tool -> coverage at saturation
    coverage_stability: Dict[str, float]  # tool -> stability score
    coverage_levels: List[str]

    def to_dataframe(self) -> pd.DataFrame:
        """Convert metrics to DataFrame."""
        return pd.DataFrame([m.to_dict() for m in self.metrics_by_coverage])

    def save(self, output_path: Path) -> None:
        """Save results to CSV."""
        self.to_dataframe().to_csv(output_path, index=False)
        logger.info(f"Saved coverage analysis to {output_path}")


def parse_coverage_value(coverage: str) -> int:
    """
    Parse coverage string to numeric value.

    Args:
        coverage: Coverage string like "100x", "500X", "1000"

    Returns:
        Numeric coverage value
    """
    if pd.isna(coverage):
        return 0

    # Remove 'x' or 'X' suffix and convert to int
    coverage_str = str(coverage).lower().rstrip('x')
    try:
        return int(coverage_str)
    except ValueError:
        logger.warning(f"Could not parse coverage value: {coverage}")
        return 0


def get_coverage_levels(df: pd.DataFrame) -> List[str]:
    """
    Get sorted list of coverage levels from DataFrame.

    Args:
        df: DataFrame with 'coverage' column

    Returns:
        List of coverage levels sorted numerically
    """
    if 'coverage' not in df.columns:
        return []

    coverages = df['coverage'].dropna().unique()
    # Sort numerically
    return sorted(coverages, key=parse_coverage_value)


def calculate_metrics_by_coverage(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    references: List[str] = None
) -> List[CoverageMetrics]:
    """
    Calculate performance metrics for each tool at each coverage level.

    Args:
        tool_outputs: Dictionary mapping tool names to DataFrames
        ground_truth: Ground truth DataFrame with known modification sites
        references: Optional list of references to filter by

    Returns:
        List of CoverageMetrics objects
    """
    from benchmark_metrics import calculate_metrics

    results = []

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty or 'coverage' not in tool_df.columns:
            continue

        coverage_levels = get_coverage_levels(tool_df)

        for coverage in coverage_levels:
            # Filter to this coverage level
            cov_df = tool_df[tool_df['coverage'] == coverage].copy()

            if cov_df.empty:
                continue

            # Filter by references if specified
            if references:
                cov_df = cov_df[cov_df['reference'].isin(references)]

            if cov_df.empty:
                continue

            # Calculate metrics
            try:
                metrics_result = calculate_metrics(
                    cov_df, ground_truth, reference=None
                )

                cov_metrics = CoverageMetrics(
                    tool=tool_name,
                    coverage=coverage,
                    coverage_numeric=parse_coverage_value(coverage),
                    n_positions=len(cov_df),
                    auprc=metrics_result.auprc,
                    auroc=metrics_result.auroc,
                    f1=metrics_result.f1_optimal,
                    precision=metrics_result.precision_optimal,
                    recall=metrics_result.recall_optimal,
                    n_true_positives=metrics_result.n_true_positives,
                    n_false_positives=metrics_result.n_false_positives,
                    n_false_negatives=metrics_result.n_false_negatives,
                )
                results.append(cov_metrics)

            except Exception as e:
                logger.warning(f"Failed to calculate metrics for {tool_name} at {coverage}: {e}")
                # Add empty metrics
                results.append(CoverageMetrics(
                    tool=tool_name,
                    coverage=coverage,
                    coverage_numeric=parse_coverage_value(coverage),
                    n_positions=len(cov_df),
                ))

    return results


def detect_saturation_point(
    metrics: List[CoverageMetrics],
    tool: str,
    metric_name: str = 'auprc',
    threshold: float = 0.01
) -> Optional[str]:
    """
    Detect the coverage level at which performance saturates.

    Saturation is defined as the point where increasing coverage
    yields less than `threshold` relative improvement.

    Args:
        metrics: List of CoverageMetrics
        tool: Tool name to analyze
        metric_name: Metric to use (auprc, f1, etc.)
        threshold: Relative improvement threshold for saturation

    Returns:
        Coverage level at saturation, or None if not reached
    """
    # Filter metrics for this tool and sort by coverage
    tool_metrics = [m for m in metrics if m.tool == tool]
    tool_metrics.sort(key=lambda m: m.coverage_numeric)

    if len(tool_metrics) < 2:
        return None

    prev_value = None
    for i, m in enumerate(tool_metrics):
        current_value = getattr(m, metric_name)

        if pd.isna(current_value):
            continue

        if prev_value is not None:
            # Calculate relative improvement
            if prev_value > 0:
                improvement = (current_value - prev_value) / prev_value
            else:
                improvement = float('inf') if current_value > 0 else 0

            # Check if improvement is below threshold
            if improvement < threshold:
                return m.coverage

        prev_value = current_value

    return None  # Saturation not reached


def calculate_coverage_stability(
    metrics: List[CoverageMetrics],
    tool: str,
    metric_name: str = 'auprc'
) -> float:
    """
    Calculate stability score for a tool across coverage levels.

    Stability is measured as the coefficient of variation (CV) of
    the metric across coverage levels. Lower CV = more stable.

    Args:
        metrics: List of CoverageMetrics
        tool: Tool name to analyze
        metric_name: Metric to use

    Returns:
        Stability score (1 - normalized CV, higher is better)
    """
    tool_metrics = [m for m in metrics if m.tool == tool]
    values = [getattr(m, metric_name) for m in tool_metrics
              if not pd.isna(getattr(m, metric_name))]

    if len(values) < 2:
        return np.nan

    mean_val = np.mean(values)
    std_val = np.std(values)

    if mean_val == 0:
        return 0.0

    cv = std_val / mean_val
    # Convert to stability score (0-1, higher is better)
    stability = 1.0 / (1.0 + cv)

    return stability


def analyze_coverage(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    references: List[str] = None
) -> CoverageAnalysisResult:
    """
    Perform complete coverage analysis.

    Args:
        tool_outputs: Dictionary mapping tool names to DataFrames
        ground_truth: Ground truth DataFrame
        references: Optional list of references to filter by

    Returns:
        CoverageAnalysisResult with all analyses
    """
    # Calculate metrics for each tool at each coverage
    metrics = calculate_metrics_by_coverage(tool_outputs, ground_truth, references)

    if not metrics:
        logger.warning("No coverage metrics calculated")
        return CoverageAnalysisResult(
            metrics_by_coverage=[],
            saturation_points={},
            coverage_stability={},
            coverage_levels=[]
        )

    # Get unique tools and coverage levels
    tools = list(set(m.tool for m in metrics))
    coverage_levels = list(set(m.coverage for m in metrics))
    coverage_levels.sort(key=parse_coverage_value)

    # Detect saturation points
    saturation_points = {}
    for tool in tools:
        saturation_points[tool] = detect_saturation_point(metrics, tool)

    # Calculate stability scores
    stability_scores = {}
    for tool in tools:
        stability_scores[tool] = calculate_coverage_stability(metrics, tool)

    return CoverageAnalysisResult(
        metrics_by_coverage=metrics,
        saturation_points=saturation_points,
        coverage_stability=stability_scores,
        coverage_levels=coverage_levels
    )


def compare_tools_by_coverage(
    coverage_result: CoverageAnalysisResult,
    metric_name: str = 'auprc'
) -> pd.DataFrame:
    """
    Create a pivot table comparing tools across coverage levels.

    Args:
        coverage_result: Coverage analysis results
        metric_name: Metric to compare

    Returns:
        DataFrame with tools as rows, coverages as columns
    """
    df = coverage_result.to_dataframe()

    if df.empty:
        return pd.DataFrame()

    # Create pivot table
    pivot = df.pivot_table(
        index='tool',
        columns='coverage',
        values=metric_name,
        aggfunc='mean'
    )

    # Sort columns by coverage value
    sorted_cols = sorted(pivot.columns, key=parse_coverage_value)
    pivot = pivot[sorted_cols]

    return pivot


def get_optimal_coverage(
    coverage_result: CoverageAnalysisResult,
    cost_weight: float = 0.5
) -> Dict[str, str]:
    """
    Suggest optimal coverage level for each tool based on
    performance/cost tradeoff.

    Args:
        coverage_result: Coverage analysis results
        cost_weight: Weight for coverage cost (0 = ignore cost, 1 = prioritize low cost)

    Returns:
        Dictionary mapping tool names to recommended coverage levels
    """
    recommendations = {}

    df = coverage_result.to_dataframe()
    if df.empty:
        return recommendations

    for tool in df['tool'].unique():
        tool_df = df[df['tool'] == tool].copy()

        if tool_df.empty:
            continue

        # Normalize AUPRC and coverage to 0-1 scale
        tool_df['auprc_norm'] = (tool_df['auprc'] - tool_df['auprc'].min()) / \
                                (tool_df['auprc'].max() - tool_df['auprc'].min() + 1e-10)
        tool_df['cov_norm'] = (tool_df['coverage_numeric'] - tool_df['coverage_numeric'].min()) / \
                              (tool_df['coverage_numeric'].max() - tool_df['coverage_numeric'].min() + 1e-10)

        # Calculate composite score (higher AUPRC, lower coverage is better)
        tool_df['score'] = (1 - cost_weight) * tool_df['auprc_norm'] - \
                           cost_weight * tool_df['cov_norm']

        # Get coverage with best score
        best_idx = tool_df['score'].idxmax()
        recommendations[tool] = tool_df.loc[best_idx, 'coverage']

    return recommendations


def generate_coverage_summary(coverage_result: CoverageAnalysisResult) -> str:
    """
    Generate a text summary of coverage analysis.

    Args:
        coverage_result: Coverage analysis results

    Returns:
        Formatted summary string
    """
    lines = [
        "=" * 60,
        "COVERAGE ANALYSIS SUMMARY",
        "=" * 60,
        "",
        f"Coverage levels analyzed: {', '.join(coverage_result.coverage_levels)}",
        "",
    ]

    # Saturation points
    lines.append("SATURATION POINTS (coverage where performance plateaus):")
    lines.append("-" * 40)
    for tool, saturation in sorted(coverage_result.saturation_points.items()):
        if saturation:
            lines.append(f"  {tool}: {saturation}")
        else:
            lines.append(f"  {tool}: Not reached")

    lines.append("")

    # Stability scores
    lines.append("COVERAGE STABILITY (higher = more consistent across coverages):")
    lines.append("-" * 40)
    for tool, stability in sorted(coverage_result.coverage_stability.items(),
                                   key=lambda x: x[1] if not pd.isna(x[1]) else 0,
                                   reverse=True):
        if not pd.isna(stability):
            lines.append(f"  {tool}: {stability:.3f}")
        else:
            lines.append(f"  {tool}: N/A")

    lines.append("")

    # Best performance by coverage
    df = coverage_result.to_dataframe()
    if not df.empty:
        lines.append("BEST AUPRC BY COVERAGE LEVEL:")
        lines.append("-" * 40)
        for coverage in coverage_result.coverage_levels:
            cov_df = df[df['coverage'] == coverage]
            if not cov_df.empty:
                best_idx = cov_df['auprc'].idxmax()
                best_row = cov_df.loc[best_idx]
                lines.append(f"  {coverage}: {best_row['tool']} (AUPRC={best_row['auprc']:.4f})")

    lines.append("")
    lines.append("=" * 60)

    return '\n'.join(lines)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze tool performance across coverage levels'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input directory with coverage subdirectories')
    parser.add_argument('--ground-truth', '-g', required=True,
                        help='Ground truth file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')

    args = parser.parse_args()

    from parse_outputs import load_all_tool_outputs
    from benchmark_metrics import load_ground_truth

    # Load data
    tool_outputs = load_all_tool_outputs(args.input, coverage_dirs=True)
    ground_truth = load_ground_truth(args.ground_truth)

    # Run analysis
    result = analyze_coverage(tool_outputs, ground_truth)

    # Save results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    result.save(output_dir / 'coverage_metrics.csv')

    # Save pivot table
    pivot = compare_tools_by_coverage(result)
    if not pivot.empty:
        pivot.to_csv(output_dir / 'coverage_comparison.csv')

    # Save summary
    summary = generate_coverage_summary(result)
    with open(output_dir / 'coverage_summary.txt', 'w') as f:
        f.write(summary)
    print(summary)
