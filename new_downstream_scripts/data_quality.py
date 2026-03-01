#!/usr/bin/env python3
"""
Data Quality Module for RNA Modification Detection

Validates tool outputs and generates quality reports including:
- Tool/replicate status tracking
- Empty file detection
- Missing replicate warnings
- Position overlap analysis
- Data completeness metrics

Usage:
    from data_quality import DataQualityReport, validate_all_outputs

    report = validate_all_outputs(tool_outputs)
    report.save("data_quality_report.csv")
    report.save_warnings("analysis_warnings.txt")
"""

import logging
from pathlib import Path
from typing import Union, Dict, List, Optional, Set, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Status(Enum):
    """Status codes for tool outputs."""
    OK = "OK"
    EMPTY = "EMPTY"
    FAILED = "FAILED"
    PARTIAL = "PARTIAL"
    MISSING = "MISSING"


class WarningLevel(Enum):
    """Warning severity levels."""
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


@dataclass
class ToolReplicateStatus:
    """Status for a single tool-replicate combination."""
    tool: str
    replicate: str
    reference: str = "all"
    coverage: str = "full"
    n_positions: int = 0
    status: Status = Status.MISSING
    completeness: float = 0.0
    warnings: List[str] = field(default_factory=list)
    parse_error: Optional[str] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame."""
        return {
            'tool': self.tool,
            'replicate': self.replicate,
            'reference': self.reference,
            'coverage': self.coverage,
            'n_positions': self.n_positions,
            'status': self.status.value,
            'completeness': f"{self.completeness:.1%}",
            'warnings': '; '.join(self.warnings) if self.warnings else '',
            'parse_error': self.parse_error or ''
        }


@dataclass
class Warning:
    """A single warning message."""
    level: WarningLevel
    tool: str
    message: str
    details: Optional[str] = None

    def __str__(self) -> str:
        detail_str = f" ({self.details})" if self.details else ""
        return f"{self.level.value}: [{self.tool}] {self.message}{detail_str}"


@dataclass
class DataQualityReport:
    """
    Comprehensive data quality report for downstream analysis.

    Tracks:
    - Status of each tool/replicate combination
    - Warnings and errors
    - Position coverage statistics
    - Replicate completeness
    """
    tool_statuses: List[ToolReplicateStatus] = field(default_factory=list)
    warnings: List[Warning] = field(default_factory=list)
    expected_replicates: int = 3
    expected_tools: List[str] = field(default_factory=list)
    total_positions_union: int = 0
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def add_status(self, status: ToolReplicateStatus) -> None:
        """Add a tool-replicate status entry."""
        self.tool_statuses.append(status)

    def add_warning(
        self,
        level: WarningLevel,
        tool: str,
        message: str,
        details: str = None
    ) -> None:
        """Add a warning message."""
        warning = Warning(level, tool, message, details)
        self.warnings.append(warning)
        # Also log it
        if level == WarningLevel.ERROR:
            logger.error(str(warning))
        elif level == WarningLevel.WARNING:
            logger.warning(str(warning))
        else:
            logger.info(str(warning))

    def get_tools_with_data(self) -> List[str]:
        """Get list of tools that have at least some data."""
        return list(set(
            s.tool for s in self.tool_statuses
            if s.status in [Status.OK, Status.PARTIAL]
        ))

    def get_complete_tools(self) -> List[str]:
        """Get list of tools with all expected replicates."""
        tool_rep_counts = {}
        for s in self.tool_statuses:
            if s.status == Status.OK:
                tool_rep_counts[s.tool] = tool_rep_counts.get(s.tool, 0) + 1

        return [
            tool for tool, count in tool_rep_counts.items()
            if count >= self.expected_replicates
        ]

    def get_incomplete_tools(self) -> Dict[str, int]:
        """Get tools with fewer than expected replicates."""
        tool_rep_counts = {}
        for s in self.tool_statuses:
            if s.status == Status.OK:
                tool_rep_counts[s.tool] = tool_rep_counts.get(s.tool, 0) + 1

        return {
            tool: count for tool, count in tool_rep_counts.items()
            if count < self.expected_replicates
        }

    def to_dataframe(self) -> pd.DataFrame:
        """Convert status entries to DataFrame."""
        if not self.tool_statuses:
            return pd.DataFrame(columns=[
                'tool', 'replicate', 'reference', 'coverage',
                'n_positions', 'status', 'completeness', 'warnings', 'parse_error'
            ])
        return pd.DataFrame([s.to_dict() for s in self.tool_statuses])

    def save(self, filepath: Union[str, Path]) -> None:
        """Save quality report to CSV."""
        filepath = Path(filepath)
        df = self.to_dataframe()
        df.to_csv(filepath, index=False)
        logger.info(f"Saved data quality report to {filepath}")

    def save_warnings(self, filepath: Union[str, Path]) -> None:
        """Save warnings to text file."""
        filepath = Path(filepath)
        with open(filepath, 'w') as f:
            f.write(f"# Data Quality Warnings\n")
            f.write(f"# Generated: {self.timestamp}\n")
            f.write(f"# Expected replicates: {self.expected_replicates}\n")
            f.write(f"# Expected tools: {', '.join(self.expected_tools)}\n")
            f.write(f"#\n\n")

            # Group by level
            errors = [w for w in self.warnings if w.level == WarningLevel.ERROR]
            warnings = [w for w in self.warnings if w.level == WarningLevel.WARNING]
            infos = [w for w in self.warnings if w.level == WarningLevel.INFO]

            if errors:
                f.write("## ERRORS\n")
                for w in errors:
                    f.write(f"  - {w}\n")
                f.write("\n")

            if warnings:
                f.write("## WARNINGS\n")
                for w in warnings:
                    f.write(f"  - {w}\n")
                f.write("\n")

            if infos:
                f.write("## INFO\n")
                for w in infos:
                    f.write(f"  - {w}\n")
                f.write("\n")

            # Summary
            f.write("## SUMMARY\n")
            f.write(f"  - Total tools with data: {len(self.get_tools_with_data())}\n")
            f.write(f"  - Complete tools (all replicates): {len(self.get_complete_tools())}\n")
            incomplete = self.get_incomplete_tools()
            if incomplete:
                f.write(f"  - Incomplete tools: {incomplete}\n")
            f.write(f"  - Total positions (union): {self.total_positions_union}\n")

        logger.info(f"Saved warnings to {filepath}")

    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            "=" * 60,
            "DATA QUALITY SUMMARY",
            "=" * 60,
            f"Timestamp: {self.timestamp}",
            f"Expected replicates: {self.expected_replicates}",
            f"Tools with data: {len(self.get_tools_with_data())}",
            f"Complete tools: {len(self.get_complete_tools())}",
            f"Total positions (union): {self.total_positions_union}",
            "",
            f"Errors: {len([w for w in self.warnings if w.level == WarningLevel.ERROR])}",
            f"Warnings: {len([w for w in self.warnings if w.level == WarningLevel.WARNING])}",
            f"Info: {len([w for w in self.warnings if w.level == WarningLevel.INFO])}",
            "=" * 60,
        ]
        return '\n'.join(lines)


def validate_tool_output(
    df: pd.DataFrame,
    tool: str,
    replicate: str = "unknown",
    reference: str = "all",
    coverage: str = "full"
) -> ToolReplicateStatus:
    """
    Validate a single tool output DataFrame.

    Args:
        df: Tool output DataFrame
        tool: Tool name
        replicate: Replicate identifier
        reference: Reference name
        coverage: Coverage level

    Returns:
        ToolReplicateStatus with validation results
    """
    status = ToolReplicateStatus(
        tool=tool,
        replicate=replicate,
        reference=reference,
        coverage=coverage
    )

    # Check if DataFrame is empty
    if df is None or df.empty:
        status.status = Status.EMPTY
        status.n_positions = 0
        status.completeness = 0.0
        status.warnings.append("No output data")
        return status

    # Check for required columns
    required_cols = ['position', 'score']
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        status.status = Status.FAILED
        status.parse_error = f"Missing columns: {missing_cols}"
        return status

    # Count positions
    status.n_positions = len(df)

    # Check for NaN scores
    nan_count = df['score'].isna().sum()
    if nan_count > 0:
        nan_pct = nan_count / len(df) * 100
        status.warnings.append(f"{nan_count} positions ({nan_pct:.1f}%) have NaN scores")

    # Check score distribution
    if 'score_type' in df.columns:
        score_type = df['score_type'].iloc[0]
        if score_type in ['pvalue', 'fdr']:
            # Check for invalid p-values
            invalid = ((df['score'] < 0) | (df['score'] > 1)).sum()
            if invalid > 0:
                status.warnings.append(f"{invalid} positions have invalid p-values (not in 0-1)")

    # Set status based on findings
    if status.warnings:
        status.status = Status.PARTIAL
    else:
        status.status = Status.OK

    status.completeness = 1.0 if status.n_positions > 0 else 0.0

    return status


def check_replicate_completeness(
    tool_outputs: Dict[str, pd.DataFrame],
    expected_replicates: int = 3,
    report: DataQualityReport = None
) -> DataQualityReport:
    """
    Check that all tools have the expected number of replicates.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        expected_replicates: Expected number of replicates
        report: Existing report to update (or create new)

    Returns:
        Updated DataQualityReport
    """
    if report is None:
        report = DataQualityReport(expected_replicates=expected_replicates)

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            report.add_warning(
                WarningLevel.ERROR,
                tool_name,
                "No data available",
                "Tool output is empty"
            )
            continue

        # Count replicates
        if 'replicate' in tool_df.columns:
            replicates = tool_df['replicate'].unique()
            n_replicates = len(replicates)
        else:
            n_replicates = 1
            replicates = ['unknown']

        # Validate each replicate
        for rep in replicates:
            if 'replicate' in tool_df.columns:
                rep_df = tool_df[tool_df['replicate'] == rep]
            else:
                rep_df = tool_df

            status = validate_tool_output(rep_df, tool_name, rep)
            report.add_status(status)

        # Check replicate count
        if n_replicates < expected_replicates:
            report.add_warning(
                WarningLevel.WARNING,
                tool_name,
                f"Only {n_replicates}/{expected_replicates} replicates have data",
                f"Found replicates: {list(replicates)}"
            )
        elif n_replicates > expected_replicates:
            report.add_warning(
                WarningLevel.INFO,
                tool_name,
                f"More replicates than expected ({n_replicates}/{expected_replicates})",
                f"Found replicates: {list(replicates)}"
            )

    return report


def check_replicate_overlap(
    tool_df: pd.DataFrame,
    tool_name: str,
    min_jaccard: float = 0.5,
    report: DataQualityReport = None
) -> DataQualityReport:
    """
    Check position overlap between replicates of the same tool.

    Args:
        tool_df: Tool output DataFrame with 'replicate' column
        tool_name: Tool name
        min_jaccard: Minimum acceptable Jaccard similarity
        report: Existing report to update

    Returns:
        Updated DataQualityReport
    """
    if report is None:
        report = DataQualityReport()

    if 'replicate' not in tool_df.columns:
        return report

    replicates = tool_df['replicate'].unique()
    if len(replicates) < 2:
        return report

    # Get positions per replicate
    rep_positions = {}
    for rep in replicates:
        rep_df = tool_df[tool_df['replicate'] == rep].copy()
        rep_df['position'] = pd.to_numeric(rep_df['position'], errors='coerce')
        rep_df = rep_df[rep_df['position'].notna()]
        positions = set(zip(rep_df['reference'].astype(str), rep_df['position'].astype(int)))
        rep_positions[rep] = positions

    # Calculate pairwise Jaccard
    from itertools import combinations
    jaccard_values = []

    for rep1, rep2 in combinations(replicates, 2):
        pos1, pos2 = rep_positions[rep1], rep_positions[rep2]
        intersection = len(pos1 & pos2)
        union = len(pos1 | pos2)
        jaccard = intersection / union if union > 0 else 0.0
        jaccard_values.append(jaccard)

        if jaccard < min_jaccard:
            report.add_warning(
                WarningLevel.WARNING,
                tool_name,
                f"Low replicate overlap between {rep1} and {rep2}",
                f"Jaccard={jaccard:.3f}, shared={intersection}, union={union}"
            )

    # Report mean Jaccard
    if jaccard_values:
        mean_jaccard = np.mean(jaccard_values)
        if mean_jaccard < min_jaccard:
            report.add_warning(
                WarningLevel.WARNING,
                tool_name,
                f"Low overall replicate reproducibility",
                f"Mean Jaccard={mean_jaccard:.3f}"
            )
        else:
            report.add_warning(
                WarningLevel.INFO,
                tool_name,
                f"Replicate reproducibility",
                f"Mean Jaccard={mean_jaccard:.3f}"
            )

    return report


def check_position_coverage(
    tool_outputs: Dict[str, pd.DataFrame],
    reference_positions: Optional[Set[Tuple[str, int]]] = None,
    report: DataQualityReport = None
) -> DataQualityReport:
    """
    Check position coverage across tools.

    Args:
        tool_outputs: Dictionary mapping tool names to DataFrames
        reference_positions: Optional set of all possible positions
        report: Existing report to update

    Returns:
        Updated DataQualityReport
    """
    if report is None:
        report = DataQualityReport()

    # Calculate union of all positions (testable universe)
    all_positions: Set[Tuple[str, int]] = set()
    tool_positions = {}

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            tool_positions[tool_name] = set()
            continue

        tmp = tool_df.copy()
        tmp['position'] = pd.to_numeric(tmp['position'], errors='coerce')
        tmp = tmp[tmp['position'].notna()]
        positions = set(zip(tmp['reference'].astype(str), tmp['position'].astype(int)))
        tool_positions[tool_name] = positions
        all_positions.update(positions)

    report.total_positions_union = len(all_positions)

    # If reference positions provided, use those as universe
    if reference_positions:
        universe = reference_positions
    else:
        universe = all_positions

    universe_size = len(universe)

    # Calculate coverage per tool
    for tool_name, positions in tool_positions.items():
        if not positions:
            continue

        coverage_pct = len(positions) / universe_size * 100 if universe_size > 0 else 0

        # Update completeness in existing statuses
        for status in report.tool_statuses:
            if status.tool == tool_name:
                status.completeness = coverage_pct / 100

        # Report low coverage
        if coverage_pct < 50:
            report.add_warning(
                WarningLevel.INFO,
                tool_name,
                f"Reports {len(positions)}/{universe_size} positions ({coverage_pct:.1f}%)",
                "May only report significant positions"
            )

    return report


def validate_all_outputs(
    tool_outputs: Dict[str, pd.DataFrame],
    expected_replicates: int = 3,
    min_jaccard: float = 0.5,
    expected_tools: List[str] = None
) -> DataQualityReport:
    """
    Comprehensive validation of all tool outputs.

    Args:
        tool_outputs: Dictionary mapping tool names to DataFrames
        expected_replicates: Expected number of replicates
        min_jaccard: Minimum acceptable Jaccard similarity
        expected_tools: List of expected tool names

    Returns:
        Complete DataQualityReport
    """
    report = DataQualityReport(
        expected_replicates=expected_replicates,
        expected_tools=expected_tools or list(tool_outputs.keys())
    )

    # Check for missing expected tools
    if expected_tools:
        for tool in expected_tools:
            if tool not in tool_outputs or tool_outputs[tool].empty:
                report.add_warning(
                    WarningLevel.ERROR,
                    tool,
                    "Expected tool has no data"
                )

    # Check replicate completeness
    report = check_replicate_completeness(
        tool_outputs, expected_replicates, report
    )

    # Check replicate overlap for each tool
    for tool_name, tool_df in tool_outputs.items():
        if not tool_df.empty:
            report = check_replicate_overlap(tool_df, tool_name, min_jaccard, report)

    # Check position coverage
    report = check_position_coverage(tool_outputs, None, report)

    return report


def generate_quality_summary(report: DataQualityReport) -> pd.DataFrame:
    """
    Generate a summary DataFrame from quality report.

    Args:
        report: DataQualityReport

    Returns:
        Summary DataFrame with one row per tool
    """
    tool_data = {}

    for status in report.tool_statuses:
        tool = status.tool
        if tool not in tool_data:
            tool_data[tool] = {
                'tool': tool,
                'n_replicates_ok': 0,
                'n_replicates_total': 0,
                'total_positions': 0,
                'status': 'OK',
                'warnings': []
            }

        tool_data[tool]['n_replicates_total'] += 1
        if status.status == Status.OK:
            tool_data[tool]['n_replicates_ok'] += 1
        tool_data[tool]['total_positions'] += status.n_positions

        if status.status in [Status.EMPTY, Status.FAILED]:
            tool_data[tool]['status'] = 'INCOMPLETE'

        if status.warnings:
            tool_data[tool]['warnings'].extend(status.warnings)

    # Create DataFrame
    rows = []
    for tool, data in tool_data.items():
        rows.append({
            'tool': data['tool'],
            'replicates_ok': f"{data['n_replicates_ok']}/{data['n_replicates_total']}",
            'total_positions': data['total_positions'],
            'status': data['status'],
            'n_warnings': len(data['warnings'])
        })

    return pd.DataFrame(rows)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Validate tool outputs and generate quality report'
    )
    parser.add_argument('--input-dir', '-i', required=True,
                        help='Directory containing tool outputs')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file prefix')
    parser.add_argument('--expected-replicates', type=int, default=3,
                        help='Expected number of replicates')
    parser.add_argument('--min-jaccard', type=float, default=0.5,
                        help='Minimum acceptable Jaccard similarity')

    args = parser.parse_args()

    # Import parser
    from parse_outputs import load_all_tool_outputs

    # Load outputs
    tool_outputs = load_all_tool_outputs(args.input_dir)

    # Validate
    report = validate_all_outputs(
        tool_outputs,
        expected_replicates=args.expected_replicates,
        min_jaccard=args.min_jaccard
    )

    # Save outputs
    report.save(f"{args.output}_quality_report.csv")
    report.save_warnings(f"{args.output}_warnings.txt")

    # Print summary
    print(report.summary())
