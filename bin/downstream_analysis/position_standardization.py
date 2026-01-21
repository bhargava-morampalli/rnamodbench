#!/usr/bin/env python3
"""
Position Standardization Module for RNA Modification Detection Pipeline

Ensures all tools report the same universe of positions, with NaN for missing
scores. This enables proper true negative calculation and fair comparisons.

Features:
- Within-file NaN handling (blank scores → NaN)
- Across-replicate position standardization (union of positions)
- Ground-truth position filling (add GT positions not in tool output)
- Position availability reporting (track what was imputed)

Usage:
    from position_standardization import (
        standardize_replicate_positions,
        fill_ground_truth_positions,
        standardize_all_tool_outputs,
        generate_availability_report,
    )

    # Standardize positions across replicates
    tool_df = standardize_replicate_positions(tool_df)

    # Fill missing ground truth positions
    tool_df = fill_ground_truth_positions(tool_df, ground_truth)
"""

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
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
class ImputationRecord:
    """Record of a single imputed position."""
    tool: str
    replicate: str
    reference: str
    position: int
    source: str  # 'within_file', 'other_rep', 'ground_truth'


@dataclass
class PositionAvailability:
    """Summary of position availability for a tool/replicate."""
    tool: str
    replicate: str
    native_positions: int
    within_file_nan: int
    added_from_reps: int
    added_from_gt: int
    total_positions: int

    def to_dict(self) -> dict:
        return {
            'tool': self.tool,
            'replicate': self.replicate,
            'native_positions': self.native_positions,
            'within_file_nan': self.within_file_nan,
            'added_from_reps': self.added_from_reps,
            'added_from_gt': self.added_from_gt,
            'total_positions': self.total_positions,
        }


@dataclass
class StandardizationReport:
    """Complete standardization report."""
    availability: List[PositionAvailability] = field(default_factory=list)
    imputation_records: List[ImputationRecord] = field(default_factory=list)

    def to_availability_dataframe(self) -> pd.DataFrame:
        """Convert availability summaries to DataFrame."""
        if not self.availability:
            return pd.DataFrame(columns=[
                'tool', 'replicate', 'native_positions', 'within_file_nan',
                'added_from_reps', 'added_from_gt', 'total_positions'
            ])
        return pd.DataFrame([a.to_dict() for a in self.availability])

    def to_imputation_dataframe(self) -> pd.DataFrame:
        """Convert imputation records to DataFrame."""
        if not self.imputation_records:
            return pd.DataFrame(columns=[
                'tool', 'replicate', 'reference', 'position', 'source'
            ])
        return pd.DataFrame([
            {
                'tool': r.tool,
                'replicate': r.replicate,
                'reference': r.reference,
                'position': r.position,
                'source': r.source,
            }
            for r in self.imputation_records
        ])

    def save(self, output_dir: Path) -> None:
        """Save reports to output directory."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save availability summary
        self.to_availability_dataframe().to_csv(
            output_dir / 'position_availability_report.csv', index=False
        )

        # Save detailed imputation records
        self.to_imputation_dataframe().to_csv(
            output_dir / 'imputed_positions_detail.csv', index=False
        )

        logger.info(f"Saved position standardization reports to {output_dir}")

    def summary(self) -> str:
        """Generate text summary."""
        lines = [
            "=" * 60,
            "POSITION STANDARDIZATION SUMMARY",
            "=" * 60,
            "",
        ]

        if not self.availability:
            lines.append("No standardization performed.")
            return '\n'.join(lines)

        # Total statistics
        total_native = sum(a.native_positions for a in self.availability)
        total_within_nan = sum(a.within_file_nan for a in self.availability)
        total_from_reps = sum(a.added_from_reps for a in self.availability)
        total_from_gt = sum(a.added_from_gt for a in self.availability)

        lines.extend([
            f"Total native positions: {total_native}",
            f"Within-file NaN values: {total_within_nan}",
            f"Positions added from other replicates: {total_from_reps}",
            f"Positions added from ground truth: {total_from_gt}",
            "",
            "Per-tool breakdown:",
            "-" * 40,
        ])

        # Group by tool
        tools = sorted(set(a.tool for a in self.availability))
        for tool in tools:
            tool_avail = [a for a in self.availability if a.tool == tool]
            native = sum(a.native_positions for a in tool_avail)
            from_reps = sum(a.added_from_reps for a in tool_avail)
            from_gt = sum(a.added_from_gt for a in tool_avail)
            lines.append(f"  {tool}: {native} native, +{from_reps} from reps, +{from_gt} from GT")

        lines.append("")
        lines.append("=" * 60)
        return '\n'.join(lines)


def standardize_replicate_positions(
    tool_df: pd.DataFrame,
    report: Optional[StandardizationReport] = None
) -> pd.DataFrame:
    """
    Ensure all replicates have the same positions (union).
    Missing positions filled with NaN score.

    Args:
        tool_df: DataFrame with tool outputs (must have replicate column)
        report: Optional report to record imputation details

    Returns:
        Standardized DataFrame with all positions across all replicates
    """
    if tool_df.empty or 'replicate' not in tool_df.columns:
        return tool_df

    replicates = tool_df['replicate'].unique()
    if len(replicates) <= 1:
        # Only one replicate, nothing to standardize
        return tool_df

    # Add _imputed column if not present
    if '_imputed' not in tool_df.columns:
        tool_df = tool_df.copy()
        tool_df['_imputed'] = False

    tool_name = tool_df['tool'].iloc[0] if 'tool' in tool_df.columns else 'unknown'
    logger.info(f"Standardizing positions across {len(replicates)} replicates for {tool_name}")

    # Build position -> reference mapping from all replicates
    position_ref_map = {}
    for _, row in tool_df[['reference', 'position']].drop_duplicates().iterrows():
        position_ref_map[row['position']] = row['reference']

    # Get positions per replicate
    rep_positions = {}
    for rep in replicates:
        rep_df = tool_df[tool_df['replicate'] == rep]
        rep_positions[rep] = set(rep_df['position'].unique())

    # Calculate union of all positions
    all_positions = set.union(*rep_positions.values())
    logger.info(f"Union of positions: {len(all_positions)}")

    # Standardize each replicate
    standardized_dfs = []
    for rep in replicates:
        rep_df = tool_df[tool_df['replicate'] == rep].copy()
        current_positions = rep_positions[rep]
        missing_positions = all_positions - current_positions

        if missing_positions:
            logger.debug(f"Adding {len(missing_positions)} positions to replicate {rep}")

            # Create rows for missing positions
            missing_rows = []
            for pos in missing_positions:
                ref = position_ref_map.get(pos, 'unknown')
                missing_rows.append({
                    'tool': tool_name,
                    'reference': ref,
                    'position': pos,
                    'score': np.nan,
                    'score_type': rep_df['score_type'].iloc[0] if not rep_df.empty else 'unknown',
                    'pvalue': np.nan,
                    'replicate': rep,
                    'sample_id': rep_df['sample_id'].iloc[0] if not rep_df.empty else '',
                    'coverage': rep_df['coverage'].iloc[0] if 'coverage' in rep_df.columns and not rep_df.empty else None,
                    '_imputed': True,
                    '_imputation_source': 'other_rep',
                })

                # Record imputation
                if report is not None:
                    report.imputation_records.append(ImputationRecord(
                        tool=tool_name,
                        replicate=rep,
                        reference=ref,
                        position=pos,
                        source='other_rep',
                    ))

            missing_df = pd.DataFrame(missing_rows)
            rep_df = pd.concat([rep_df, missing_df], ignore_index=True)

        standardized_dfs.append(rep_df)

    return pd.concat(standardized_dfs, ignore_index=True)


def fill_ground_truth_positions(
    tool_df: pd.DataFrame,
    ground_truth: pd.DataFrame,
    report: Optional[StandardizationReport] = None
) -> pd.DataFrame:
    """
    Add ground truth positions that tool didn't report (with NaN score).

    This is CRITICAL for accurate metrics because:
    - If tool reports position 500 with high score → can be True Positive
    - If tool reports position 500 with low score → can be True Negative
    - If tool DOESN'T report position 500 → we don't know!

    By filling with NaN, we mark these as "not reported" and downstream
    analysis can handle them appropriately (exclude from curve calculation).

    Args:
        tool_df: DataFrame with tool outputs
        ground_truth: Ground truth DataFrame with reference and position columns
        report: Optional report to record imputation details

    Returns:
        DataFrame with ground truth positions added (NaN score where missing)
    """
    if tool_df.empty or ground_truth.empty:
        return tool_df

    # Add _imputed column if not present
    if '_imputed' not in tool_df.columns:
        tool_df = tool_df.copy()
        tool_df['_imputed'] = False

    tool_name = tool_df['tool'].iloc[0] if 'tool' in tool_df.columns else 'unknown'

    # Get ground truth positions
    gt_positions = set(zip(ground_truth['reference'], ground_truth['position']))

    # Get tool positions
    tool_positions = set(zip(tool_df['reference'], tool_df['position']))

    # Find missing ground truth positions
    missing_gt_positions = gt_positions - tool_positions

    if not missing_gt_positions:
        logger.info(f"{tool_name}: All ground truth positions present")
        return tool_df

    logger.info(f"{tool_name}: Adding {len(missing_gt_positions)} ground truth positions")

    # Get replicates
    replicates = tool_df['replicate'].unique() if 'replicate' in tool_df.columns else ['rep1']

    # Create rows for missing positions (for each replicate)
    missing_rows = []
    for ref, pos in missing_gt_positions:
        for rep in replicates:
            rep_df = tool_df[tool_df['replicate'] == rep]

            missing_rows.append({
                'tool': tool_name,
                'reference': ref,
                'position': pos,
                'score': np.nan,
                'score_type': rep_df['score_type'].iloc[0] if not rep_df.empty else 'unknown',
                'pvalue': np.nan,
                'replicate': rep,
                'sample_id': rep_df['sample_id'].iloc[0] if not rep_df.empty else '',
                'coverage': rep_df['coverage'].iloc[0] if 'coverage' in rep_df.columns and not rep_df.empty else None,
                '_imputed': True,
                '_imputation_source': 'ground_truth',
            })

            # Record imputation
            if report is not None:
                report.imputation_records.append(ImputationRecord(
                    tool=tool_name,
                    replicate=rep,
                    reference=ref,
                    position=pos,
                    source='ground_truth',
                ))

    missing_df = pd.DataFrame(missing_rows)
    return pd.concat([tool_df, missing_df], ignore_index=True)


def count_within_file_nan(tool_df: pd.DataFrame) -> int:
    """Count positions with NaN score that are not from imputation."""
    if tool_df.empty:
        return 0

    # NaN scores that are not imputed = within-file NaN
    if '_imputed' in tool_df.columns:
        native_rows = tool_df[~tool_df['_imputed']]
    else:
        native_rows = tool_df

    return native_rows['score'].isna().sum()


def standardize_all_tool_outputs(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: Optional[pd.DataFrame] = None,
) -> Tuple[Dict[str, pd.DataFrame], StandardizationReport]:
    """
    Standardize positions across all tools.

    This function:
    1. Standardizes positions across replicates (within each tool)
    2. Fills missing ground truth positions (if ground truth provided)
    3. Generates availability report

    Args:
        tool_outputs: Dictionary mapping tool names to DataFrames
        ground_truth: Optional ground truth DataFrame

    Returns:
        Tuple of (standardized tool outputs, standardization report)
    """
    report = StandardizationReport()
    standardized_outputs = {}

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            standardized_outputs[tool_name] = tool_df
            continue

        logger.info(f"Standardizing {tool_name}...")

        # Count within-file NaN before standardization
        within_file_nan = count_within_file_nan(tool_df)

        # Count native positions per replicate BEFORE standardization
        native_counts = {}
        if 'replicate' in tool_df.columns:
            for rep in tool_df['replicate'].unique():
                rep_df = tool_df[tool_df['replicate'] == rep]
                native_counts[rep] = len(rep_df)

        # Step 1: Standardize across replicates
        standardized_df = standardize_replicate_positions(tool_df, report)

        # Count positions added from other replicates
        added_from_reps = {}
        if '_imputation_source' in standardized_df.columns:
            for rep in standardized_df['replicate'].unique():
                rep_df = standardized_df[standardized_df['replicate'] == rep]
                from_reps = (rep_df['_imputation_source'] == 'other_rep').sum()
                added_from_reps[rep] = from_reps

        # Step 2: Fill ground truth positions (if provided)
        pre_gt_len = len(standardized_df)
        if ground_truth is not None and not ground_truth.empty:
            standardized_df = fill_ground_truth_positions(standardized_df, ground_truth, report)

        # Count positions added from ground truth
        added_from_gt = {}
        if '_imputation_source' in standardized_df.columns:
            for rep in standardized_df['replicate'].unique():
                rep_df = standardized_df[standardized_df['replicate'] == rep]
                from_gt = (rep_df['_imputation_source'] == 'ground_truth').sum()
                added_from_gt[rep] = from_gt

        # Record availability for each replicate
        if 'replicate' in standardized_df.columns:
            for rep in standardized_df['replicate'].unique():
                rep_df = standardized_df[standardized_df['replicate'] == rep]
                report.availability.append(PositionAvailability(
                    tool=tool_name,
                    replicate=rep,
                    native_positions=native_counts.get(rep, len(rep_df)),
                    within_file_nan=within_file_nan // len(native_counts) if native_counts else within_file_nan,
                    added_from_reps=added_from_reps.get(rep, 0),
                    added_from_gt=added_from_gt.get(rep, 0),
                    total_positions=len(rep_df),
                ))

        standardized_outputs[tool_name] = standardized_df

    logger.info(f"Standardization complete. {len(report.imputation_records)} positions imputed.")
    return standardized_outputs, report


def generate_availability_report(
    tool_outputs: Dict[str, pd.DataFrame],
    output_dir: Path
) -> StandardizationReport:
    """
    Generate position availability report without modifying data.

    This analyzes the current state of tool outputs and reports
    on position coverage and NaN values.

    Args:
        tool_outputs: Dictionary mapping tool names to DataFrames
        output_dir: Directory to save report

    Returns:
        StandardizationReport with availability information
    """
    report = StandardizationReport()

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            continue

        if 'replicate' not in tool_df.columns:
            # Single replicate case
            within_nan = tool_df['score'].isna().sum()
            imputed = tool_df['_imputed'].sum() if '_imputed' in tool_df.columns else 0

            report.availability.append(PositionAvailability(
                tool=tool_name,
                replicate='all',
                native_positions=len(tool_df) - imputed,
                within_file_nan=within_nan - imputed,  # NaN that weren't imputed
                added_from_reps=0,
                added_from_gt=imputed,  # Assuming imputed = from GT
                total_positions=len(tool_df),
            ))
        else:
            for rep in tool_df['replicate'].unique():
                rep_df = tool_df[tool_df['replicate'] == rep]
                within_nan = rep_df['score'].isna().sum()
                imputed = rep_df['_imputed'].sum() if '_imputed' in rep_df.columns else 0

                report.availability.append(PositionAvailability(
                    tool=tool_name,
                    replicate=rep,
                    native_positions=len(rep_df) - imputed,
                    within_file_nan=max(0, within_nan - imputed),
                    added_from_reps=0,  # Can't distinguish without full standardization
                    added_from_gt=imputed,
                    total_positions=len(rep_df),
                ))

    report.save(output_dir)
    return report


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Standardize positions across tool outputs'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input directory with tool outputs')
    parser.add_argument('--ground-truth', '-g',
                        help='Ground truth file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output directory')

    args = parser.parse_args()

    from parse_outputs import load_all_tool_outputs
    from benchmark_metrics import load_ground_truth

    # Load data
    tool_outputs = load_all_tool_outputs(args.input)

    ground_truth = None
    if args.ground_truth:
        ground_truth = load_ground_truth(args.ground_truth)

    # Standardize
    standardized, report = standardize_all_tool_outputs(tool_outputs, ground_truth)

    # Save report
    output_dir = Path(args.output)
    report.save(output_dir)

    # Print summary
    print(report.summary())
