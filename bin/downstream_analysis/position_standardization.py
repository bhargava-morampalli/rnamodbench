#!/usr/bin/env python3
"""
Position Standardization Module for RNA Modification Detection Pipeline

Implements reference-universe standardization per tool/replicate while keeping
raw parsed outputs intact. The evaluation universe can be the full reference or
an explicitly defined sub-interval within padded coordinates.
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
    added_from_universe: int
    full_reference_length: int
    eval_start: int
    eval_end: int
    total_positions: int

    def to_dict(self) -> dict:
        return {
            "tool": self.tool,
            "replicate": self.replicate,
            "reference": self.reference,
            "native_positions": self.native_positions,
            "added_from_universe": self.added_from_universe,
            "full_reference_length": self.full_reference_length,
            "eval_start": self.eval_start,
            "eval_end": self.eval_end,
            "total_positions": self.total_positions,
        }


@dataclass
class StandardizationReport:
    """Complete standardization report."""
    availability: List[PositionAvailability] = field(default_factory=list)
    imputation_records: List[ImputationRecord] = field(default_factory=list)

    def to_availability_dataframe(self) -> pd.DataFrame:
        """Convert availability summaries to DataFrame."""
        if not self.availability:
            return pd.DataFrame(
                columns=[
                    "tool",
                    "replicate",
                    "reference",
                    "native_positions",
                    "added_from_universe",
                    "full_reference_length",
                    "eval_start",
                    "eval_end",
                    "total_positions",
                ]
            )
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


@dataclass
class ReferenceCatalogEntry:
    target: str
    reference_path: str
    reference_id: str
    length: int
    eval_start: int
    eval_end: int
    eval_length: int


def _default_eval_interval(target: str, reference_id: str, reference_length: int) -> Tuple[int, int]:
    key = str(target).strip().lower()
    ref_l = str(reference_id).strip().lower()
    if key == "16s" or ref_l.startswith("16s"):
        return 201, 1742
    if key == "23s" or ref_l.startswith("23s"):
        return 201, 3104
    return 1, int(reference_length)


def _read_fasta_lengths(fasta_path: Path) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    seq_len = 0

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = seq_len
                name = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)

    if name is not None:
        lengths[name] = seq_len

    return lengths


def load_reference_catalog(references_csv: Path) -> pd.DataFrame:
    """
    Load reference mapping and lengths from references CSV.

    Expected columns: target,reference
    Optional columns: eval_start,eval_end
    """
    references_csv = Path(references_csv)
    if not references_csv.exists():
        raise FileNotFoundError(f"references CSV not found: {references_csv}")

    df = pd.read_csv(references_csv)
    if "target" not in df.columns or "reference" not in df.columns:
        raise ValueError("references CSV must have columns: target,reference")

    rows: List[ReferenceCatalogEntry] = []

    for _, row in df.iterrows():
        target = str(row["target"]).strip().lower()
        ref_path = Path(str(row["reference"]).strip())
        if not ref_path.exists():
            raise FileNotFoundError(f"reference FASTA not found: {ref_path}")

        lengths = _read_fasta_lengths(ref_path)
        if not lengths:
            raise ValueError(f"No FASTA records found in {ref_path}")

        # Pipeline references are single-target FASTAs; use first sequence as canonical ID
        ref_id, ref_len = next(iter(lengths.items()))
        default_start, default_end = _default_eval_interval(target, ref_id, int(ref_len))
        eval_start = (
            int(row["eval_start"])
            if "eval_start" in df.columns and pd.notna(row["eval_start"])
            else default_start
        )
        eval_end = (
            int(row["eval_end"])
            if "eval_end" in df.columns and pd.notna(row["eval_end"])
            else default_end
        )
        if eval_start < 1 or eval_end > int(ref_len) or eval_start > eval_end:
            raise ValueError(
                f"Invalid evaluation interval for {target}: "
                f"eval_start={eval_start}, eval_end={eval_end}, reference_length={ref_len}"
            )
        rows.append(
            ReferenceCatalogEntry(
                target=target,
                reference_path=str(ref_path),
                reference_id=ref_id,
                length=int(ref_len),
                eval_start=eval_start,
                eval_end=eval_end,
                eval_length=int(eval_end - eval_start + 1),
            )
        )

    out = pd.DataFrame([r.__dict__ for r in rows])
    return out


def canonicalize_tool_references(
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

        # loose match for prefixes (e.g., 16s_...)
        for target, ref_id in alias_map.items():
            if target in {"16s", "23s", "5s"} and rl.startswith(target):
                return ref_id

        return r

    out["reference"] = out["reference"].astype(str).map(_canon)
    return out


def get_reference_lengths(reference_catalog: pd.DataFrame) -> Dict[str, int]:
    return {str(r["reference_id"]): int(r["length"]) for _, r in reference_catalog.iterrows()}


def get_reference_eval_regions(reference_catalog: pd.DataFrame) -> Dict[str, Tuple[int, int]]:
    return {
        str(r["reference_id"]): (int(r["eval_start"]), int(r["eval_end"]))
        for _, r in reference_catalog.iterrows()
    }


def standardize_to_reference_universe(
    tool_df: pd.DataFrame,
    reference_id: str,
    reference_length: int,
    universe_positions: Optional[List[int]] = None,
    eval_start: int = 1,
    eval_end: Optional[int] = None,
    replicates: Optional[List[str]] = None,
    report: Optional[StandardizationReport] = None,
) -> pd.DataFrame:
    """
    Standardize one tool to the evaluation universe for one reference.

    Output includes:
    - score_raw
    - pvalue_raw
    - is_reported
    - is_imputed
    - imputation_source
    """
    if universe_positions is None:
        eval_end = int(reference_length) if eval_end is None else int(eval_end)
        eval_start = int(eval_start)
        if eval_start < 1 or eval_end > int(reference_length) or eval_start > eval_end:
            raise ValueError(
                f"Invalid evaluation interval for {reference_id}: "
                f"eval_start={eval_start}, eval_end={eval_end}, reference_length={reference_length}"
            )
        universe = list(range(eval_start, eval_end + 1))
    else:
        universe = sorted({int(pos) for pos in universe_positions})
        if not universe:
            raise ValueError(f"Empty evaluation universe for {reference_id}")
        if universe[0] < 1 or universe[-1] > int(reference_length):
            raise ValueError(
                f"Invalid evaluation universe for {reference_id}: "
                f"min={universe[0]}, max={universe[-1]}, reference_length={reference_length}"
            )
        eval_start = universe[0]
        eval_end = universe[-1]

    if replicates is None:
        if tool_df.empty or "replicate" not in tool_df.columns:
            replicates = ["rep1"]
        else:
            replicates = sorted(tool_df["replicate"].astype(str).unique().tolist())

    tool_name = (
        str(tool_df["tool"].iloc[0])
        if (not tool_df.empty and "tool" in tool_df.columns)
        else "unknown"
    )

    df_ref = tool_df.copy()
    if not df_ref.empty:
        df_ref = df_ref[df_ref["reference"].astype(str) == str(reference_id)].copy()
        df_ref["position"] = pd.to_numeric(df_ref["position"], errors="coerce")
        df_ref = df_ref[df_ref["position"].notna()].copy()
        df_ref["position"] = df_ref["position"].astype(int)

    all_rows: List[dict] = []

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

        for pos in universe:
            if pos in reported_by_pos:
                row = reported_by_pos[pos]
                all_rows.append(
                    {
                        "tool": tool_name,
                        "reference": reference_id,
                        "position": pos,
                        "replicate": str(rep),
                        "sample_id": row.get("sample_id", f"{tool_name}_{rep}"),
                        "score_type": row.get("score_type", "unknown"),
                        "score_raw": pd.to_numeric(row.get("score", np.nan), errors="coerce"),
                        "pvalue_raw": pd.to_numeric(row.get("pvalue", np.nan), errors="coerce"),
                        "coverage": row.get("coverage", None),
                        "is_reported": True,
                        "is_imputed": False,
                        "imputation_source": "reported",
                    }
                )
            else:
                added_from_universe += 1
                all_rows.append(
                    {
                        "tool": tool_name,
                        "reference": reference_id,
                        "position": pos,
                        "replicate": str(rep),
                        "sample_id": f"{tool_name}_{rep}",
                        "score_type": rep_df["score_type"].iloc[0] if not rep_df.empty else "unknown",
                        "score_raw": np.nan,
                        "pvalue_raw": np.nan,
                        "coverage": rep_df["coverage"].iloc[0] if (not rep_df.empty and "coverage" in rep_df.columns) else None,
                        "is_reported": False,
                        "is_imputed": True,
                        "imputation_source": "reference_universe",
                    }
                )

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
                    replicate=str(rep),
                    reference=reference_id,
                    native_positions=native_positions,
                    added_from_universe=added_from_universe,
                    full_reference_length=int(reference_length),
                    eval_start=eval_start,
                    eval_end=eval_end,
                    total_positions=len(universe),
                )
            )

    missing_df = pd.DataFrame(missing_rows)
    return pd.concat([tool_df, missing_df], ignore_index=True)


def count_within_file_nan(tool_df: pd.DataFrame) -> int:
    """Count positions with NaN score that are not from imputation."""
    if tool_df.empty:
        return 0

    st = str(score_type).lower()
    value = float(score_raw)

    if st in {"pvalue", "fdr"}:
        value = max(min(value, 1.0), 1e-300)
        return float(-np.log10(value))

    if st == "neglog10_fdr":
        return float(value)

    if st == "zscore":
        return float(abs(value))

    return float(value)


def _sanitize_score_metric(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    arr = values.to_numpy(copy=True)

    pos_inf = np.isposinf(arr)
    neg_inf = np.isneginf(arr)
    if not (pos_inf.any() or neg_inf.any()):
        return values

    finite = arr[np.isfinite(arr)]
    if finite.size:
        upper = float(finite.max() + 1.0)
        lower = float(finite.min() - 1.0)
    else:
        upper = 1000.0
        lower = -1000.0

    arr[pos_inf] = upper
    arr[neg_inf] = lower
    return pd.Series(arr, index=series.index, dtype=float)


def build_metric_ready_table(
    imputed_df: pd.DataFrame,
    ground_truth_positions: Optional[Set[int]] = None,
) -> pd.DataFrame:
    """
    Add score_metric + label columns to imputed table.

    ground_truth_positions must be position integers already filtered to one reference.
    """
    if imputed_df.empty:
        return imputed_df

    out = imputed_df.copy()
    out["score_metric"] = [
        _score_metric_from_raw(st, sc)
        for st, sc in zip(out["score_type"].values, out["score_raw"].values)
    ]
    out["score_metric"] = _sanitize_score_metric(out["score_metric"])

    if ground_truth_positions is None:
        out["label"] = 0
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
