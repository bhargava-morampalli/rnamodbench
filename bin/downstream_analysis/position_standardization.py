#!/usr/bin/env python3
"""
Position standardization utilities.

Implements full-reference universe standardization (1..L) per tool/replicate,
while keeping raw parsed outputs intact.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


@dataclass
class ImputationRecord:
    tool: str
    replicate: str
    reference: str
    position: int
    source: str


@dataclass
class PositionAvailability:
    tool: str
    replicate: str
    reference: str
    native_positions: int
    added_from_universe: int
    total_positions: int

    def to_dict(self) -> dict:
        return {
            "tool": self.tool,
            "replicate": self.replicate,
            "reference": self.reference,
            "native_positions": self.native_positions,
            "added_from_universe": self.added_from_universe,
            "total_positions": self.total_positions,
        }


@dataclass
class StandardizationReport:
    availability: List[PositionAvailability] = field(default_factory=list)
    imputation_records: List[ImputationRecord] = field(default_factory=list)

    def to_availability_dataframe(self) -> pd.DataFrame:
        if not self.availability:
            return pd.DataFrame(
                columns=[
                    "tool",
                    "replicate",
                    "reference",
                    "native_positions",
                    "added_from_universe",
                    "total_positions",
                ]
            )
        return pd.DataFrame([a.to_dict() for a in self.availability])

    def to_imputation_dataframe(self) -> pd.DataFrame:
        if not self.imputation_records:
            return pd.DataFrame(columns=["tool", "replicate", "reference", "position", "source"])
        return pd.DataFrame(
            [
                {
                    "tool": rec.tool,
                    "replicate": rec.replicate,
                    "reference": rec.reference,
                    "position": rec.position,
                    "source": rec.source,
                }
                for rec in self.imputation_records
            ]
        )

    def save(self, output_dir: Path) -> None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        self.to_availability_dataframe().to_csv(
            output_dir / "position_availability_report.csv", index=False
        )
        self.to_imputation_dataframe().to_csv(
            output_dir / "imputed_positions_detail.csv", index=False
        )

    def summary(self) -> str:
        if not self.availability:
            return "No standardization records"
        df = self.to_availability_dataframe()
        total_native = int(df["native_positions"].sum())
        total_added = int(df["added_from_universe"].sum())
        total_rows = int(df["total_positions"].sum())
        return (
            "POSITION STANDARDIZATION SUMMARY\n"
            f"native_positions={total_native}\n"
            f"added_from_universe={total_added}\n"
            f"total_positions={total_rows}"
        )


@dataclass
class ReferenceCatalogEntry:
    target: str
    reference_path: str
    reference_id: str
    length: int


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
        rows.append(
            ReferenceCatalogEntry(
                target=target,
                reference_path=str(ref_path),
                reference_id=ref_id,
                length=int(ref_len),
            )
        )

    out = pd.DataFrame([r.__dict__ for r in rows])
    return out


def canonicalize_tool_references(
    tool_df: pd.DataFrame,
    reference_catalog: pd.DataFrame,
) -> pd.DataFrame:
    """
    Canonicalize references to FASTA header IDs.

    Maps common aliases (e.g. 16s -> 16s_88_rrsE) while preserving already
    canonical IDs.
    """
    if tool_df.empty:
        return tool_df

    ref_df = reference_catalog.copy()
    alias_map = {}
    for _, row in ref_df.iterrows():
        target = str(row["target"]).lower()
        ref_id = str(row["reference_id"])
        alias_map[target] = ref_id
        alias_map[ref_id.lower()] = ref_id

    out = tool_df.copy()

    def _canon(ref: str) -> str:
        r = str(ref).strip()
        rl = r.lower()
        if rl in alias_map:
            return alias_map[rl]

        # loose match for prefixes (e.g., 16s_...)
        for target, ref_id in alias_map.items():
            if target in {"16s", "23s", "5s"} and rl.startswith(target):
                return ref_id

        return r

    out["reference"] = out["reference"].astype(str).map(_canon)
    return out


def get_reference_lengths(reference_catalog: pd.DataFrame) -> Dict[str, int]:
    return {str(r["reference_id"]): int(r["length"]) for _, r in reference_catalog.iterrows()}


def standardize_to_reference_universe(
    tool_df: pd.DataFrame,
    reference_id: str,
    reference_length: int,
    replicates: Optional[List[str]] = None,
    report: Optional[StandardizationReport] = None,
) -> pd.DataFrame:
    """
    Standardize one tool to full position universe (1..reference_length) for one reference.

    Output includes:
    - score_raw
    - pvalue_raw
    - is_reported
    - is_imputed
    - imputation_source
    """
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

    for rep in replicates:
        rep_df = df_ref[df_ref["replicate"].astype(str) == str(rep)].copy() if not df_ref.empty else pd.DataFrame()

        # One row per position for reported calls (first row kept on duplicates)
        reported_by_pos = {}
        if not rep_df.empty:
            rep_df = rep_df.sort_values("position")
            for _, row in rep_df.iterrows():
                pos = int(row["position"])
                if pos not in reported_by_pos:
                    reported_by_pos[pos] = row

        native_positions = len(reported_by_pos)
        added_from_universe = 0

        for pos in range(1, int(reference_length) + 1):
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

                if report is not None:
                    report.imputation_records.append(
                        ImputationRecord(
                            tool=tool_name,
                            replicate=str(rep),
                            reference=reference_id,
                            position=pos,
                            source="reference_universe",
                        )
                    )

        if report is not None:
            report.availability.append(
                PositionAvailability(
                    tool=tool_name,
                    replicate=str(rep),
                    reference=reference_id,
                    native_positions=native_positions,
                    added_from_universe=added_from_universe,
                    total_positions=int(reference_length),
                )
            )

    out = pd.DataFrame(all_rows)
    return out


def _score_metric_from_raw(score_type: str, score_raw: float) -> float:
    if pd.isna(score_raw):
        return 0.0

    st = str(score_type).lower()
    value = float(score_raw)

    if st in {"pvalue", "fdr"}:
        value = max(min(value, 1.0), 1e-300)
        return float(-np.log10(value))

    if st == "zscore":
        return float(abs(value))

    return float(value)


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

    if ground_truth_positions is None:
        out["label"] = 0
    else:
        gt = {int(p) for p in ground_truth_positions}
        out["label"] = out["position"].astype(int).isin(gt).astype(int)

    return out


def compute_no_call_diagnostics(
    metric_ready_df: pd.DataFrame,
    ground_truth_positions: Optional[Set[int]] = None,
) -> pd.DataFrame:
    """Compute no-call diagnostics per tool/replicate/reference."""
    if metric_ready_df.empty:
        return pd.DataFrame(
            columns=[
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
        )

    gt = {int(p) for p in (ground_truth_positions or set())}

    rows = []
    grp_cols = ["reference", "tool", "replicate"]
    for (ref, tool, rep), sub in metric_ready_df.groupby(grp_cols, dropna=False):
        n_universe = len(sub)
        n_reported = int(sub["is_reported"].sum())
        n_no_call = int(n_universe - n_reported)
        no_call_rate = float(n_no_call / n_universe) if n_universe else np.nan

        n_gt_total = len(gt)
        if n_gt_total > 0:
            gt_rows = sub[sub["position"].astype(int).isin(gt)]
            n_gt_reported = int(gt_rows["is_reported"].sum())
            gt_recall_raw = float(n_gt_reported / n_gt_total)
        else:
            n_gt_reported = 0
            gt_recall_raw = np.nan

        rows.append(
            {
                "reference": ref,
                "tool": tool,
                "replicate": rep,
                "n_universe": n_universe,
                "n_reported": n_reported,
                "n_no_call": n_no_call,
                "no_call_rate": no_call_rate,
                "n_gt_total": n_gt_total,
                "n_gt_reported": n_gt_reported,
                "gt_recall_raw": gt_recall_raw,
            }
        )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Backwards-compatible wrappers used by earlier code paths
# -----------------------------------------------------------------------------

def standardize_all_tool_outputs(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: Optional[pd.DataFrame] = None,
):
    """Legacy wrapper retained for compatibility with older run_analysis."""
    del ground_truth
    report = StandardizationReport()
    return tool_outputs, report


def generate_availability_report(
    tool_outputs: Dict[str, pd.DataFrame],
    output_dir: Path,
) -> StandardizationReport:
    report = StandardizationReport()
    for tool, df in tool_outputs.items():
        if df.empty:
            continue
        reps = sorted(df["replicate"].astype(str).unique()) if "replicate" in df.columns else ["rep1"]
        refs = sorted(df["reference"].astype(str).unique()) if "reference" in df.columns else ["unknown"]
        for rep in reps:
            for ref in refs:
                sub = df[(df["replicate"].astype(str) == rep) & (df["reference"].astype(str) == ref)]
                report.availability.append(
                    PositionAvailability(
                        tool=tool,
                        replicate=rep,
                        reference=ref,
                        native_positions=len(sub),
                        added_from_universe=0,
                        total_positions=len(sub),
                    )
                )
    report.save(Path(output_dir))
    return report


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Position standardization helpers")
    parser.add_argument("--references-csv", required=True)
    args = parser.parse_args()

    catalog = load_reference_catalog(Path(args.references_csv))
    print(catalog.to_string(index=False))
