#!/usr/bin/env python3
"""
Tool comparison for RNA modification detection outputs.

All overlap calculations are reference-aware using (reference, position) tuples.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

Key = Tuple[str, int]


@dataclass
class ToolComparisonResult:
    n_tools: int
    tool_names: List[str]
    n_positions_per_tool: Dict[str, int]
    n_union: int
    n_intersection: int
    pairwise_overlap: Dict[Tuple[str, str], int]
    pairwise_jaccard: Dict[Tuple[str, str], float]
    positions_by_n_tools: Dict[int, int]
    consensus_positions: Set[Key]

    def to_dict(self) -> Dict:
        out = {
            "n_tools": self.n_tools,
            "tool_names": ",".join(self.tool_names),
            "n_union": self.n_union,
            "n_intersection": self.n_intersection,
        }
        for tool, n in self.n_positions_per_tool.items():
            out[f"n_{tool}"] = n
        for n_tools, count in self.positions_by_n_tools.items():
            out[f"positions_by_{n_tools}_tools"] = count
        return out


def _filter_reference(df: pd.DataFrame, reference: Optional[str]) -> pd.DataFrame:
    if reference is None or df.empty or "reference" not in df.columns:
        return df
    ref_l = reference.strip().lower()
    return df[df["reference"].astype(str).str.lower() == ref_l].copy()


def get_tool_positions(
    tool_output: pd.DataFrame,
    threshold: float = None,
    reference: str = None,
    include_imputed: bool = False,
) -> Set[Key]:
    df = _filter_reference(tool_output.copy(), reference)
    if df.empty:
        return set()

    if not include_imputed and "_imputed" in df.columns:
        df = df[~df["_imputed"]]

    df = df[df["reference"].notna()]
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df = df[df["position"].notna()]

    if threshold is not None and "score" in df.columns:
        score_type = str(df["score_type"].iloc[0]).lower() if "score_type" in df.columns else "pvalue"
        minimize_score_types = {"pvalue", "fdr"}
        score_vals = pd.to_numeric(df["score"], errors="coerce")
        if score_type in minimize_score_types:
            df = df[score_vals <= threshold]
        else:
            df = df[score_vals >= threshold]

    return set(zip(df["reference"].astype(str), df["position"].astype(int)))


def compare_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    reference: str = None,
    min_tools: int = 2,
) -> ToolComparisonResult:
    logger.info("Comparing %d tools", len(tool_outputs))

    tool_positions: Dict[str, Set[Key]] = {}
    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            continue
        positions = get_tool_positions(tool_df, threshold, reference)
        if positions:
            tool_positions[tool_name] = positions

    if not tool_positions:
        return ToolComparisonResult(
            n_tools=0,
            tool_names=[],
            n_positions_per_tool={},
            n_union=0,
            n_intersection=0,
            pairwise_overlap={},
            pairwise_jaccard={},
            positions_by_n_tools={},
            consensus_positions=set(),
        )

    tool_names = sorted(tool_positions.keys())
    n_tools = len(tool_names)

    all_positions: Set[Key] = set().union(*tool_positions.values())
    common_positions = set.intersection(*tool_positions.values()) if tool_positions else set()

    pairwise_overlap = {}
    pairwise_jaccard = {}
    for t1, t2 in combinations(tool_names, 2):
        p1, p2 = tool_positions[t1], tool_positions[t2]
        overlap = len(p1 & p2)
        union = len(p1 | p2)
        pairwise_overlap[(t1, t2)] = overlap
        pairwise_jaccard[(t1, t2)] = overlap / union if union else 0.0

    positions_by_n_tools = defaultdict(int)
    consensus_positions: Set[Key] = set()

    for key in all_positions:
        n_detecting = sum(1 for pos_set in tool_positions.values() if key in pos_set)
        positions_by_n_tools[n_detecting] += 1
        if n_detecting >= min_tools:
            consensus_positions.add(key)

    return ToolComparisonResult(
        n_tools=n_tools,
        tool_names=tool_names,
        n_positions_per_tool={k: len(v) for k, v in tool_positions.items()},
        n_union=len(all_positions),
        n_intersection=len(common_positions),
        pairwise_overlap=pairwise_overlap,
        pairwise_jaccard=pairwise_jaccard,
        positions_by_n_tools=dict(positions_by_n_tools),
        consensus_positions=consensus_positions,
    )


def generate_upset_data(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    reference: str = None,
) -> pd.DataFrame:
    tool_positions = {
        tool: get_tool_positions(df, threshold, reference)
        for tool, df in tool_outputs.items()
        if not df.empty
    }

    if not tool_positions:
        return pd.DataFrame()

    all_positions: Set[Key] = set().union(*tool_positions.values())
    tool_names = sorted(tool_positions.keys())

    rows = []
    for ref, pos in sorted(all_positions):
        row = {"reference": ref, "position": pos}
        for tool in tool_names:
            row[tool] = int((ref, pos) in tool_positions[tool])
        rows.append(row)

    return pd.DataFrame(rows)


def generate_venn_data(
    tool_outputs: Dict[str, pd.DataFrame],
    tools: List[str] = None,
    threshold: float = None,
    reference: str = None,
) -> Dict[str, int]:
    available = [t for t, df in tool_outputs.items() if not df.empty]

    if tools is None:
        tools = available[:3]
    else:
        tools = [t for t in tools if t in available]

    if len(tools) < 2 or len(tools) > 4:
        raise ValueError(f"Venn diagram requires 2-4 tools, got {len(tools)}")

    sets = {t: get_tool_positions(tool_outputs[t], threshold, reference) for t in tools}

    if len(tools) == 2:
        a, b = tools
        return {
            "10": len(sets[a] - sets[b]),
            "01": len(sets[b] - sets[a]),
            "11": len(sets[a] & sets[b]),
            "labels": (a, b),
        }

    if len(tools) == 3:
        a, b, c = tools
        pa, pb, pc = sets[a], sets[b], sets[c]
        return {
            "100": len(pa - pb - pc),
            "010": len(pb - pa - pc),
            "001": len(pc - pa - pb),
            "110": len((pa & pb) - pc),
            "101": len((pa & pc) - pb),
            "011": len((pb & pc) - pa),
            "111": len(pa & pb & pc),
            "labels": (a, b, c),
        }

    a, b, c, d = tools
    pa, pb, pc, pd_ = sets[a], sets[b], sets[c], sets[d]
    return {
        "sets": {t: len(sets[t]) for t in tools},
        "intersection_all": len(pa & pb & pc & pd_),
        "labels": tuple(tools),
    }


def _import_calculate_metrics():
    try:
        from .benchmark_metrics import calculate_metrics
    except Exception:
        from benchmark_metrics import calculate_metrics
    return calculate_metrics


def rank_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    metric: str = "auprc",
    reference: str = None,
) -> pd.DataFrame:
    calculate_metrics = _import_calculate_metrics()
    rows = []

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            continue
        try:
            metrics = calculate_metrics(tool_df, ground_truth, reference)
            rows.append(
                {
                    "tool": tool_name,
                    "auprc": metrics.auprc,
                    "auroc": metrics.auroc,
                    "f1": metrics.f1_optimal,
                    "precision": metrics.precision_optimal,
                    "recall": metrics.recall_optimal,
                    "n_predictions": metrics.n_predictions,
                }
            )
        except Exception as exc:
            logger.warning("Could not rank %s: %s", tool_name, exc)

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out = out.sort_values(metric, ascending=False).reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def calculate_pairwise_agreement(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    reference: str = None,
) -> pd.DataFrame:
    tool_positions = {
        t: get_tool_positions(df, threshold, reference)
        for t, df in tool_outputs.items()
        if not df.empty
    }
    tool_positions = {k: v for k, v in tool_positions.items() if v}

    if len(tool_positions) < 2:
        return pd.DataFrame()

    rows = []
    for t1, t2 in combinations(sorted(tool_positions.keys()), 2):
        p1, p2 = tool_positions[t1], tool_positions[t2]
        inter = len(p1 & p2)
        union = len(p1 | p2)
        rows.append(
            {
                "tool_1": t1,
                "tool_2": t2,
                "n_tool_1": len(p1),
                "n_tool_2": len(p2),
                "n_intersection": inter,
                "n_union": union,
                "jaccard": inter / union if union else 0.0,
                "pct_1_in_2": inter / len(p1) if p1 else 0.0,
                "pct_2_in_1": inter / len(p2) if p2 else 0.0,
            }
        )

    return pd.DataFrame(rows)


def get_tool_unique_positions(
    tool_outputs: Dict[str, pd.DataFrame],
    tool_name: str,
    threshold: float = None,
    reference: str = None,
) -> pd.DataFrame:
    if tool_name not in tool_outputs:
        raise ValueError(f"Tool '{tool_name}' not found")

    target = get_tool_positions(tool_outputs[tool_name], threshold, reference)
    others = set()
    for other, df in tool_outputs.items():
        if other == tool_name or df.empty:
            continue
        others.update(get_tool_positions(df, threshold, reference))

    unique_keys = target - others
    if not unique_keys:
        return pd.DataFrame(columns=tool_outputs[tool_name].columns)

    mask = [
        (str(ref), int(pos)) in unique_keys
        for ref, pos in zip(
            tool_outputs[tool_name]["reference"].astype(str),
            pd.to_numeric(tool_outputs[tool_name]["position"], errors="coerce").fillna(-1).astype(int),
        )
    ]
    return tool_outputs[tool_name][mask].copy()


def get_consensus_positions(
    tool_outputs: Dict[str, pd.DataFrame],
    min_tools: int = 2,
    threshold: float = None,
    reference: str = None,
) -> pd.DataFrame:
    tool_positions = {
        t: get_tool_positions(df, threshold, reference)
        for t, df in tool_outputs.items()
        if not df.empty
    }
    if not tool_positions:
        return pd.DataFrame()

    all_keys: Set[Key] = set().union(*tool_positions.values())

    rows = []
    for ref, pos in sorted(all_keys):
        det = sorted([t for t, keys in tool_positions.items() if (ref, pos) in keys])
        if len(det) >= min_tools:
            rows.append(
                {
                    "reference": ref,
                    "position": pos,
                    "n_tools": len(det),
                    "tools": ",".join(det),
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["n_tools", "reference", "position"], ascending=[False, True, True])
    return out


def summarize_tool_comparison(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame = None,
    threshold: float = None,
    reference: str = None,
) -> Dict:
    comparison = compare_tools(tool_outputs, threshold, reference)

    summary = {
        "n_tools": comparison.n_tools,
        "tools": comparison.tool_names,
        "n_union": comparison.n_union,
        "n_intersection": comparison.n_intersection,
    }

    for tool, n in comparison.n_positions_per_tool.items():
        summary[f"n_{tool}"] = n

    for n, cnt in sorted(comparison.positions_by_n_tools.items()):
        summary[f"detected_by_{n}_tools"] = cnt

    if comparison.pairwise_jaccard:
        summary["mean_pairwise_jaccard"] = float(np.mean(list(comparison.pairwise_jaccard.values())))

    if ground_truth is not None:
        rank_df = rank_tools(tool_outputs, ground_truth, "auprc", reference)
        if not rank_df.empty:
            summary["best_tool"] = rank_df.iloc[0]["tool"]
            summary["best_auprc"] = float(rank_df.iloc[0]["auprc"])
            summary["tool_ranking"] = rank_df[["tool", "auprc", "f1"]].to_dict("records")

    return summary


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compare tools")
    parser.add_argument("--tools", nargs="+", required=True, help="tool:path")
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--reference")
    parser.add_argument("--output", "-o")
    args = parser.parse_args()

    from parse_outputs import parse_tool_output

    outputs = {}
    for item in args.tools:
        tool, path = item.split(":", 1)
        outputs[tool] = parse_tool_output(tool, path)

    comp = compare_tools(outputs, args.threshold, args.reference)
    print(f"Tools: {comp.tool_names}")
    print(f"Union: {comp.n_union}")
    print(f"Intersection: {comp.n_intersection}")

    if args.output:
        generate_upset_data(outputs, args.threshold, args.reference).to_csv(
            f"{args.output}_upset.csv", index=False
        )
        calculate_pairwise_agreement(outputs, args.threshold, args.reference).to_csv(
            f"{args.output}_pairwise.csv", index=False
        )
