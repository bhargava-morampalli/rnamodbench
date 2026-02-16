#!/usr/bin/env python3
"""
Replicate analysis for RNA modification tool outputs.

All concordance/consensus calculations are reference-aware:
keys are (reference, position), not position-only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

Key = Tuple[str, int]


@dataclass
class ReplicateStats:
    tool: str
    reference: str
    n_replicates: int
    n_sites_per_replicate: List[int]
    n_consensus_sites: int
    mean_jaccard: float
    std_jaccard: float
    min_jaccard: float
    max_jaccard: float
    pairwise_jaccard: Dict[Tuple[str, str], float]

    def to_dict(self) -> Dict:
        return {
            "tool": self.tool,
            "reference": self.reference,
            "n_replicates": self.n_replicates,
            "n_sites_rep1": self.n_sites_per_replicate[0] if self.n_sites_per_replicate else 0,
            "n_sites_rep2": self.n_sites_per_replicate[1] if len(self.n_sites_per_replicate) > 1 else 0,
            "n_sites_rep3": self.n_sites_per_replicate[2] if len(self.n_sites_per_replicate) > 2 else 0,
            "n_consensus_sites": self.n_consensus_sites,
            "mean_jaccard": self.mean_jaccard,
            "std_jaccard": self.std_jaccard,
            "min_jaccard": self.min_jaccard,
            "max_jaccard": self.max_jaccard,
        }


def _filter_reference(df: pd.DataFrame, reference: Optional[str]) -> pd.DataFrame:
    if reference is None or df.empty or "reference" not in df.columns:
        return df
    ref_l = reference.strip().lower()
    return df[df["reference"].astype(str).str.lower() == ref_l].copy()


def _to_key_set(df: pd.DataFrame) -> Set[Key]:
    if df.empty:
        return set()
    tmp = df.copy()
    tmp["position"] = pd.to_numeric(tmp["position"], errors="coerce")
    tmp = tmp[tmp["position"].notna()]
    return set(zip(tmp["reference"].astype(str), tmp["position"].astype(int)))


def get_significant_positions(
    df: pd.DataFrame,
    threshold: float = None,
    score_type: str = None,
    top_n: int = None,
) -> Set[Key]:
    if df.empty:
        return set()

    sub = df.copy()

    if "_imputed" in sub.columns:
        sub = sub[~sub["_imputed"]]

    if score_type is None and "score_type" in sub.columns and not sub["score_type"].empty:
        score_type = str(sub["score_type"].iloc[0]).lower()

    if threshold is not None:
        scores = pd.to_numeric(sub["score"], errors="coerce")
        if score_type in {"pvalue", "fdr"}:
            sub = sub[scores <= threshold]
        else:
            sub = sub[scores >= threshold]

    if top_n is not None and len(sub) > top_n:
        score_num = pd.to_numeric(sub["score"], errors="coerce")
        if score_type in {"pvalue", "fdr"}:
            sub = sub.loc[score_num.nsmallest(top_n).index]
        else:
            sub = sub.loc[score_num.nlargest(top_n).index]

    return _to_key_set(sub)


def jaccard_similarity(set1: Set[Key], set2: Set[Key]) -> float:
    if not set1 and not set2:
        return 1.0
    if not set1 or not set2:
        return 0.0
    inter = len(set1 & set2)
    union = len(set1 | set2)
    return inter / union if union else 0.0


def consensus_calling(
    tool_output: pd.DataFrame,
    min_replicates: int = 2,
    threshold: float = None,
    reference: str = None,
) -> pd.DataFrame:
    df = _filter_reference(tool_output.copy(), reference)

    if df.empty:
        return pd.DataFrame(columns=["reference", "position", "n_replicates", "mean_score", "replicates"])

    reps = sorted(df["replicate"].astype(str).unique()) if "replicate" in df.columns else ["rep1"]
    score_type = str(df["score_type"].iloc[0]).lower() if "score_type" in df.columns and not df.empty else None

    rep_positions: Dict[str, Set[Key]] = {}
    for rep in reps:
        rep_df = df[df["replicate"].astype(str) == rep] if "replicate" in df.columns else df
        rep_positions[rep] = get_significant_positions(rep_df, threshold, score_type)

    all_keys: Set[Key] = set().union(*rep_positions.values()) if rep_positions else set()
    rows = []

    for ref, pos in sorted(all_keys):
        hit_reps = [rep for rep, keys in rep_positions.items() if (ref, pos) in keys]
        if len(hit_reps) < min_replicates:
            continue

        hit_mask = (
            (df["reference"].astype(str) == ref)
            & (pd.to_numeric(df["position"], errors="coerce") == pos)
            & (df["replicate"].astype(str).isin(hit_reps) if "replicate" in df.columns else True)
        )
        scores = pd.to_numeric(df.loc[hit_mask, "score"], errors="coerce")

        rows.append(
            {
                "reference": ref,
                "position": int(pos),
                "n_replicates": len(hit_reps),
                "mean_score": float(np.nanmean(scores)) if not scores.empty else np.nan,
                "replicates": ",".join(hit_reps),
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["reference", "position"]).reset_index(drop=True)
    return out


def calculate_concordance(
    tool_output: pd.DataFrame,
    threshold: float = None,
    reference: str = None,
) -> ReplicateStats:
    tool_name = (
        tool_output["tool"].iloc[0] if "tool" in tool_output.columns and not tool_output.empty else "unknown"
    )

    df = _filter_reference(tool_output.copy(), reference)
    if df.empty:
        return ReplicateStats(
            tool=tool_name,
            reference=reference or "all",
            n_replicates=0,
            n_sites_per_replicate=[],
            n_consensus_sites=0,
            mean_jaccard=0.0,
            std_jaccard=0.0,
            min_jaccard=0.0,
            max_jaccard=0.0,
            pairwise_jaccard={},
        )

    reps = sorted(df["replicate"].astype(str).unique()) if "replicate" in df.columns else ["rep1"]
    score_type = str(df["score_type"].iloc[0]).lower() if "score_type" in df.columns else None

    rep_keys: Dict[str, Set[Key]] = {}
    n_sites = []
    for rep in reps:
        rep_df = df[df["replicate"].astype(str) == rep] if "replicate" in df.columns else df
        keys = get_significant_positions(rep_df, threshold, score_type)
        rep_keys[rep] = keys
        n_sites.append(len(keys))

    pairwise = {}
    jvals = []
    for r1, r2 in combinations(reps, 2):
        j = jaccard_similarity(rep_keys[r1], rep_keys[r2])
        pairwise[(r1, r2)] = j
        jvals.append(j)

    if jvals:
        mean_j = float(np.mean(jvals))
        std_j = float(np.std(jvals, ddof=1)) if len(jvals) > 1 else 0.0
        min_j = float(np.min(jvals))
        max_j = float(np.max(jvals))
    else:
        mean_j = std_j = min_j = max_j = 0.0

    consensus = consensus_calling(df, min_replicates=2, threshold=threshold, reference=reference)

    return ReplicateStats(
        tool=tool_name,
        reference=reference or "all",
        n_replicates=len(reps),
        n_sites_per_replicate=n_sites,
        n_consensus_sites=len(consensus),
        mean_jaccard=mean_j,
        std_jaccard=std_j,
        min_jaccard=min_j,
        max_jaccard=max_j,
        pairwise_jaccard=pairwise,
    )


def calculate_concordance_all_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    references: List[str] = None,
) -> pd.DataFrame:
    rows = []

    for tool_name, df in tool_outputs.items():
        if df.empty:
            continue

        refs = references
        if refs is None:
            refs = sorted(df["reference"].dropna().astype(str).unique().tolist())

        for ref in refs:
            try:
                rows.append(calculate_concordance(df, threshold, ref).to_dict())
            except Exception as exc:
                logger.error("Concordance failed for %s/%s: %s", tool_name, ref, exc)

    return pd.DataFrame(rows)


def aggregate_replicates(
    tool_output: pd.DataFrame,
    method: str = "mean",
    reference: str = None,
) -> pd.DataFrame:
    df = _filter_reference(tool_output.copy(), reference)
    if df.empty:
        return pd.DataFrame()

    agg_map = {"mean": "mean", "median": "median", "min": "min", "max": "max"}
    if method not in agg_map:
        raise ValueError(f"Unknown aggregation method: {method}")

    out = (
        df.groupby(["reference", "position"], as_index=False)
        .agg(score=("score", agg_map[method]), n_replicates=("replicate", "count"))
        .rename(columns={"score": f"score_{method}"})
    )

    out["tool"] = df["tool"].iloc[0] if "tool" in df.columns else "unknown"
    out["score_type"] = df["score_type"].iloc[0] if "score_type" in df.columns else "unknown"
    out["aggregation"] = method
    return out


def filter_by_replicate_count(
    tool_output: pd.DataFrame,
    min_replicates: int = 2,
    threshold: float = None,
) -> pd.DataFrame:
    consensus = consensus_calling(tool_output, min_replicates=min_replicates, threshold=threshold)
    if consensus.empty:
        return pd.DataFrame(columns=tool_output.columns)

    allowed = set(zip(consensus["reference"].astype(str), consensus["position"].astype(int)))
    mask = [
        (str(ref), int(pos)) in allowed
        for ref, pos in zip(
            tool_output["reference"].astype(str),
            pd.to_numeric(tool_output["position"], errors="coerce").fillna(-1).astype(int),
        )
    ]
    return tool_output[mask].copy()


def calculate_replicate_correlation(
    tool_output: pd.DataFrame,
    reference: str = None,
    method: str = "spearman",
) -> pd.DataFrame:
    df = _filter_reference(tool_output.copy(), reference)
    if df.empty or "replicate" not in df.columns:
        return pd.DataFrame()

    reps = sorted(df["replicate"].astype(str).unique())
    if len(reps) < 2:
        return pd.DataFrame()

    pivot = df.pivot_table(
        index=["reference", "position"],
        columns="replicate",
        values="score",
        aggfunc="first",
    )

    rows = []
    for r1, r2 in combinations(reps, 2):
        if r1 not in pivot.columns or r2 not in pivot.columns:
            continue

        shared = pivot[[r1, r2]].dropna()
        if len(shared) < 2:
            continue

        if method == "pearson":
            corr, pval = stats.pearsonr(shared[r1], shared[r2])
        else:
            corr, pval = stats.spearmanr(shared[r1], shared[r2])

        rows.append(
            {
                "replicate_1": r1,
                "replicate_2": r2,
                "correlation": corr,
                "p_value": pval,
                "n_shared": len(shared),
                "method": method,
            }
        )

    return pd.DataFrame(rows)


def _import_calculate_metrics():
    try:
        from .benchmark_metrics import calculate_metrics
    except Exception:
        from benchmark_metrics import calculate_metrics
    return calculate_metrics


def summarize_replicates(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame = None,
    threshold: float = None,
) -> Dict:
    tool_name = tool_output["tool"].iloc[0] if "tool" in tool_output.columns and not tool_output.empty else "unknown"
    reps = sorted(tool_output["replicate"].astype(str).unique()) if "replicate" in tool_output.columns else ["rep1"]

    summary = {"tool": tool_name, "n_replicates": len(reps), "replicates": reps}

    for rep in reps:
        rep_df = tool_output[tool_output["replicate"].astype(str) == rep] if "replicate" in tool_output.columns else tool_output
        summary[f"{rep}_n_positions"] = len(rep_df)
        summary[f"{rep}_score_mean"] = pd.to_numeric(rep_df["score"], errors="coerce").mean()
        summary[f"{rep}_score_std"] = pd.to_numeric(rep_df["score"], errors="coerce").std()

    conc = calculate_concordance(tool_output, threshold)
    summary["mean_jaccard"] = conc.mean_jaccard
    summary["n_consensus_2"] = len(consensus_calling(tool_output, 2, threshold))
    summary["n_consensus_3"] = len(consensus_calling(tool_output, 3, threshold))

    corr_df = calculate_replicate_correlation(tool_output)
    if not corr_df.empty:
        summary["mean_correlation"] = corr_df["correlation"].mean()

    if ground_truth is not None:
        calculate_metrics = _import_calculate_metrics()
        for rep in reps:
            rep_df = tool_output[tool_output["replicate"].astype(str) == rep] if "replicate" in tool_output.columns else tool_output
            m = calculate_metrics(rep_df, ground_truth)
            summary[f"{rep}_auprc"] = m.auprc
            summary[f"{rep}_f1"] = m.f1_optimal

    return summary


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Replicate analysis")
    parser.add_argument("input")
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--threshold", type=float)
    parser.add_argument("--output", "-o")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    stats = calculate_concordance(df, args.threshold)
    consensus = consensus_calling(df, args.min_replicates, args.threshold)

    print(pd.DataFrame([stats.to_dict()]).to_string(index=False))
    print(f"Consensus rows: {len(consensus)}")

    if args.output:
        pd.DataFrame([stats.to_dict()]).to_csv(f"{args.output}_concordance.csv", index=False)
        consensus.to_csv(f"{args.output}_consensus.csv", index=False)
