#!/usr/bin/env python3
"""
Replicate Analysis Module for RNA Modification Detection

Analyzes modification calls across biological replicates, including:
- Consensus calling (sites detected in >=N replicates)
- Inter-replicate concordance (Jaccard similarity)
- Reproducibility statistics
- Per-replicate and combined metrics

Usage:
    from replicate_analysis import consensus_calling, calculate_concordance

    consensus_sites = consensus_calling(tool_output, min_replicates=2)
    concordance = calculate_concordance(tool_output)
"""

import logging
from pathlib import Path
from typing import Union, Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ReplicateStats:
    """Statistics for replicate analysis."""
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
        """Convert to dictionary."""
        return {
            'tool': self.tool,
            'reference': self.reference,
            'n_replicates': self.n_replicates,
            'n_sites_rep1': self.n_sites_per_replicate[0] if self.n_sites_per_replicate else 0,
            'n_sites_rep2': self.n_sites_per_replicate[1] if len(self.n_sites_per_replicate) > 1 else 0,
            'n_sites_rep3': self.n_sites_per_replicate[2] if len(self.n_sites_per_replicate) > 2 else 0,
            'n_consensus_sites': self.n_consensus_sites,
            'mean_jaccard': self.mean_jaccard,
            'std_jaccard': self.std_jaccard,
            'min_jaccard': self.min_jaccard,
            'max_jaccard': self.max_jaccard,
        }


def get_significant_positions(
    df: pd.DataFrame,
    threshold: float = None,
    score_type: str = None,
    top_n: int = None
) -> Set[int]:
    """
    Get significant positions from tool output.

    Args:
        df: Tool output DataFrame with 'position' and 'score' columns
        threshold: Score threshold for significance (None = use all)
        score_type: Type of score ('pvalue', 'fdr', 'zscore', etc.)
        top_n: Return only top N positions by score

    Returns:
        Set of significant positions
    """
    if df.empty:
        return set()

    # Get score type from data if not specified
    if score_type is None and 'score_type' in df.columns:
        score_type = df['score_type'].iloc[0]

    positions = df.copy()

    # Apply threshold if specified
    if threshold is not None:
        if score_type in ['pvalue', 'fdr']:
            # Lower is better for p-values
            positions = positions[positions['score'] <= threshold]
        else:
            # Higher is better for other scores
            positions = positions[positions['score'] >= threshold]

    # Apply top_n filter
    if top_n is not None and len(positions) > top_n:
        if score_type in ['pvalue', 'fdr']:
            positions = positions.nsmallest(top_n, 'score')
        else:
            positions = positions.nlargest(top_n, 'score')

    return set(positions['position'].values)


def jaccard_similarity(set1: Set, set2: Set) -> float:
    """
    Calculate Jaccard similarity between two sets.

    Jaccard = |A ∩ B| / |A ∪ B|

    Args:
        set1: First set
        set2: Second set

    Returns:
        Jaccard similarity coefficient (0 to 1)
    """
    if not set1 and not set2:
        return 1.0  # Both empty = perfect agreement
    if not set1 or not set2:
        return 0.0  # One empty = no agreement

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    return intersection / union if union > 0 else 0.0


def consensus_calling(
    tool_output: pd.DataFrame,
    min_replicates: int = 2,
    threshold: float = None,
    reference: str = None
) -> pd.DataFrame:
    """
    Perform consensus calling across replicates.

    A site is called as consensus if it is detected in at least
    `min_replicates` replicates.

    Args:
        tool_output: Tool output DataFrame with 'replicate' column
        min_replicates: Minimum number of replicates for consensus
        threshold: Optional score threshold for significance
        reference: Filter to specific reference (optional)

    Returns:
        DataFrame with consensus sites and their replicate counts
    """
    logger.info(f"Performing consensus calling (min_replicates={min_replicates})")

    df = tool_output.copy()

    # Filter by reference if specified
    if reference:
        df = df[df['reference'].str.contains(reference, case=False, na=False)]

    if df.empty:
        logger.warning("No data for consensus calling")
        return pd.DataFrame(columns=['position', 'reference', 'n_replicates', 'mean_score', 'replicates'])

    # Get unique replicates
    replicates = df['replicate'].unique()
    logger.info(f"Found {len(replicates)} replicates: {replicates}")

    # Get score type for threshold application
    score_type = df['score_type'].iloc[0] if 'score_type' in df.columns else None

    # Collect positions per replicate
    replicate_positions = {}
    for rep in replicates:
        rep_df = df[df['replicate'] == rep]
        positions = get_significant_positions(rep_df, threshold, score_type)
        replicate_positions[rep] = positions
        logger.debug(f"Replicate {rep}: {len(positions)} positions")

    # Count occurrences across replicates
    all_positions = set()
    for positions in replicate_positions.values():
        all_positions.update(positions)

    # Build consensus data
    consensus_data = []
    for pos in all_positions:
        reps_with_pos = [rep for rep, positions in replicate_positions.items() if pos in positions]
        n_reps = len(reps_with_pos)

        if n_reps >= min_replicates:
            # Get scores from replicates that have this position
            scores = df[(df['position'] == pos) & (df['replicate'].isin(reps_with_pos))]['score'].values
            mean_score = np.nanmean(scores) if len(scores) > 0 else np.nan

            # Get reference for this position
            ref = df[df['position'] == pos]['reference'].iloc[0]

            consensus_data.append({
                'position': pos,
                'reference': ref,
                'n_replicates': n_reps,
                'mean_score': mean_score,
                'replicates': ','.join(sorted(reps_with_pos))
            })

    consensus_df = pd.DataFrame(consensus_data)

    if not consensus_df.empty:
        consensus_df = consensus_df.sort_values('position')
        logger.info(f"Found {len(consensus_df)} consensus sites (in >= {min_replicates} replicates)")
    else:
        logger.warning("No consensus sites found")

    return consensus_df


def calculate_concordance(
    tool_output: pd.DataFrame,
    threshold: float = None,
    reference: str = None
) -> ReplicateStats:
    """
    Calculate inter-replicate concordance statistics.

    Args:
        tool_output: Tool output DataFrame with 'replicate' column
        threshold: Optional score threshold for significance
        reference: Filter to specific reference (optional)

    Returns:
        ReplicateStats object with concordance metrics
    """
    tool_name = tool_output['tool'].iloc[0] if 'tool' in tool_output.columns else 'unknown'

    logger.info(f"Calculating concordance for {tool_name}")

    df = tool_output.copy()

    # Filter by reference if specified
    if reference:
        df = df[df['reference'].str.contains(reference, case=False, na=False)]

    if df.empty:
        return ReplicateStats(
            tool=tool_name,
            reference=reference or 'all',
            n_replicates=0,
            n_sites_per_replicate=[],
            n_consensus_sites=0,
            mean_jaccard=0.0,
            std_jaccard=0.0,
            min_jaccard=0.0,
            max_jaccard=0.0,
            pairwise_jaccard={}
        )

    # Get score type
    score_type = df['score_type'].iloc[0] if 'score_type' in df.columns else None

    # Get unique replicates
    replicates = sorted(df['replicate'].unique())
    n_replicates = len(replicates)

    # Get positions per replicate
    replicate_positions = {}
    n_sites_per_replicate = []
    for rep in replicates:
        rep_df = df[df['replicate'] == rep]
        positions = get_significant_positions(rep_df, threshold, score_type)
        replicate_positions[rep] = positions
        n_sites_per_replicate.append(len(positions))

    # Calculate pairwise Jaccard similarities
    pairwise_jaccard = {}
    jaccard_values = []

    for rep1, rep2 in combinations(replicates, 2):
        j = jaccard_similarity(replicate_positions[rep1], replicate_positions[rep2])
        pairwise_jaccard[(rep1, rep2)] = j
        jaccard_values.append(j)

    # Calculate statistics
    if jaccard_values:
        mean_jaccard = np.mean(jaccard_values)
        std_jaccard = np.std(jaccard_values)
        min_jaccard = np.min(jaccard_values)
        max_jaccard = np.max(jaccard_values)
    else:
        mean_jaccard = std_jaccard = min_jaccard = max_jaccard = 0.0

    # Get consensus count
    consensus_df = consensus_calling(df, min_replicates=2, threshold=threshold)
    n_consensus = len(consensus_df)

    return ReplicateStats(
        tool=tool_name,
        reference=reference or 'all',
        n_replicates=n_replicates,
        n_sites_per_replicate=n_sites_per_replicate,
        n_consensus_sites=n_consensus,
        mean_jaccard=mean_jaccard,
        std_jaccard=std_jaccard,
        min_jaccard=min_jaccard,
        max_jaccard=max_jaccard,
        pairwise_jaccard=pairwise_jaccard
    )


def calculate_concordance_all_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    references: List[str] = None
) -> pd.DataFrame:
    """
    Calculate concordance for all tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        threshold: Optional score threshold
        references: List of references to analyze

    Returns:
        DataFrame with concordance metrics for all tools
    """
    results = []

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            continue

        refs_to_analyze = references or [None]
        if references is None:
            # Add per-reference analysis
            refs_to_analyze = tool_df['reference'].unique().tolist() + [None]

        for ref in refs_to_analyze:
            try:
                stats = calculate_concordance(tool_df, threshold, ref)
                results.append(stats.to_dict())
            except Exception as e:
                logger.error(f"Error calculating concordance for {tool_name}/{ref}: {e}")

    return pd.DataFrame(results)


def aggregate_replicates(
    tool_output: pd.DataFrame,
    method: str = 'mean',
    reference: str = None
) -> pd.DataFrame:
    """
    Aggregate scores across replicates.

    Args:
        tool_output: Tool output DataFrame with 'replicate' column
        method: Aggregation method ('mean', 'median', 'min', 'max')
        reference: Filter to specific reference (optional)

    Returns:
        DataFrame with aggregated scores per position
    """
    df = tool_output.copy()

    if reference:
        df = df[df['reference'].str.contains(reference, case=False, na=False)]

    if df.empty:
        return pd.DataFrame()

    # Group by position and reference
    agg_funcs = {
        'mean': 'mean',
        'median': 'median',
        'min': 'min',
        'max': 'max'
    }

    if method not in agg_funcs:
        raise ValueError(f"Unknown aggregation method: {method}")

    # Aggregate
    aggregated = df.groupby(['reference', 'position']).agg({
        'score': agg_funcs[method],
        'replicate': 'count'  # Count number of replicates
    }).reset_index()

    aggregated = aggregated.rename(columns={
        'score': f'score_{method}',
        'replicate': 'n_replicates'
    })

    # Add metadata
    aggregated['tool'] = df['tool'].iloc[0]
    aggregated['score_type'] = df['score_type'].iloc[0] if 'score_type' in df.columns else 'unknown'
    aggregated['aggregation'] = method

    return aggregated


def filter_by_replicate_count(
    tool_output: pd.DataFrame,
    min_replicates: int = 2,
    threshold: float = None
) -> pd.DataFrame:
    """
    Filter tool output to only include positions found in minimum replicates.

    Args:
        tool_output: Tool output DataFrame
        min_replicates: Minimum number of replicates
        threshold: Optional score threshold

    Returns:
        Filtered DataFrame
    """
    # Get consensus sites
    consensus = consensus_calling(tool_output, min_replicates, threshold)

    if consensus.empty:
        return pd.DataFrame(columns=tool_output.columns)

    # Filter original data to consensus positions
    consensus_positions = set(consensus['position'].values)

    filtered = tool_output[tool_output['position'].isin(consensus_positions)].copy()

    return filtered


def calculate_replicate_correlation(
    tool_output: pd.DataFrame,
    reference: str = None,
    method: str = 'spearman'
) -> pd.DataFrame:
    """
    Calculate score correlation between replicates.

    Args:
        tool_output: Tool output DataFrame
        reference: Filter to specific reference (optional)
        method: Correlation method ('pearson' or 'spearman')

    Returns:
        DataFrame with pairwise correlations
    """
    df = tool_output.copy()

    if reference:
        df = df[df['reference'].str.contains(reference, case=False, na=False)]

    if df.empty:
        return pd.DataFrame()

    replicates = sorted(df['replicate'].unique())

    if len(replicates) < 2:
        return pd.DataFrame()

    # Pivot to wide format
    pivot = df.pivot_table(
        index='position',
        columns='replicate',
        values='score',
        aggfunc='first'
    )

    # Calculate correlations
    correlations = []
    for rep1, rep2 in combinations(replicates, 2):
        if rep1 in pivot.columns and rep2 in pivot.columns:
            # Get shared positions
            shared = pivot[[rep1, rep2]].dropna()

            if len(shared) > 1:
                if method == 'pearson':
                    corr, pval = stats.pearsonr(shared[rep1], shared[rep2])
                else:
                    corr, pval = stats.spearmanr(shared[rep1], shared[rep2])

                correlations.append({
                    'replicate_1': rep1,
                    'replicate_2': rep2,
                    'correlation': corr,
                    'p_value': pval,
                    'n_shared': len(shared),
                    'method': method
                })

    return pd.DataFrame(correlations)


def summarize_replicates(
    tool_output: pd.DataFrame,
    ground_truth: pd.DataFrame = None,
    threshold: float = None
) -> Dict:
    """
    Generate comprehensive replicate summary.

    Args:
        tool_output: Tool output DataFrame
        ground_truth: Optional ground truth for benchmark metrics
        threshold: Score threshold

    Returns:
        Dictionary with summary statistics
    """
    tool_name = tool_output['tool'].iloc[0] if 'tool' in tool_output.columns else 'unknown'
    replicates = sorted(tool_output['replicate'].unique())

    summary = {
        'tool': tool_name,
        'n_replicates': len(replicates),
        'replicates': replicates,
    }

    # Per-replicate stats
    for rep in replicates:
        rep_df = tool_output[tool_output['replicate'] == rep]
        summary[f'{rep}_n_positions'] = len(rep_df)
        summary[f'{rep}_score_mean'] = rep_df['score'].mean()
        summary[f'{rep}_score_std'] = rep_df['score'].std()

    # Concordance
    concordance = calculate_concordance(tool_output, threshold)
    summary['mean_jaccard'] = concordance.mean_jaccard
    summary['n_consensus_2'] = len(consensus_calling(tool_output, min_replicates=2, threshold=threshold))
    summary['n_consensus_3'] = len(consensus_calling(tool_output, min_replicates=3, threshold=threshold))

    # Correlation
    correlations = calculate_replicate_correlation(tool_output)
    if not correlations.empty:
        summary['mean_correlation'] = correlations['correlation'].mean()

    # Benchmark metrics if ground truth provided
    if ground_truth is not None:
        from .benchmark_metrics import calculate_metrics

        # Per-replicate metrics
        for rep in replicates:
            rep_df = tool_output[tool_output['replicate'] == rep]
            metrics = calculate_metrics(rep_df, ground_truth)
            summary[f'{rep}_auprc'] = metrics.auprc
            summary[f'{rep}_f1'] = metrics.f1_optimal

        # Consensus metrics
        consensus_df = consensus_calling(tool_output, min_replicates=2, threshold=threshold)
        if not consensus_df.empty:
            # Create pseudo-output for consensus
            consensus_output = pd.DataFrame({
                'position': consensus_df['position'],
                'score': consensus_df['mean_score'],
                'reference': consensus_df['reference'],
                'tool': tool_name,
                'score_type': tool_output['score_type'].iloc[0] if 'score_type' in tool_output.columns else 'unknown'
            })
            consensus_metrics = calculate_metrics(consensus_output, ground_truth)
            summary['consensus_auprc'] = consensus_metrics.auprc
            summary['consensus_f1'] = consensus_metrics.f1_optimal

    return summary


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze replicate data for RNA modification detection'
    )
    parser.add_argument('input', help='Path to combined tool output file')
    parser.add_argument('--tool', required=True, help='Tool name')
    parser.add_argument('--min-replicates', type=int, default=2,
                        help='Minimum replicates for consensus')
    parser.add_argument('--threshold', type=float, help='Score threshold')
    parser.add_argument('--output', '-o', help='Output file prefix')

    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input)

    # Calculate concordance
    stats = calculate_concordance(df, args.threshold)
    print(f"\nReplicate Concordance for {stats.tool}:")
    print(f"  N replicates: {stats.n_replicates}")
    print(f"  Sites per replicate: {stats.n_sites_per_replicate}")
    print(f"  Consensus sites: {stats.n_consensus_sites}")
    print(f"  Mean Jaccard: {stats.mean_jaccard:.4f} +/- {stats.std_jaccard:.4f}")

    # Consensus calling
    consensus = consensus_calling(df, args.min_replicates, args.threshold)
    print(f"\nConsensus sites (>= {args.min_replicates} replicates): {len(consensus)}")

    if args.output:
        # Save results
        pd.DataFrame([stats.to_dict()]).to_csv(f"{args.output}_concordance.csv", index=False)
        consensus.to_csv(f"{args.output}_consensus.csv", index=False)
        print(f"\nSaved results to {args.output}_*.csv")
