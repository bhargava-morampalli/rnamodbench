#!/usr/bin/env python3
"""
Tool Comparison Module for RNA Modification Detection

Compares modification calls across different detection tools, including:
- UpSet plots showing tool agreement
- Venn diagrams for multi-tool overlaps
- Tool ranking by benchmark metrics
- Pairwise tool agreement statistics

Usage:
    from tool_comparison import compare_tools, generate_upset_data

    comparison = compare_tools(tool_outputs, ground_truth)
    upset_data = generate_upset_data(tool_outputs)
"""

import logging
from pathlib import Path
from typing import Union, Dict, List, Tuple, Set, Optional
from dataclasses import dataclass
from itertools import combinations
from collections import defaultdict

import numpy as np
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class ToolComparisonResult:
    """Results from tool comparison."""
    n_tools: int
    tool_names: List[str]
    n_positions_per_tool: Dict[str, int]
    n_union: int
    n_intersection: int
    pairwise_overlap: Dict[Tuple[str, str], int]
    pairwise_jaccard: Dict[Tuple[str, str], float]
    positions_by_n_tools: Dict[int, int]  # Number of positions detected by N tools
    consensus_positions: Set[int]  # Positions detected by all tools

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        result = {
            'n_tools': self.n_tools,
            'tool_names': ','.join(self.tool_names),
            'n_union': self.n_union,
            'n_intersection': self.n_intersection,
        }
        for tool, n in self.n_positions_per_tool.items():
            result[f'n_{tool}'] = n
        for n_tools, count in self.positions_by_n_tools.items():
            result[f'positions_by_{n_tools}_tools'] = count
        return result


def get_tool_positions(
    tool_output: pd.DataFrame,
    threshold: float = None,
    reference: str = None
) -> Set[int]:
    """
    Extract significant positions from tool output.

    Args:
        tool_output: Tool output DataFrame
        threshold: Optional score threshold
        reference: Filter to specific reference

    Returns:
        Set of significant positions
    """
    df = tool_output.copy()

    if df.empty:
        return set()

    # Filter by reference
    if reference:
        df = df[df['reference'].str.contains(reference, case=False, na=False)]

    # Apply threshold
    if threshold is not None:
        score_type = df['score_type'].iloc[0] if 'score_type' in df.columns else 'pvalue'
        if score_type in ['pvalue', 'fdr']:
            df = df[df['score'] <= threshold]
        else:
            df = df[df['score'] >= threshold]

    return set(df['position'].values)


def compare_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    reference: str = None,
    min_tools: int = 2
) -> ToolComparisonResult:
    """
    Compare modification calls across multiple tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        threshold: Optional score threshold for significance
        reference: Filter to specific reference
        min_tools: Minimum tools for consensus

    Returns:
        ToolComparisonResult with comparison statistics
    """
    logger.info(f"Comparing {len(tool_outputs)} tools")

    # Get positions for each tool
    tool_positions = {}
    for tool_name, tool_df in tool_outputs.items():
        if not tool_df.empty:
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
            consensus_positions=set()
        )

    tool_names = sorted(tool_positions.keys())
    n_tools = len(tool_names)

    # Calculate union and intersection
    all_positions = set()
    for positions in tool_positions.values():
        all_positions.update(positions)

    common_positions = set.intersection(*tool_positions.values()) if tool_positions else set()

    # Calculate pairwise overlaps and Jaccard
    pairwise_overlap = {}
    pairwise_jaccard = {}

    for tool1, tool2 in combinations(tool_names, 2):
        pos1 = tool_positions[tool1]
        pos2 = tool_positions[tool2]

        overlap = len(pos1 & pos2)
        union = len(pos1 | pos2)
        jaccard = overlap / union if union > 0 else 0.0

        pairwise_overlap[(tool1, tool2)] = overlap
        pairwise_jaccard[(tool1, tool2)] = jaccard

    # Count positions by number of tools
    positions_by_n_tools = defaultdict(int)
    for pos in all_positions:
        n_detecting = sum(1 for positions in tool_positions.values() if pos in positions)
        positions_by_n_tools[n_detecting] += 1

    # Get consensus positions (detected by >= min_tools)
    consensus_positions = set()
    for pos in all_positions:
        n_detecting = sum(1 for positions in tool_positions.values() if pos in positions)
        if n_detecting >= min_tools:
            consensus_positions.add(pos)

    return ToolComparisonResult(
        n_tools=n_tools,
        tool_names=tool_names,
        n_positions_per_tool={t: len(p) for t, p in tool_positions.items()},
        n_union=len(all_positions),
        n_intersection=len(common_positions),
        pairwise_overlap=pairwise_overlap,
        pairwise_jaccard=pairwise_jaccard,
        positions_by_n_tools=dict(positions_by_n_tools),
        consensus_positions=consensus_positions
    )


def generate_upset_data(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    reference: str = None
) -> pd.DataFrame:
    """
    Generate data for UpSet plot.

    Returns a DataFrame where each row represents a position and
    columns indicate which tools detected it.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        threshold: Optional score threshold
        reference: Filter to specific reference

    Returns:
        DataFrame suitable for UpSet plotting
    """
    # Get positions for each tool
    tool_positions = {}
    for tool_name, tool_df in tool_outputs.items():
        if not tool_df.empty:
            positions = get_tool_positions(tool_df, threshold, reference)
            tool_positions[tool_name] = positions

    if not tool_positions:
        return pd.DataFrame()

    # Get all unique positions
    all_positions = set()
    for positions in tool_positions.values():
        all_positions.update(positions)

    # Build membership matrix
    tool_names = sorted(tool_positions.keys())
    data = []

    for pos in sorted(all_positions):
        row = {'position': pos}
        for tool in tool_names:
            row[tool] = pos in tool_positions[tool]
        data.append(row)

    df = pd.DataFrame(data)

    # Convert boolean columns to int for UpSet
    for tool in tool_names:
        df[tool] = df[tool].astype(int)

    return df


def generate_venn_data(
    tool_outputs: Dict[str, pd.DataFrame],
    tools: List[str] = None,
    threshold: float = None,
    reference: str = None
) -> Dict[str, int]:
    """
    Generate data for Venn diagram (2-4 tools).

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        tools: List of 2-4 tool names to compare (default: first 3)
        threshold: Optional score threshold
        reference: Filter to specific reference

    Returns:
        Dictionary with set sizes for Venn diagram
    """
    # Select tools
    available_tools = [t for t, df in tool_outputs.items() if not df.empty]

    if tools is None:
        tools = available_tools[:3]  # Default to first 3
    else:
        tools = [t for t in tools if t in available_tools]

    if len(tools) < 2 or len(tools) > 4:
        raise ValueError(f"Venn diagram requires 2-4 tools, got {len(tools)}")

    # Get positions
    tool_positions = {}
    for tool in tools:
        positions = get_tool_positions(tool_outputs[tool], threshold, reference)
        tool_positions[tool] = positions

    # Calculate set sizes for matplotlib-venn format
    if len(tools) == 2:
        a, b = tools
        return {
            '10': len(tool_positions[a] - tool_positions[b]),
            '01': len(tool_positions[b] - tool_positions[a]),
            '11': len(tool_positions[a] & tool_positions[b]),
            'labels': (a, b)
        }

    elif len(tools) == 3:
        a, b, c = tools
        pa, pb, pc = tool_positions[a], tool_positions[b], tool_positions[c]
        return {
            '100': len(pa - pb - pc),
            '010': len(pb - pa - pc),
            '001': len(pc - pa - pb),
            '110': len((pa & pb) - pc),
            '101': len((pa & pc) - pb),
            '011': len((pb & pc) - pa),
            '111': len(pa & pb & pc),
            'labels': (a, b, c)
        }

    else:  # 4 tools
        a, b, c, d = tools
        pa, pb, pc, pd = [tool_positions[t] for t in tools]
        # This gets complex - return simplified version
        return {
            'sets': {t: len(tool_positions[t]) for t in tools},
            'intersection_all': len(pa & pb & pc & pd),
            'labels': tuple(tools)
        }


def rank_tools(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame,
    metric: str = 'auprc',
    reference: str = None
) -> pd.DataFrame:
    """
    Rank tools by benchmark metric.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Ground truth DataFrame
        metric: Metric to rank by ('auprc', 'auroc', 'f1')
        reference: Filter to specific reference

    Returns:
        DataFrame with tool rankings
    """
    from .benchmark_metrics import calculate_metrics

    rankings = []

    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            continue

        try:
            metrics = calculate_metrics(tool_df, ground_truth, reference)
            rankings.append({
                'tool': tool_name,
                'auprc': metrics.auprc,
                'auroc': metrics.auroc,
                'f1': metrics.f1_optimal,
                'precision': metrics.precision_optimal,
                'recall': metrics.recall_optimal,
                'n_predictions': metrics.n_predictions
            })
        except Exception as e:
            logger.warning(f"Could not calculate metrics for {tool_name}: {e}")

    df = pd.DataFrame(rankings)

    if not df.empty:
        df = df.sort_values(metric, ascending=False)
        df['rank'] = range(1, len(df) + 1)

    return df


def calculate_pairwise_agreement(
    tool_outputs: Dict[str, pd.DataFrame],
    threshold: float = None,
    reference: str = None
) -> pd.DataFrame:
    """
    Calculate pairwise agreement statistics between all tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        threshold: Optional score threshold
        reference: Filter to specific reference

    Returns:
        DataFrame with pairwise agreement metrics
    """
    # Get positions for each tool
    tool_positions = {}
    for tool_name, tool_df in tool_outputs.items():
        if not tool_df.empty:
            positions = get_tool_positions(tool_df, threshold, reference)
            if positions:
                tool_positions[tool_name] = positions

    if len(tool_positions) < 2:
        return pd.DataFrame()

    tool_names = sorted(tool_positions.keys())
    results = []

    for tool1, tool2 in combinations(tool_names, 2):
        pos1 = tool_positions[tool1]
        pos2 = tool_positions[tool2]

        intersection = len(pos1 & pos2)
        union = len(pos1 | pos2)
        jaccard = intersection / union if union > 0 else 0.0

        # Calculate agreement percentages
        pct_1_in_2 = intersection / len(pos1) if pos1 else 0.0
        pct_2_in_1 = intersection / len(pos2) if pos2 else 0.0

        results.append({
            'tool_1': tool1,
            'tool_2': tool2,
            'n_tool_1': len(pos1),
            'n_tool_2': len(pos2),
            'n_intersection': intersection,
            'n_union': union,
            'jaccard': jaccard,
            'pct_1_in_2': pct_1_in_2,
            'pct_2_in_1': pct_2_in_1
        })

    return pd.DataFrame(results)


def get_tool_unique_positions(
    tool_outputs: Dict[str, pd.DataFrame],
    tool_name: str,
    threshold: float = None,
    reference: str = None
) -> pd.DataFrame:
    """
    Get positions unique to a specific tool.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        tool_name: Tool to get unique positions for
        threshold: Optional score threshold
        reference: Filter to specific reference

    Returns:
        DataFrame with unique positions and their scores
    """
    if tool_name not in tool_outputs:
        raise ValueError(f"Tool '{tool_name}' not found in outputs")

    # Get positions for target tool
    target_positions = get_tool_positions(tool_outputs[tool_name], threshold, reference)

    # Get positions from all other tools
    other_positions = set()
    for other_name, other_df in tool_outputs.items():
        if other_name != tool_name and not other_df.empty:
            positions = get_tool_positions(other_df, threshold, reference)
            other_positions.update(positions)

    # Find unique positions
    unique_positions = target_positions - other_positions

    # Get full data for unique positions
    tool_df = tool_outputs[tool_name]
    unique_df = tool_df[tool_df['position'].isin(unique_positions)].copy()

    return unique_df


def get_consensus_positions(
    tool_outputs: Dict[str, pd.DataFrame],
    min_tools: int = 2,
    threshold: float = None,
    reference: str = None
) -> pd.DataFrame:
    """
    Get positions detected by at least min_tools.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        min_tools: Minimum number of tools
        threshold: Optional score threshold
        reference: Filter to specific reference

    Returns:
        DataFrame with consensus positions
    """
    # Get positions for each tool
    tool_positions = {}
    for tool_name, tool_df in tool_outputs.items():
        if not tool_df.empty:
            positions = get_tool_positions(tool_df, threshold, reference)
            tool_positions[tool_name] = positions

    if not tool_positions:
        return pd.DataFrame()

    # Get all unique positions
    all_positions = set()
    for positions in tool_positions.values():
        all_positions.update(positions)

    # Find consensus positions
    consensus_data = []
    for pos in all_positions:
        detecting_tools = [t for t, p in tool_positions.items() if pos in p]
        n_tools = len(detecting_tools)

        if n_tools >= min_tools:
            consensus_data.append({
                'position': pos,
                'n_tools': n_tools,
                'tools': ','.join(sorted(detecting_tools))
            })

    df = pd.DataFrame(consensus_data)
    if not df.empty:
        df = df.sort_values(['n_tools', 'position'], ascending=[False, True])

    return df


def summarize_tool_comparison(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: pd.DataFrame = None,
    threshold: float = None,
    reference: str = None
) -> Dict:
    """
    Generate comprehensive tool comparison summary.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Optional ground truth DataFrame
        threshold: Optional score threshold
        reference: Filter to specific reference

    Returns:
        Dictionary with summary statistics
    """
    comparison = compare_tools(tool_outputs, threshold, reference)

    summary = {
        'n_tools': comparison.n_tools,
        'tools': comparison.tool_names,
        'n_union': comparison.n_union,
        'n_intersection': comparison.n_intersection,
    }

    # Per-tool counts
    for tool, n in comparison.n_positions_per_tool.items():
        summary[f'n_{tool}'] = n

    # Positions by number of tools
    for n, count in sorted(comparison.positions_by_n_tools.items()):
        summary[f'detected_by_{n}_tools'] = count

    # Mean pairwise Jaccard
    if comparison.pairwise_jaccard:
        summary['mean_pairwise_jaccard'] = np.mean(list(comparison.pairwise_jaccard.values()))

    # Rankings if ground truth provided
    if ground_truth is not None:
        rankings = rank_tools(tool_outputs, ground_truth, 'auprc', reference)
        if not rankings.empty:
            summary['best_tool'] = rankings.iloc[0]['tool']
            summary['best_auprc'] = rankings.iloc[0]['auprc']
            summary['tool_ranking'] = rankings[['tool', 'auprc', 'f1']].to_dict('records')

    return summary


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Compare RNA modification detection tools'
    )
    parser.add_argument('--tools', nargs='+', required=True,
                        help='Tool names and their output files (format: tool:file)')
    parser.add_argument('--threshold', type=float, help='Score threshold')
    parser.add_argument('--reference', help='Filter to specific reference')
    parser.add_argument('--output', '-o', help='Output file prefix')

    args = parser.parse_args()

    # Parse tool:file pairs
    from parse_outputs import parse_tool_output

    tool_outputs = {}
    for item in args.tools:
        tool, filepath = item.split(':')
        tool_outputs[tool] = parse_tool_output(tool, filepath)

    # Compare tools
    comparison = compare_tools(tool_outputs, args.threshold, args.reference)

    print(f"\nTool Comparison Summary:")
    print(f"  Tools: {', '.join(comparison.tool_names)}")
    print(f"  Union: {comparison.n_union} positions")
    print(f"  Intersection: {comparison.n_intersection} positions")
    print(f"\nPositions per tool:")
    for tool, n in comparison.n_positions_per_tool.items():
        print(f"  {tool}: {n}")
    print(f"\nPositions by number of detecting tools:")
    for n, count in sorted(comparison.positions_by_n_tools.items()):
        print(f"  Detected by {n} tools: {count}")

    if args.output:
        # Save upset data
        upset_df = generate_upset_data(tool_outputs, args.threshold, args.reference)
        upset_df.to_csv(f"{args.output}_upset.csv", index=False)

        # Save pairwise agreement
        pairwise = calculate_pairwise_agreement(tool_outputs, args.threshold, args.reference)
        pairwise.to_csv(f"{args.output}_pairwise.csv", index=False)

        print(f"\nSaved results to {args.output}_*.csv")
