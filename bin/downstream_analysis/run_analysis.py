#!/usr/bin/env python3
"""
Main Analysis Script for RNA Modification Detection Pipeline

Orchestrates all downstream analyses including:
- Output parsing from 9 modification detection tools
- Benchmark metrics calculation (AUPRC, AUROC, F1)
- Cross-replicate analysis and consensus calling
- Cross-tool comparison
- Visualization and report generation

Usage:
    python run_analysis.py \\
        --input-dir /path/to/modifications \\
        --ground-truth /path/to/known_sites.csv \\
        --output-dir /path/to/downstream_analysis

    Or provide individual tool outputs:
    python run_analysis.py \\
        --tombo /path/to/tombo.csv \\
        --drummer /path/to/drummer_dir \\
        --ground-truth /path/to/known_sites.csv \\
        --output-dir /path/to/downstream_analysis
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from parse_outputs import (
    parse_tool_output,
    load_all_tool_outputs,
)
from benchmark_metrics import (
    load_ground_truth,
    calculate_metrics,
    calculate_metrics_all_tools,
)
from replicate_analysis import (
    consensus_calling,
    calculate_concordance,
    calculate_concordance_all_tools,
)
from tool_comparison import (
    compare_tools,
    generate_upset_data,
    calculate_pairwise_agreement,
    rank_tools,
    get_consensus_positions,
)
from visualization import (
    plot_pr_curves,
    plot_roc_curves,
    plot_score_distributions,
    plot_upset,
    plot_venn,
    plot_metrics_comparison,
    plot_replicate_concordance,
    plot_correlation_heatmap,
    plot_coverage_performance,
    plot_coverage_heatmap,
    plot_coverage_stability,
    plot_saturation_curves,
    plot_multi_metric_coverage,
    generate_html_report,
)
from data_quality import (
    DataQualityReport,
    validate_all_outputs,
    generate_quality_summary,
)
from coverage_analysis import (
    analyze_coverage,
    compare_tools_by_coverage,
    generate_coverage_summary,
    get_coverage_levels,
)
from position_standardization import (
    standardize_all_tool_outputs,
    StandardizationReport,
    generate_availability_report,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Tool names
SUPPORTED_TOOLS = [
    'tombo', 'yanocomp', 'nanocompore', 'xpore', 'eligos',
    'epinano', 'differr', 'drummer', 'jacusa2'
]


def load_tool_outputs_from_args(args: argparse.Namespace) -> Dict[str, pd.DataFrame]:
    """Load tool outputs from command line arguments."""
    tool_outputs = {}

    for tool in SUPPORTED_TOOLS:
        arg_name = tool.replace('-', '_')
        filepath = getattr(args, arg_name, None)

        if filepath:
            filepath = Path(filepath)
            if filepath.exists():
                try:
                    df = parse_tool_output(tool, filepath)
                    if not df.empty:
                        tool_outputs[tool] = df
                        logger.info(f"Loaded {len(df)} positions from {tool}")
                except Exception as e:
                    logger.error(f"Failed to load {tool}: {e}")

    return tool_outputs


def run_analysis(
    tool_outputs: Dict[str, pd.DataFrame],
    ground_truth: Optional[pd.DataFrame],
    output_dir: Path,
    references: list = None,
    min_replicates: int = 2,
    score_threshold: float = None,
    expected_replicates: int = 3
) -> None:
    """
    Run complete downstream analysis.

    Args:
        tool_outputs: Dictionary mapping tool names to output DataFrames
        ground_truth: Ground truth DataFrame (optional)
        output_dir: Output directory
        references: List of references to analyze
        min_replicates: Minimum replicates for consensus
        score_threshold: Score threshold for significance
        expected_replicates: Expected number of replicates per tool (for quality checks)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create subdirectories
    metrics_dir = output_dir / 'metrics'
    comparison_dir = output_dir / 'tool_comparison'
    replicate_dir = output_dir / 'replicate_analysis'
    viz_dir = output_dir / 'visualizations'

    for d in [metrics_dir, comparison_dir, replicate_dir, viz_dir]:
        d.mkdir(exist_ok=True)

    logger.info(f"Running analysis with {len(tool_outputs)} tools")
    logger.info(f"Output directory: {output_dir}")

    # =========================================================================
    # 0. Data Quality Validation
    # =========================================================================
    logger.info("Validating data quality...")

    quality_report = validate_all_outputs(tool_outputs, expected_replicates=expected_replicates)

    # Save quality report and warnings
    quality_report.save(output_dir / 'data_quality_report.csv')
    quality_report.save_warnings(output_dir / 'analysis_warnings.txt')

    # Generate and save quality summary
    quality_summary = generate_quality_summary(quality_report)
    with open(output_dir / 'quality_summary.txt', 'w') as f:
        f.write(quality_summary)

    # Log quality summary
    logger.info(f"Data Quality Summary:\n{quality_report.summary()}")

    # Log any critical warnings
    critical_warnings = [w for w in quality_report.warnings if w.level.value == 'ERROR']
    if critical_warnings:
        for warning in critical_warnings:
            logger.error(f"[{warning.tool}] {warning.message}")

    # Filter out tools with FAILED status (keep PARTIAL and OK)
    valid_tool_outputs = {}
    for tool_name, tool_df in tool_outputs.items():
        tool_status = next(
            (s for s in quality_report.tool_statuses if s.tool == tool_name),
            None
        )
        if tool_status and tool_status.status.value == 'FAILED':
            logger.warning(f"Excluding {tool_name} from analysis due to FAILED status")
        else:
            valid_tool_outputs[tool_name] = tool_df

    if not valid_tool_outputs:
        logger.error("No valid tool outputs to analyze after quality filtering")
        return

    tool_outputs = valid_tool_outputs
    logger.info(f"Proceeding with {len(tool_outputs)} valid tools")

    # =========================================================================
    # 0.5. Position Standardization
    # =========================================================================
    logger.info("Standardizing positions across replicates and tools...")

    standardized_outputs, std_report = standardize_all_tool_outputs(
        tool_outputs, ground_truth
    )

    # Save standardization reports
    std_report.save(output_dir)

    # Log summary
    logger.info(f"Position Standardization Summary:\n{std_report.summary()}")

    # Use standardized outputs for all downstream analysis
    tool_outputs = standardized_outputs

    # =========================================================================
    # 1. Benchmark Metrics (if ground truth provided)
    # =========================================================================
    if ground_truth is not None and not ground_truth.empty:
        logger.info("Calculating benchmark metrics...")

        # Calculate metrics for all tools
        metrics_df = calculate_metrics_all_tools(tool_outputs, ground_truth, references)

        if not metrics_df.empty:
            metrics_df.to_csv(metrics_dir / 'all_metrics.csv', index=False)
            logger.info(f"Saved metrics to {metrics_dir / 'all_metrics.csv'}")

            # Generate metrics comparison plot
            plot_metrics_comparison(metrics_df, viz_dir / 'metrics_comparison.png')

            # PR and ROC curves
            for ref in (references or [None]):
                ref_suffix = f"_{ref}" if ref else ""
                plot_pr_curves(tool_outputs, ground_truth,
                              viz_dir / f'pr_curves{ref_suffix}.png', ref)
                plot_roc_curves(tool_outputs, ground_truth,
                               viz_dir / f'roc_curves{ref_suffix}.png', ref)

            # Tool ranking
            ranking = rank_tools(tool_outputs, ground_truth)
            ranking.to_csv(metrics_dir / 'tool_ranking.csv', index=False)

    # =========================================================================
    # 2. Replicate Analysis
    # =========================================================================
    logger.info("Analyzing replicates...")

    concordance_results = []
    for tool_name, tool_df in tool_outputs.items():
        if tool_df.empty:
            continue

        # Check if there are multiple replicates
        replicates = tool_df['replicate'].unique()
        if len(replicates) > 1:
            # Calculate concordance
            stats = calculate_concordance(tool_df, score_threshold)
            concordance_results.append(stats.to_dict())

            # Consensus calling
            consensus = consensus_calling(tool_df, min_replicates, score_threshold)
            if not consensus.empty:
                consensus.to_csv(
                    replicate_dir / f'{tool_name}_consensus.csv',
                    index=False
                )

    if concordance_results:
        concordance_df = pd.DataFrame(concordance_results)
        concordance_df.to_csv(replicate_dir / 'concordance.csv', index=False)
        plot_replicate_concordance(concordance_df, viz_dir / 'replicate_concordance.png')

    # =========================================================================
    # 3. Tool Comparison
    # =========================================================================
    logger.info("Comparing tools...")

    # Generate UpSet data
    upset_df = generate_upset_data(tool_outputs, score_threshold)
    if not upset_df.empty:
        upset_df.to_csv(comparison_dir / 'upset_data.csv', index=False)
        try:
            plot_upset(tool_outputs, viz_dir / 'upset_plot.png',
                      threshold=score_threshold)
        except Exception as e:
            logger.warning(f"Could not generate UpSet plot: {e}")

    # Pairwise agreement
    pairwise = calculate_pairwise_agreement(tool_outputs, score_threshold)
    if not pairwise.empty:
        pairwise.to_csv(comparison_dir / 'pairwise_agreement.csv', index=False)
        plot_correlation_heatmap(tool_outputs, viz_dir / 'tool_correlation.png',
                                threshold=score_threshold)

    # Consensus positions across tools
    for min_tools in [2, 3, len(tool_outputs)]:
        if min_tools <= len(tool_outputs):
            consensus = get_consensus_positions(
                tool_outputs, min_tools=min_tools, threshold=score_threshold
            )
            if not consensus.empty:
                consensus.to_csv(
                    comparison_dir / f'consensus_{min_tools}tools.csv',
                    index=False
                )

    # Tool comparison summary
    comparison = compare_tools(tool_outputs, score_threshold)
    summary_data = [{
        'metric': 'n_tools',
        'value': comparison.n_tools
    }, {
        'metric': 'n_union',
        'value': comparison.n_union
    }, {
        'metric': 'n_intersection',
        'value': comparison.n_intersection
    }]
    for tool, n in comparison.n_positions_per_tool.items():
        summary_data.append({'metric': f'n_{tool}', 'value': n})
    pd.DataFrame(summary_data).to_csv(comparison_dir / 'summary.csv', index=False)

    # =========================================================================
    # 4. Score Distributions
    # =========================================================================
    logger.info("Generating score distributions...")
    plot_score_distributions(tool_outputs, viz_dir / 'score_distributions.png')

    # =========================================================================
    # 5. Venn Diagrams (for top tools)
    # =========================================================================
    if len(tool_outputs) >= 2:
        tools_to_compare = list(tool_outputs.keys())[:3]
        try:
            plot_venn(tool_outputs, viz_dir / 'venn_diagram.png',
                     tools=tools_to_compare, threshold=score_threshold)
        except Exception as e:
            logger.warning(f"Could not generate Venn diagram: {e}")

    # =========================================================================
    # 6. Coverage Analysis (if coverage data available)
    # =========================================================================
    coverage_dir = output_dir / 'coverage_analysis'
    has_coverage_data = False

    # Check if any tool has coverage data
    for tool_df in tool_outputs.values():
        if 'coverage' in tool_df.columns and tool_df['coverage'].notna().any():
            has_coverage_data = True
            break

    if has_coverage_data and ground_truth is not None:
        logger.info("Running coverage analysis...")
        coverage_dir.mkdir(exist_ok=True)

        try:
            # Analyze coverage
            coverage_result = analyze_coverage(tool_outputs, ground_truth, references)

            # Save metrics
            coverage_result.save(coverage_dir / 'coverage_metrics.csv')

            # Save comparison table
            comparison_table = compare_tools_by_coverage(coverage_result)
            if not comparison_table.empty:
                comparison_table.to_csv(coverage_dir / 'coverage_comparison.csv')

            # Save summary
            cov_summary = generate_coverage_summary(coverage_result)
            with open(coverage_dir / 'coverage_summary.txt', 'w') as f:
                f.write(cov_summary)

            # Generate visualizations
            cov_metrics_df = coverage_result.to_dataframe()

            if not cov_metrics_df.empty:
                # Performance vs coverage
                plot_coverage_performance(
                    cov_metrics_df,
                    viz_dir / 'coverage_performance.png'
                )

                # Coverage heatmap
                plot_coverage_heatmap(
                    cov_metrics_df,
                    viz_dir / 'coverage_heatmap.png'
                )

                # Saturation curves with marked saturation points
                plot_saturation_curves(
                    cov_metrics_df,
                    viz_dir / 'saturation_curves.png',
                    saturation_points=coverage_result.saturation_points
                )

                # Stability comparison
                plot_coverage_stability(
                    coverage_result.coverage_stability,
                    viz_dir / 'coverage_stability.png'
                )

                # Multi-metric coverage plot
                plot_multi_metric_coverage(
                    cov_metrics_df,
                    viz_dir / 'multi_metric_coverage.png'
                )

            logger.info(f"Coverage analysis complete: {coverage_dir}")

        except Exception as e:
            logger.warning(f"Coverage analysis failed: {e}")
    elif has_coverage_data:
        logger.info("Coverage data found but no ground truth - skipping coverage analysis")

    # =========================================================================
    # 7. HTML Report
    # =========================================================================
    logger.info("Generating HTML report...")
    if ground_truth is not None:
        report_path = generate_html_report(
            tool_outputs, ground_truth, output_dir / 'report'
        )
        logger.info(f"Generated HTML report: {report_path}")

    # =========================================================================
    # Summary
    # =========================================================================
    logger.info("=" * 60)
    logger.info("Analysis Complete!")
    logger.info("=" * 60)
    logger.info(f"Tools analyzed: {len(tool_outputs)}")
    logger.info(f"Output directory: {output_dir}")

    # Quality summary
    n_warnings = len(quality_report.warnings)
    n_errors = len([w for w in quality_report.warnings if w.level.value == 'ERROR'])
    if n_warnings > 0:
        logger.info(f"Data quality: {n_warnings} warnings ({n_errors} errors)")
        logger.info(f"See {output_dir / 'analysis_warnings.txt'} for details")

    if ground_truth is not None and 'metrics_df' in dir() and not metrics_df.empty:
        best_tool = metrics_df.loc[metrics_df['auprc'].idxmax()]
        logger.info(f"Best performing tool: {best_tool['tool']} (AUPRC={best_tool['auprc']:.4f})")

    logger.info(f"Total unique positions: {comparison.n_union}")
    logger.info(f"Positions detected by all tools: {comparison.n_intersection}")


def main():
    parser = argparse.ArgumentParser(
        description='Run downstream analysis for RNA modification detection',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input options
    input_group = parser.add_argument_group('Input Options')
    input_group.add_argument(
        '--input-dir', '-i',
        help='Directory containing all tool outputs (standard pipeline output structure)'
    )

    # Individual tool inputs
    for tool in SUPPORTED_TOOLS:
        input_group.add_argument(
            f'--{tool}',
            help=f'Path to {tool} output file or directory'
        )

    # Ground truth
    parser.add_argument(
        '--ground-truth', '-g',
        help='Path to ground truth CSV/TSV file with known modification sites'
    )

    # Output options
    parser.add_argument(
        '--output-dir', '-o',
        required=True,
        help='Output directory for analysis results'
    )

    # Analysis options
    parser.add_argument(
        '--references', '-r',
        nargs='+',
        help='Filter to specific references (e.g., 16s 23s)'
    )
    parser.add_argument(
        '--min-replicates',
        type=int,
        default=2,
        help='Minimum replicates for consensus calling (default: 2)'
    )
    parser.add_argument(
        '--expected-replicates',
        type=int,
        default=3,
        help='Expected number of replicates per tool for quality checks (default: 3)'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        help='Score threshold for significance (p-value or FDR cutoff)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '--coverage-dirs',
        action='store_true',
        help='Expect {coverage}/tool/ directory structure for coverage analysis'
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Load tool outputs
    tool_outputs = {}

    if args.input_dir:
        # Load from standard directory structure
        input_dir = Path(args.input_dir)
        if not input_dir.exists():
            logger.error(f"Input directory not found: {input_dir}")
            sys.exit(1)

        tool_outputs = load_all_tool_outputs(
            input_dir, SUPPORTED_TOOLS, coverage_dirs=args.coverage_dirs
        )
    else:
        # Load from individual arguments
        tool_outputs = load_tool_outputs_from_args(args)

    if not tool_outputs:
        logger.error("No tool outputs loaded. Provide --input-dir or individual tool paths.")
        sys.exit(1)

    # Load ground truth
    ground_truth = None
    if args.ground_truth:
        ground_truth_path = Path(args.ground_truth)
        if ground_truth_path.exists():
            ground_truth = load_ground_truth(ground_truth_path)
        else:
            logger.warning(f"Ground truth file not found: {ground_truth_path}")

    # Run analysis
    run_analysis(
        tool_outputs=tool_outputs,
        ground_truth=ground_truth,
        output_dir=Path(args.output_dir),
        references=args.references,
        min_replicates=args.min_replicates,
        score_threshold=args.threshold,
        expected_replicates=args.expected_replicates
    )


if __name__ == '__main__':
    main()
