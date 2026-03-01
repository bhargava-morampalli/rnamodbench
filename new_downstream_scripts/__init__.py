"""
Downstream Analysis Module for RNA Modification Detection Pipeline

This module provides tools for analyzing and comparing outputs from multiple
RNA modification detection tools, including:
- Output parsing and standardization
- Benchmark metrics (AUPRC, AUROC, F1, precision, recall)
- Cross-replicate analysis and consensus calling
- Cross-tool comparison (UpSet plots, Venn diagrams)
- Visualization suite
"""

from .parse_outputs import (
    parse_tombo,
    parse_yanocomp,
    parse_nanocompore,
    parse_xpore,
    parse_eligos,
    parse_epinano,
    parse_differr,
    parse_drummer,
    parse_jacusa2,
    parse_tool_output,
    load_all_tool_outputs,
)

from .benchmark_metrics import (
    load_ground_truth,
    calculate_metrics,
    calculate_metrics_all_tools,
    get_precision_recall_curve,
    get_roc_curve,
    find_optimal_threshold,
    MetricsResult,
)

from .replicate_analysis import (
    consensus_calling,
    calculate_concordance,
    calculate_concordance_all_tools,
    aggregate_replicates,
    calculate_replicate_correlation,
    ReplicateStats,
)

from .tool_comparison import (
    compare_tools,
    generate_upset_data,
    generate_venn_data,
    rank_tools,
    calculate_pairwise_agreement,
    get_consensus_positions,
    ToolComparisonResult,
)

from .visualization import (
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

from .data_quality import (
    DataQualityReport,
    ToolReplicateStatus,
    Status,
    WarningLevel,
    validate_tool_output,
    check_replicate_completeness,
    check_replicate_overlap,
    check_position_coverage,
    validate_all_outputs,
    generate_quality_summary,
)

from .coverage_analysis import (
    CoverageMetrics,
    CoverageAnalysisResult,
    parse_coverage_value,
    get_coverage_levels,
    calculate_metrics_by_coverage,
    detect_saturation_point,
    calculate_coverage_stability,
    analyze_coverage,
    compare_tools_by_coverage,
    get_optimal_coverage,
    generate_coverage_summary,
)

from .position_standardization import (
    ImputationRecord,
    PositionAvailability,
    StandardizationReport,
    standardize_replicate_positions,
    fill_ground_truth_positions,
    count_within_file_nan,
    standardize_all_tool_outputs,
    generate_availability_report,
)

__version__ = "1.0.0"

__all__ = [
    # Parsing
    "parse_tombo",
    "parse_yanocomp",
    "parse_nanocompore",
    "parse_xpore",
    "parse_eligos",
    "parse_epinano",
    "parse_differr",
    "parse_drummer",
    "parse_jacusa2",
    "parse_tool_output",
    "load_all_tool_outputs",
    # Benchmark metrics
    "load_ground_truth",
    "calculate_metrics",
    "calculate_metrics_all_tools",
    "get_precision_recall_curve",
    "get_roc_curve",
    "find_optimal_threshold",
    "MetricsResult",
    # Replicate analysis
    "consensus_calling",
    "calculate_concordance",
    "calculate_concordance_all_tools",
    "aggregate_replicates",
    "calculate_replicate_correlation",
    "ReplicateStats",
    # Tool comparison
    "compare_tools",
    "generate_upset_data",
    "generate_venn_data",
    "rank_tools",
    "calculate_pairwise_agreement",
    "get_consensus_positions",
    "ToolComparisonResult",
    # Visualization
    "plot_pr_curves",
    "plot_roc_curves",
    "plot_score_distributions",
    "plot_upset",
    "plot_venn",
    "plot_metrics_comparison",
    "plot_replicate_concordance",
    "plot_correlation_heatmap",
    "plot_coverage_performance",
    "plot_coverage_heatmap",
    "plot_coverage_stability",
    "plot_saturation_curves",
    "plot_multi_metric_coverage",
    "generate_html_report",
    # Data Quality
    "DataQualityReport",
    "ToolReplicateStatus",
    "Status",
    "WarningLevel",
    "validate_tool_output",
    "check_replicate_completeness",
    "check_replicate_overlap",
    "check_position_coverage",
    "validate_all_outputs",
    "generate_quality_summary",
    # Coverage Analysis
    "CoverageMetrics",
    "CoverageAnalysisResult",
    "parse_coverage_value",
    "get_coverage_levels",
    "calculate_metrics_by_coverage",
    "detect_saturation_point",
    "calculate_coverage_stability",
    "analyze_coverage",
    "compare_tools_by_coverage",
    "get_optimal_coverage",
    "generate_coverage_summary",
    # Position Standardization
    "ImputationRecord",
    "PositionAvailability",
    "StandardizationReport",
    "standardize_replicate_positions",
    "fill_ground_truth_positions",
    "count_within_file_nan",
    "standardize_all_tool_outputs",
    "generate_availability_report",
]
