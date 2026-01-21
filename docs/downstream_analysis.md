# Downstream Analysis Module Documentation

## Overview

The downstream analysis module provides comprehensive tools for evaluating and comparing outputs from RNA modification detection tools. It calculates benchmark metrics (AUPRC, F1, precision, recall), performs cross-replicate analysis, and generates publication-quality visualizations.

---

## Module Structure

```
bin/downstream_analysis/
├── __init__.py                  # Package exports
├── parse_outputs.py             # Parsers for all 9 tool outputs
├── benchmark_metrics.py         # AUPRC, AUROC, F1, precision, recall
├── replicate_analysis.py        # Consensus calling, concordance
├── tool_comparison.py           # UpSet plots, Venn diagrams, rankings
├── visualization.py             # All plotting functions
├── data_quality.py              # Data quality validation & error handling
├── coverage_analysis.py         # Coverage vs performance analysis
├── position_standardization.py  # Position standardization & NaN filling
└── run_analysis.py              # Main orchestration script

modules/local/downstream_analysis/
├── main.nf                  # Nextflow process
└── environment.yml          # Conda dependencies
```

---

## Supported Tools & Output Formats

| Tool | Output Format | Key Columns | Score Type |
|------|---------------|-------------|------------|
| **TOMBO** | CSV | position, p-value, KS statistic | p-value |
| **Yanocomp** | BED + JSON | position, FDR, KS statistic | FDR |
| **Nanocompore** | TSV (outSampComp_results.tsv) | pos, GMM_logit_pvalue, KS_pvalue | p-value |
| **xPore** | diffmod.table | position, diff_mod_rate, pval | probability |
| **ELIGOS** | TXT (test_paired_diff_mod_result.txt) | start_loc, pval, oddR | p-value |
| **EpiNano** | CSV (*_diff_err.csv) | pos, z-score, diff_error | z-score |
| **DiffErr** | BED | start, pval, fdr, odds_ratio | FDR |
| **DRUMMER** | summary.txt | position, pval, odds_ratio | p-value |
| **JACUSA2** | BED | start, score | score |

---

## Components

### 1. Output Parsing (`parse_outputs.py`)

Standardizes outputs from all 9 tools into a common DataFrame format:

```python
{
    "tool": str,           # Tool name
    "reference": str,      # Reference/chromosome name
    "position": int,       # 1-based position
    "score": float,        # Primary score for ranking (NaN for missing)
    "score_type": str,     # pvalue, fdr, zscore, probability, score
    "pvalue": float,       # P-value if available (NaN for missing)
    "replicate": str,      # Replicate identifier
    "sample_id": str,      # Full sample identifier
    "coverage": str,       # Coverage level (e.g., "100x", "500x") if applicable
    "_imputed": bool       # True if position was added during standardization
}
```

**Key Functions:**
- `parse_tombo()`, `parse_yanocomp()`, `parse_nanocompore()`, etc.
- `parse_tool_output(tool, filepath)` - Universal parser
- `load_all_tool_outputs(output_dir)` - Load from standard directory structure

### 2. Benchmark Metrics (`benchmark_metrics.py`)

Calculates performance metrics against known modification sites (ground truth).

**Metrics:**
- **AUPRC** (Area Under Precision-Recall Curve) - preferred for imbalanced data
- **AUROC** (Area Under ROC Curve)
- **F1 Score** at optimal threshold
- **Precision/Recall** at various thresholds
- **Optimal threshold** determination (maximizes F1)

**Key Functions:**
```python
# Load ground truth (CSV/TSV with reference, position columns)
ground_truth = load_ground_truth("known_sites.csv")

# Calculate metrics for one tool
metrics = calculate_metrics(tool_output, ground_truth, reference="16s")
# Returns: MetricsResult with auprc, auroc, f1_optimal, etc.

# Calculate for all tools
metrics_df = calculate_metrics_all_tools(tool_outputs, ground_truth)

# Get curve data for plotting
precision, recall, thresholds = get_precision_recall_curve(tool_output, ground_truth)
fpr, tpr, thresholds = get_roc_curve(tool_output, ground_truth)
```

**Ground Truth Format:**
```csv
reference,position,modification_type,strand
16S_rRNA,1402,m6A,+
23S_rRNA,2445,Nm,+
```

### 3. Replicate Analysis (`replicate_analysis.py`)

Analyzes reproducibility across biological replicates.

**Features:**
- **Consensus calling**: Sites detected in ≥N replicates
- **Jaccard concordance**: Pairwise similarity between replicates
- **Replicate correlation**: Score correlation across replicates

**Key Functions:**
```python
# Consensus calling (sites in ≥2 of 3 replicates)
consensus = consensus_calling(tool_output, min_replicates=2)
# Returns: DataFrame with position, n_replicates, mean_score, replicates

# Calculate concordance
stats = calculate_concordance(tool_output, threshold=0.05)
# Returns: ReplicateStats with mean_jaccard, std_jaccard, pairwise_jaccard

# Aggregate scores across replicates
aggregated = aggregate_replicates(tool_output, method='mean')
```

### 4. Tool Comparison (`tool_comparison.py`)

Compares modification calls across different detection tools.

**Features:**
- UpSet plot data generation
- Venn diagram data (2-4 tools)
- Pairwise agreement statistics
- Tool ranking by metrics
- Consensus positions across tools

**Key Functions:**
```python
# Compare all tools
comparison = compare_tools(tool_outputs, threshold=0.05)
# Returns: ToolComparisonResult with n_union, n_intersection, pairwise_jaccard

# Generate UpSet data
upset_df = generate_upset_data(tool_outputs, threshold=0.05)

# Pairwise agreement
pairwise = calculate_pairwise_agreement(tool_outputs)

# Rank tools by metric
ranking = rank_tools(tool_outputs, ground_truth, metric='auprc')

# Get consensus positions (detected by ≥N tools)
consensus = get_consensus_positions(tool_outputs, min_tools=3)
```

### 5. Visualization (`visualization.py`)

Generates publication-quality plots.

**Available Plots:**
- PR curves (per tool, overlaid)
- ROC curves (per tool, overlaid)
- Score distributions (histograms per tool)
- UpSet plots (tool agreement)
- Venn diagrams (2-3 tools)
- Metrics comparison bar charts
- Replicate concordance bar charts
- Tool correlation heatmaps
- Coverage vs performance curves
- Coverage heatmaps (tools × coverage levels)
- Saturation curves (with saturation point markers)
- Coverage stability charts
- Multi-metric coverage plots

**Key Functions:**
```python
# PR/ROC curves
plot_pr_curves(tool_outputs, ground_truth, "pr_curves.png")
plot_roc_curves(tool_outputs, ground_truth, "roc_curves.png")

# Tool comparison
plot_upset(tool_outputs, "upset_plot.png")
plot_venn(tool_outputs, "venn.png", tools=['tombo', 'drummer', 'eligos'])
plot_correlation_heatmap(tool_outputs, "heatmap.png")

# Coverage analysis plots
plot_coverage_performance(coverage_metrics, "coverage_performance.png")
plot_coverage_heatmap(coverage_metrics, "coverage_heatmap.png")
plot_saturation_curves(coverage_metrics, "saturation.png", saturation_points=sat_points)
plot_coverage_stability(stability_scores, "stability.png")
plot_multi_metric_coverage(coverage_metrics, "multi_metric.png")

# Generate comprehensive HTML report
report_path = generate_html_report(tool_outputs, ground_truth, "output_dir")
```

### 6. Data Quality Validation (`data_quality.py`)

Validates tool outputs and handles error conditions gracefully.

**Features:**
- Validates file existence and content
- Checks replicate completeness (expected vs. actual)
- Analyzes position overlap between replicates
- Generates warnings and error reports
- Supports "warn but continue" policy for missing data

**Status Types:**
- `OK`: All replicates present and valid
- `PARTIAL`: Some replicates missing (1-2 of 3 expected)
- `EMPTY`: File exists but no data
- `FAILED`: File missing or unparseable
- `MISSING`: Expected file not found

**Key Functions:**
```python
from data_quality import (
    validate_all_outputs,
    DataQualityReport,
    generate_quality_summary
)

# Validate all tool outputs
quality_report = validate_all_outputs(tool_outputs, expected_replicates=3)

# Check status
for status in quality_report.tool_statuses:
    print(f"{status.tool}: {status.status} ({status.n_replicates} replicates)")

# Save reports
quality_report.save("data_quality_report.csv")
quality_report.save_warnings("analysis_warnings.txt")

# Get summary
print(quality_report.summary())
```

**Output Files:**
- `data_quality_report.csv`: Per-tool status report
- `analysis_warnings.txt`: All warnings and errors
- `quality_summary.txt`: Human-readable summary

### 7. Coverage Analysis (`coverage_analysis.py`)

Analyzes tool performance across different sequencing coverage levels.

**Features:**
- Calculates metrics at each coverage level
- Detects saturation points (where performance plateaus)
- Measures coverage stability (consistency across coverages)
- Recommends optimal coverage based on performance/cost tradeoff

**Expected Directory Structure for Coverage Analysis:**
```
output_dir/
├── 100x/
│   ├── tombo/
│   ├── yanocomp/
│   └── ...
├── 500x/
│   ├── tombo/
│   └── ...
└── 1000x/
    └── ...
```

**Key Functions:**
```python
from coverage_analysis import (
    analyze_coverage,
    compare_tools_by_coverage,
    generate_coverage_summary,
    detect_saturation_point,
    get_optimal_coverage
)

# Run complete coverage analysis
coverage_result = analyze_coverage(tool_outputs, ground_truth)

# Get metrics by coverage
metrics_df = coverage_result.to_dataframe()

# Compare tools across coverages (pivot table)
comparison = compare_tools_by_coverage(coverage_result)
print(comparison)
#               100x    500x   1000x
# tombo        0.72    0.85    0.89
# drummer      0.68    0.82    0.88

# Check saturation points
for tool, sat_cov in coverage_result.saturation_points.items():
    print(f"{tool}: saturates at {sat_cov}")

# Check stability scores
for tool, stability in coverage_result.coverage_stability.items():
    print(f"{tool}: stability={stability:.3f}")

# Get optimal coverage recommendations
optimal = get_optimal_coverage(coverage_result, cost_weight=0.5)
print(f"Recommended coverage: {optimal}")
```

**Output Files:**
- `coverage_metrics.csv`: Metrics at each coverage level
- `coverage_comparison.csv`: Pivot table of metrics × coverage
- `coverage_summary.txt`: Text summary with saturation points

### 8. Position Standardization (`position_standardization.py`)

Ensures all tools report the same universe of positions for fair comparison.

**Problem Solved:**
- Different tools may report different positions (some miss low-confidence sites)
- Missing positions create bias in metrics (can't distinguish "not detected" from "detected as negative")
- Replicates may have different position coverage

**Standardization Scenarios:**
1. **Within-file NaN**: Blank or invalid scores converted to NaN
2. **Across-replicate filling**: Positions from any replicate added to all replicates (NaN score)
3. **Ground-truth filling**: Known modification sites added if tool didn't report them (NaN score)

**Tracking & Reporting:**
- `_imputed` column: `True` for positions added during standardization, `False` for native
- Position availability report: Counts native vs. imputed positions per tool/replicate
- Detailed imputation records: Every imputed position logged with source

**Key Functions:**
```python
from position_standardization import (
    standardize_replicate_positions,
    fill_ground_truth_positions,
    standardize_all_tool_outputs,
    generate_availability_report,
    StandardizationReport,
)

# Standardize a single tool's replicates
tool_df = standardize_replicate_positions(tool_df)

# Fill missing ground truth positions
tool_df = fill_ground_truth_positions(tool_df, ground_truth)

# Standardize all tools at once (recommended)
standardized_outputs, report = standardize_all_tool_outputs(
    tool_outputs, ground_truth
)

# Check what was imputed
print(report.summary())
# ============================================================
# POSITION STANDARDIZATION SUMMARY
# ============================================================
# Total native positions: 4500
# Within-file NaN values: 12
# Positions added from other replicates: 234
# Positions added from ground truth: 89
# ...

# Access availability data
availability_df = report.to_availability_dataframe()
imputation_df = report.to_imputation_dataframe()

# Save reports
report.save(output_dir)
```

**Output Files:**
- `position_availability_report.csv`: Per-tool/replicate position counts
- `imputed_positions_detail.csv`: Every imputed position with source

**NaN Handling in Downstream Metrics:**
- NaN scores are **excluded** from AUPRC/AUROC calculations
- Metrics include `n_positions_analyzed` to show how many were used
- This "Exclude + Report" approach ensures metrics aren't skewed by imputed positions

---

## Usage

### Command Line

```bash
# Using standard directory structure
python bin/downstream_analysis/run_analysis.py \
    --input-dir results/modifications \
    --ground-truth known_sites.csv \
    --output-dir downstream_analysis \
    --min-replicates 2 \
    --expected-replicates 3

# Using individual tool outputs
python bin/downstream_analysis/run_analysis.py \
    --tombo results/tombo/output.csv \
    --drummer results/drummer_dir \
    --eligos results/eligos_dir \
    --ground-truth known_sites.csv \
    --output-dir downstream_analysis \
    --threshold 0.05

# With coverage analysis (coverage directory structure)
python bin/downstream_analysis/run_analysis.py \
    --input-dir results/coverage_experiment \
    --coverage-dirs \
    --ground-truth known_sites.csv \
    --output-dir downstream_analysis
```

**Command Line Options:**
| Option | Description |
|--------|-------------|
| `--input-dir, -i` | Directory with tool outputs |
| `--ground-truth, -g` | CSV/TSV with known modification sites |
| `--output-dir, -o` | Output directory for results |
| `--min-replicates` | Minimum replicates for consensus (default: 2) |
| `--expected-replicates` | Expected replicates for quality check (default: 3) |
| `--threshold` | Score threshold for significance |
| `--coverage-dirs` | Enable {coverage}/tool/ directory structure |
| `--references, -r` | Filter to specific references |
| `--verbose, -v` | Enable verbose logging |

### Python API

```python
from downstream_analysis import (
    parse_tool_output,
    load_ground_truth,
    calculate_metrics,
    consensus_calling,
    compare_tools,
    plot_pr_curves,
    generate_html_report
)

# Load data
tombo_df = parse_tool_output('tombo', 'results/tombo.csv')
drummer_df = parse_tool_output('drummer', 'results/drummer_dir')
ground_truth = load_ground_truth('known_sites.csv')

tool_outputs = {'tombo': tombo_df, 'drummer': drummer_df}

# Calculate metrics
metrics = calculate_metrics(tombo_df, ground_truth)
print(f"AUPRC: {metrics.auprc:.4f}")

# Compare tools
comparison = compare_tools(tool_outputs)
print(f"Overlap: {comparison.n_intersection} positions")

# Generate visualizations
plot_pr_curves(tool_outputs, ground_truth, 'pr_curves.png')
generate_html_report(tool_outputs, ground_truth, 'report/')
```

### Nextflow Integration

```nextflow
include { DOWNSTREAM_ANALYSIS } from './modules/local/downstream_analysis/main'

DOWNSTREAM_ANALYSIS(
    tombo_outputs,
    yanocomp_outputs,
    nanocompore_outputs,
    xpore_outputs,
    eligos_outputs,
    epinano_outputs,
    differr_outputs,
    drummer_outputs,
    jacusa2_outputs,
    ground_truth_file
)
```

---

## Output Structure

```
downstream_analysis/
├── data_quality_report.csv           # Per-tool validation status
├── analysis_warnings.txt             # All warnings and errors
├── quality_summary.txt               # Human-readable quality summary
├── position_availability_report.csv  # Native vs imputed positions per tool
├── imputed_positions_detail.csv      # Every imputed position with source
├── metrics/
│   ├── all_metrics.csv           # All tools, all references
│   └── tool_ranking.csv          # Tools ranked by AUPRC
├── replicate_analysis/
│   ├── concordance.csv           # Jaccard concordance per tool
│   └── {tool}_consensus.csv      # Consensus sites per tool
├── tool_comparison/
│   ├── summary.csv               # Comparison summary
│   ├── upset_data.csv            # Data for UpSet plot
│   ├── pairwise_agreement.csv    # Pairwise Jaccard
│   └── consensus_{N}tools.csv    # Positions by N+ tools
├── coverage_analysis/            # (If coverage data available)
│   ├── coverage_metrics.csv      # Metrics per coverage level
│   ├── coverage_comparison.csv   # Tools × coverages pivot table
│   └── coverage_summary.txt      # Text summary
├── visualizations/
│   ├── pr_curves.png
│   ├── roc_curves.png
│   ├── score_distributions.png
│   ├── upset_plot.png
│   ├── venn_diagram.png
│   ├── metrics_comparison.png
│   ├── replicate_concordance.png
│   ├── tool_correlation.png
│   ├── coverage_performance.png  # (If coverage data)
│   ├── coverage_heatmap.png      # (If coverage data)
│   ├── saturation_curves.png     # (If coverage data)
│   ├── coverage_stability.png    # (If coverage data)
│   └── multi_metric_coverage.png # (If coverage data)
└── report/
    ├── report.html               # Comprehensive HTML report
    └── *.png                     # Embedded plots
```

---

## Dependencies

```yaml
# environment.yml
dependencies:
  - python>=3.10
  - pandas>=1.5.0
  - numpy>=1.24.0
  - scikit-learn>=1.3.0
  - matplotlib>=3.7.0
  - seaborn>=0.12.0
  - pip:
    - upsetplot>=0.8.0
    - matplotlib-venn>=0.11.0
```

---

## Key Design Decisions

1. **Score Standardization**: Each tool outputs different score types (p-value, FDR, z-score). For metrics calculation, scores are transformed:
   - P-values/FDR: `-log10(score)` (higher = more significant)
   - Z-scores: `abs(score)` (higher = more significant)

2. **Position Matching**: Default is exact position match. Window-based matching (±N bp) available via `window` parameter.

3. **Replicate Handling**:
   - Per-replicate metrics calculated separately
   - Consensus = sites in ≥2/3 replicates by default
   - Supports both analyses simultaneously

4. **Ground Truth Format**: Flexible CSV/TSV with auto-detection of column names (reference/chr/chrom, position/pos/start).

5. **Error Handling ("Warn but Continue")**:
   - Empty files: Return empty DataFrame, log warning, track in quality report
   - Missing replicates: Warn but continue with available data
   - Missing positions: Use union of all positions as testable universe
   - Different positions between replicates: Report Jaccard similarity, flag if <0.5
   - Failed tools: Exclude from analysis, continue with valid tools

6. **Coverage Analysis**:
   - Supports directory structure: `{coverage}/tool/output.csv`
   - Saturation defined as <1% relative improvement
   - Stability = 1 / (1 + CV) where CV is coefficient of variation
   - Auto-detects coverage from path (100x, 500X, coverage_500, etc.)

7. **Position Standardization**:
   - All tools standardized to report the same universe of positions
   - Missing positions filled with NaN (not zero) to distinguish from detected negatives
   - Three scenarios: within-file NaN, across-replicate filling, ground-truth filling
   - `_imputed` column tracks which positions were added vs. native
   - Position availability report shows per-tool statistics
   - NaN scores excluded from metrics calculation (not treated as 0)

---

## References

- [Nature Methods - Systematic evaluation (2025)](https://www.nature.com/articles/s41592-025-02974-y)
- [sklearn metrics documentation](https://scikit-learn.org/stable/modules/model_evaluation.html)
- [DRS_vRNA_modification_content](https://github.com/lrslab/DRS_vRNA_modification_content)
