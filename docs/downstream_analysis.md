# Downstream Analysis

## Overview

The downstream analysis module evaluates and compares RNA modification detection tools against a ground truth. It calculates benchmark metrics (AUPRC, AUROC, F1, precision, recall), performs cross-replicate concordance analysis, and generates publication-quality visualizations.

## Supported Tools

| Tool | Score Type | Output Format |
|------|-----------|---------------|
| Tombo | p-value | CSV |
| Yanocomp | FDR | BED |
| Nanocompore | p-value | TSV |
| xPore | probability | diffmod.table |
| ELIGOS | p-value | TXT |
| EpiNano | z-score | CSV |
| DiffErr | neglog10_fdr or g_stat | BED |
| DRUMMER | p-value | TSV |
| JACUSA2 | score | BED |

## Usage

### Via Nextflow pipeline

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --references references.csv \
    --outdir results \
    --run_downstream true \
    --ground_truth known_sites.csv \
    --differr_score_field both \
    -profile singularity
```

### Standalone

```bash
python bin/downstream_analysis/run_analysis.py \
    --input-dir results/modifications \
    --ground-truth known_sites.csv \
    --output-dir downstream_analysis \
    --differr-score-field both
```

### Command-line options

| Option | Description | Default |
|--------|-------------|---------|
| `--input-dir, -i` | Directory with tool outputs | Required |
| `--ground-truth, -g` | CSV/TSV with known modification sites | Optional |
| `--output-dir, -o` | Output directory | Required |
| `--min-replicates` | Minimum replicates for consensus | `2` |
| `--expected-replicates` | Expected replicates for quality check | `3` |
| `--threshold` | Score threshold for significance | None |
| `--differr-score-field` | DiffErr mode: `both`, `g_fdr_neglog10`, or `g_stat` | `both` |
| `--references-csv` | References CSV for position standardization | Optional |
| `--run-id` | Run identifier for collation | Optional |
| `--coverage-label` | Coverage label for collation | Optional |
| `--quality-label` | Quality label for collation | Optional |

### Multi-run collation

Merge results from multiple downstream analysis runs (e.g., different coverage levels):

```bash
# Explicit runs
python bin/downstream_analysis/collate_runs.py \
    --runs results_5x/downstream_analysis results_10x/downstream_analysis \
    --output collated_output

# Auto-discovery
python bin/downstream_analysis/collate_runs.py \
    --runs-root covbench_results \
    --run-glob 'results_*/downstream_analysis' \
    --output collated_output
```

| Option | Description |
|--------|-------------|
| `--runs` | Explicit paths to downstream analysis directories |
| `--runs-root` | Root directory for auto-discovery |
| `--run-glob` | Discovery pattern under `--runs-root` |
| `--output, -o` | Output directory for merged CSVs |
| `--tables` | Subset of tables to collate |
| `--strict` | Non-zero exit on soft issues |

## Output Structure

```
downstream_analysis/
├── metadata/
│   ├── data_quality_report.csv
│   ├── analysis_warnings.txt
│   └── quality_summary.txt
├── by_reference/<ref>/
│   ├── 01_metrics/
│   ├── 02_replicate_analysis/
│   ├── 03_tool_comparison/
│   └── 08_visualizations/
├── collation/
│   ├── metrics_long.csv
│   ├── metrics_summary_long.csv
│   ├── window_metrics_long.csv
│   ├── window_metrics_summary_long.csv
│   ├── lag_metrics_long.csv
│   └── lag_metrics_summary_long.csv
└── differr_gstat/              # When --differr-score-field both
    └── (same structure as above)
```

## Ground Truth Format

CSV or TSV with at minimum `reference` and `position` columns:

```csv
reference,position,modification_type
16S_rRNA,1402,m6A
23S_rRNA,2445,Nm
```

Column names are auto-detected (`reference`/`chr`/`chrom`, `position`/`pos`/`start`).
