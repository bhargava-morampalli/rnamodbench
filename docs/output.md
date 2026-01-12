# nf-core/rnamodifications: Output

## Introduction

This document describes the output produced by the pipeline. The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

1. [Mapping](#mapping) - Align reads to 16S and 23S rRNA references
2. [QC Statistics](#qc-statistics) - Generate alignment quality metrics
3. [Signal Processing](#signal-processing) - Prepare signal-level data for modification calling
4. [Modification Calling](#modification-calling) - Detect RNA modifications using multiple tools
5. [Pipeline Information](#pipeline-information) - Report metrics generated during the workflow execution

## Mapping

<details markdown="1">
<summary>Output files</summary>

- `mapping/`
  - `*/*.sorted.bam`: Sorted BAM files aligned to 16S/23S references
  - `*/*.sorted.bam.bai`: BAM index files

</details>

The pipeline aligns reads to both 16S and 23S rRNA references using minimap2, then extracts mapped reads and re-aligns for improved quality.

## QC Statistics

<details markdown="1">
<summary>Output files</summary>

- `qc/flagstat/`
  - `*_flagstat.txt`: Samtools flagstat output with alignment statistics
- `qc/depth/`
  - `*_depth.txt`: Per-position read depth
- `qc/coverage/`
  - `*_coverage.png`: Coverage plots
- `qc/nanoplot/`
  - `*/`: NanoPlot QC reports for each sample

</details>

Quality control metrics are generated for all mapped samples including alignment rates, coverage depth, and read length distributions.

## Signal Processing

<details markdown="1">
<summary>Output files</summary>

- `signal_data/`
  - `fast5_subset/`: Subset of FAST5 files containing mapped reads
  - `single_fast5/`: FAST5 files converted to single-read format
  - `eventalign/`: F5C eventalign output for downstream tools
  - `tombo_resquiggled/`: Tombo-resquiggled FAST5 files

</details>

Signal-level data is prepared for modification calling tools that require raw nanopore signals.

## Modification Calling

The pipeline runs multiple RNA modification detection tools in parallel:

### Tombo

<details markdown="1">
<summary>Output files</summary>

- `modifications/tombo/`
  - `*_de_novo.bed`: BED file with modification statistics
  - `*.tombo.stats`: Raw Tombo statistics files

</details>

Tombo performs de novo detection by comparing native signals to expected unmodified signals.

### Yanocomp

<details markdown="1">
<summary>Output files</summary>

- `modifications/yanocomp/`
  - `*_results.bed`: BED file with modification calls
  - `*_analysis.hdf5`: HDF5 files with detailed results

</details>

Yanocomp uses statistical comparison of current distributions between native and IVT samples.

### Nanocompore

<details markdown="1">
<summary>Output files</summary>

- `modifications/nanocompore/`
  - `*/outSampComp*.tsv`: Tab-separated results with modification statistics
  - `*/outSampComp*.bed`: BED-formatted results

</details>

Nanocompore performs sample comparison using signal-level data from f5c eventalign.

### Xpore

<details markdown="1">
<summary>Output files</summary>

- `modifications/xpore/`
  - `*_diffmod.table`: Differential modification results
  - `dataprep/`: Prepared data for analysis

</details>

Xpore uses a Bayesian approach to identify differential RNA modifications.

### ELIGOS

<details markdown="1">
<summary>Output files</summary>

- `modifications/eligos/`
  - `*_pair_diff_mod.csv`: Pairwise differential modification results
  - `*_combine_per_site.csv`: Per-site combined statistics

</details>

ELIGOS identifies modifications based on error patterns in aligned reads.

### EpiNano-Error

<details markdown="1">
<summary>Output files</summary>

- `modifications/epinano/`
  - `*_error_report.csv`: Error-based modification calls

</details>

EpiNano-Error detects modifications based on systematic base-calling errors.

### DiffErr

<details markdown="1">
<summary>Output files</summary>

- `modifications/differr/`
  - `*_differr.bed`: BED file with differential error rates

</details>

DiffErr identifies modifications by comparing error rates between conditions.

### DRUMMER

<details markdown="1">
<summary>Output files</summary>

- `modifications/drummer/`
  - `*_drummer.tsv`: DRUMMER results with odds ratios

</details>

DRUMMER calculates odds ratios for modification detection.

### JACUSA2

<details markdown="1">
<summary>Output files</summary>

- `modifications/jacusa2/`
  - `*_jacusa2.bed`: BED file with variant calls

</details>

JACUSA2 detects modifications using a variant-calling approach.

## Pipeline Information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html`: Nextflow execution report
  - `execution_timeline_*.html`: Timeline of task execution
  - `execution_trace_*.txt`: Trace file with detailed task information
  - `pipeline_dag_*.html`: DAG visualization of the pipeline
  - `software_versions.yml`: Versions of all software used

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Interpreting results

### Modification significance

Each tool uses different statistical approaches:

| Tool        | Key Metric            | Significance Threshold |
| ----------- | --------------------- | ---------------------- |
| Tombo       | Statistic             | Tool-specific          |
| Yanocomp    | FDR                   | < 0.05                 |
| Nanocompore | GMM_logit_pvalue      | < 0.05                 |
| Xpore       | pval_A_vs_B           | < 0.05                 |
| ELIGOS      | pval, oddR            | pval < 0.05            |
| DiffErr     | FDR                   | < 0.05                 |
| DRUMMER     | pvalue, odds_ratio    | pval < 0.05, OR > 1.5  |
| JACUSA2     | score                 | Tool-specific          |

### Recommended workflow

1. Start with tools that have the best coverage (Tombo, Nanocompore)
2. Cross-reference positions called by multiple tools
3. Apply appropriate significance thresholds
4. Validate top candidates experimentally
