# Methods

## Computational pipeline

All analyses were performed using a custom Nextflow (version >=25.04.0) pipeline (rnamodbench v1.0.0dev) implementing ten RNA modification detection tools for systematic benchmarking of Oxford Nanopore direct RNA sequencing (dRNA-seq) data. The pipeline was built using the Nextflow DSL2 framework with the nf-schema plugin (v2.6.1) for parameter validation and samplesheet parsing, and the nf-wave plugin (v1.16.1) for dynamic container provisioning. Software dependencies were managed through Conda environments and pre-built BioContainers, with Seqera Wave used as a fallback to build containers from Conda specifications when pre-built images were unavailable.

## Read alignment and BAM processing

Direct RNA sequencing reads in FASTQ format were aligned to rRNA reference sequences using minimap2 (v2.24; Li, 2018) with the parameters `-ax splice -uf -k14 --secondary=no`, optimised for spliced long-read RNA alignment with secondary alignments disabled. Aligned reads were processed using SAMtools (v1.17; Danecek et al., 2021): SAM-to-BAM conversion was performed with `samtools view` (filtering unmapped reads with flag `-F 4`), followed by coordinate sorting (`samtools sort`), BAM indexing (`samtools index`), and reference indexing (`samtools faidx`). The MD tag was added to BAM files using `samtools calmd` where required by downstream tools. Mapping statistics were generated using `samtools flagstat`, and per-position coverage profiles were calculated using `samtools depth`.

## Raw signal preprocessing

Multi-read FAST5 files were converted to single-read format using the ONT-Fast5-API (v4.1.0) `multi_to_single_fast5` utility, and reads mapping to target references were extracted as subsets using `fast5_subset`. Mapped read identifiers were extracted from BAM files using SeqKit (v2.3.1; Shen et al., 2016). For signal-level analyses requiring event-to-reference alignment, f5c (v1.5; Gamaarachchi et al., 2022) was used to index FAST5 and FASTQ files (`f5c index`) and to generate signal-level event alignments (`f5c eventalign`) with the parameters `--rna --scale-events --signal-index --print-read-names --samples`. A separate f5c eventalign run with the parameters `--rna --scale-events --signal-index` (without `--print-read-names`) was performed to produce output compatible with xPore.

## Quality control

Sequencing quality metrics were assessed using NanoPlot (v1.41.0; De Coster and Rademakers, 2023) applied to aligned BAM files. Per-position coverage profiles were visualised using custom Python scripts with matplotlib (v3.4.3) and seaborn (v0.11.2).

## RNA modification detection

Ten tools were evaluated for their ability to detect RNA modifications by comparing native (modified) dRNA-seq samples against in vitro transcribed (IVT; unmodified) control samples. These tools fall into two broad categories based on their detection strategy: signal-level approaches that analyse raw nanopore current signals, and error-rate-based approaches that exploit systematic basecalling errors introduced by RNA modifications.

### Signal-level modification detection

**Tombo** (v1.5.1; Stoiber et al., 2017) was run using the `level_sample_compare` method, which performs a two-sample Kolmogorov-Smirnov (KS) test at each reference position comparing the distribution of normalised current levels between native and IVT samples. Raw signal data were first re-squiggled to the reference using `tombo resquiggle` with the parameters `--rna --overwrite --num-most-common-errors 5`. Modification detection was performed with `tombo detect_modifications level_sample_compare` using the KS test statistic (`--statistic-type ks`) with a minimum of one test read (`--minimum-test-reads 1`) and p-values stored (`--store-p-value`). Per-position statistics were exported using `tombo text_output`.

**Nanocompore** (v1.0.4; Leger et al., 2021) was applied to collapsed f5c eventalign output. Signal-level event alignments were first collapsed to per-position summary statistics using `nanocompore eventalign_collapse`. Differential modification was then assessed using `nanocompore sampcomp` with logistic regression (`--logit`), a sequence context of 2 bases (5-mer context, `--sequence_context 2`), and a minimum coverage of 1 read per position (`--min_coverage 1`).

**xPore** (v2.1; Pratanwanich et al., 2021) detects differential RNA modifications using a Bayesian Gaussian mixture model (GMM). Eventalign data from f5c were first preprocessed using `xpore dataprep`, and differential modification rates were estimated using `xpore diffmod`.

**Yanocomp** (v0.2; Parker et al., 2022) employs a GMM-based approach with a KS test to identify modified positions. Eventalign data from f5c were converted to HDF5 format using `yanocomp prepare`, and modification detection was performed using `yanocomp analysis`.

### Error-rate-based modification detection

**ELIGOS2** (Jenjaroenpun et al., 2021) identifies RNA modifications by comparing per-position basecalling error patterns between native and IVT BAM files using `eligos2 pair_diff_mod`. The analysis was run with a minimum read depth of 1 (`--min_depth 1`) and maximum depth of 10,000 (`--max_depth 10000`).

**EpiNano** (Liu et al., 2019) detects modifications through differential error analysis. Per-position error metrics (mismatch, insertion, and deletion frequencies) were extracted from native and IVT BAM files using `Epinano_Variants.py`, pinned to GitHub commit `eba4700`. Differential error analysis was performed using `Epinano_DiffErr.R` in two modes: mismatch-only (`-f mis`) and sum-of-errors (`-f sum_err`), the latter using preprocessed data from `Epinano_sumErr.py`. An error rate difference threshold of 0.1 was applied (`-d 0.1`).

**DiffErr** (Parker et al., 2022) compares read-level error patterns between conditions to detect RNA modifications, reporting per-position odds ratios, G-statistics, and p-values. The tool was installed from GitHub (commit `a554d00`) and run with expression thresholds disabled (`--median-expr-threshold 0 --min-expr-threshold 0`).

**DRUMMER** (Depledge et al., 2021) identifies modifications by comparing per-position base fraction distributions between native and IVT samples using a G-test. The tool was installed from GitHub (commit `6683822`) and run in exome mode (`-a exome`), suitable for targeted rRNA sequencing data.

**JACUSA2** (v2.0.4; Piechotta et al., 2022) uses a pileup-based approach to detect positions with significant differences in base composition between conditions. The tool was run in pairwise comparison mode (`call-2`) with minimum coverage of 1 (`-c 1`), minimum score of 0 (`-m 0`), and minimum base quality of 0 (`-q 0`).

**nanoRMS** (Begik et al., 2022) employs a K-nearest neighbours (KNN) classifier to predict modification status from EpiNano per-site error metrics, comparing mismatch frequencies between native and IVT samples. nanoRMS was not included in the benchmarking analysis as it requires a pre-trained model and was not compatible with the comparative framework used here.

### Benchmarking parameter configuration

To enable unbiased benchmarking across all tools, filtering thresholds were set to maximally permissive values so that all tested positions were reported in each tool's output. Specifically, p-value and false discovery rate (FDR) thresholds were set to 1.0, minimum coverage thresholds were set to 0 or 1, and odds ratio thresholds were set to 0 for tools that support these parameters (Tombo, Nanocompore, xPore, Yanocomp, ELIGOS2, EpiNano, DiffErr, DRUMMER, and JACUSA2). nanoDoc reports a continuous score without user-adjustable thresholds and was included using its default output. nanoDoc produced results at 15 of 25 coverage levels (5x–1000x) owing to its internal coverage requirements. This approach allowed receiver operating characteristic (ROC) and precision-recall (PR) curves to be constructed across the full range of score thresholds for each tool.

## Downstream benchmarking analysis

A custom Python analysis module (v1.0.0) was developed to standardise, compare, and benchmark the outputs of all ten modification detection tools. The module requires Python (>=3.10), pandas (>=1.5.0), NumPy (>=1.24.0), scikit-learn (>=1.3.0), matplotlib (>=3.7.0), and seaborn (>=0.12.0).

**Output parsing and standardisation.** Each tool's native output format was parsed into a common schema containing tool name, reference, position (1-based), primary score, score type (p-value, FDR, z-score, probability, or score), and replicate identifier. Positions were standardised against the full reference universe (all positions from 1 to the reference length), with imputed positions tracked via an `_imputed` flag to ensure transparent handling of missing data. Positions reported by a tool but lacking a computable score (e.g. Nanocompore positions where the GMM logistic regression did not converge) were assigned a score of zero, equivalent to unreported positions; all performance metrics were computed at universe scope (all reference positions) so this assignment does not affect the reported results.

**Score harmonisation for AUROC/AUPRC.** Tool-specific output statistics were transformed into a common ranking metric before AUROC/AUPRC computation. For p-value/FDR-like outputs, the transformed score was `-log10(value)`; for z-scores, `abs(z)` was used; score-based outputs were used as reported. The exact per-tool mapping used in this coverage analysis is shown below.

| Tool | Exact statistic used | Transformation for AUROC/AUPRC | Status in this coverage analysis |
| --- | --- | --- | --- |
| Tombo | `stat` (KS-test p-value) | `-log10(p)` | Used |
| Nanocompore | `GMM_logit_pvalue` | `-log10(p)` | Used |
| xPore | `pval_*` (for example, `pval_ivt_vs_native`) | `-log10(p)` | Used |
| nanoDoc | `scoreTotal` | Used as reported | Used (15 coverages only) |
| Yanocomp | BED FDR score | `-log10(FDR)` | Used |
| ELIGOS2 | `adjPval` | `-log10(adjPval)` | Used |
| EpiNano | `delta_sum_err` from delta-sum_err prediction file | Used as reported | Used |
| DiffErr | `g_stat` (G-test statistic) | Used as reported | Used |
| DRUMMER | `G_padj` (fallback `OR_padj`) | `-log10(adjusted p)` | Used |
| JACUSA2 | BED score column | Used as reported | Used |

For EpiNano, both mismatch (`mis`) and sum-of-errors (`sum_err`) outputs were generated upstream. The downstream parser prioritised delta-sum_err prediction files (`*delta-sum_err*.prediction.csv`), falling back to sum-err and then mismatch files when the preferred output was absent. The `delta_sum_err` column was used as the ranking score (higher values indicating greater modification signal). Therefore, AUROC/AUPRC curves were generated for a single EpiNano track in this run.

**Benchmark metrics.** Classification performance was evaluated using known modification sites as ground truth. Area under the precision-recall curve (AUPRC) and area under the receiver operating characteristic curve (AUROC) were calculated using scikit-learn's `average_precision_score` and `roc_auc_score` functions, respectively. F1 score, precision, and recall were computed at optimal thresholds determined from the precision-recall curve. Confusion matrices (true positives, false positives, false negatives, true negatives) were generated at each tool's optimal operating point.

**Replicate analysis.** Concordance between biological or technical replicates was assessed through consensus calling (sites detected in a minimum number of replicates), pairwise Jaccard similarity indices, and replicate score correlation analysis.

**Tool comparison.** Multi-tool agreement was visualised using UpSet plots (upsetplot >=0.8.0) for multi-set intersection analysis and Venn diagrams (matplotlib-venn >=0.11.0) for pairwise and three-way comparisons. Pairwise agreement between tools was quantified using Jaccard similarity indices, and consensus positions detected by a minimum number of tools were identified.

**Coverage analysis.** The relationship between sequencing depth and detection performance was assessed by running the pipeline across multiple coverage levels. Saturation points were identified where performance improvement fell below 1% with additional coverage, and coverage stability was quantified as the inverse of the coefficient of variation across coverage levels.
