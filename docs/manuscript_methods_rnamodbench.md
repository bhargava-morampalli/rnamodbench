# RNAModBench Manuscript Methods Package (Core Pipeline + Downstream)

## Main Methods

### 1. Workflow engine and runtime profile
RNAModBench was executed as a Nextflow DSL2 workflow (`main.nf`) that orchestrates mapping, signal processing, multi-tool RNA modification calling, and optional downstream benchmarking. The representative execution context used for this methods description is the February 17, 2026 run recorded in `results_rasusa/results_100x/pipeline_info/execution_report_2026-02-17_13-26-10.html`, launched from `${PROJECT_DIR}` with profile `singularity`, work directory `/mnt/nvme_work`, and command-line parameters:

```bash
nextflow run main.nf \
  --input ${PROJECT_DIR}/samplesheets_rasusa/samplesheet_rasusa_100x.csv \
  --references ${PROJECT_DIR}/references.csv \
  --outdir ${PROJECT_DIR}/results_rasusa/results_100x \
  --run_downstream true \
  --ground_truth ${PROJECT_DIR}/ground_truth_mod_positions.csv \
  --downstream_run_id run_rasusa_100x \
  --downstream_coverage_label 100x \
  --downstream_quality_label base \
  -w /mnt/nvme_work \
  -profile singularity
```

The run report records Nextflow `25.10.4` (build 11173), Wave enabled (`true`), Fusion disabled (`false`), and successful workflow completion with process-level ignored failures permitted where configured.

### 2. Input design and pairing logic (target + replicate)
The workflow consumes (i) a samplesheet CSV (`sample,fastq,type,replicate,fast5_dir,target`) and (ii) a references CSV (`target,reference`). Samples are dynamically assigned to references by `target` (stored as `meta.rrna`), and each sample is mapped exactly once to its designated reference (no dual mapping).

Input validation enforces:
- every `target` appearing in the samplesheet must be present in `--references`;
- every `target + replicate` must contain both `native` and `ivt` entries;
- pooled FAST5 mode must resolve to exactly one canonical (real-path) FAST5 directory per type (`native`, `ivt`).

These constraints are implemented in `workflows/rnamodbench.nf` before mapping and before any modification caller execution.

### 3. Mapping and QC
Reads are aligned using minimap2, converted/sorted/indexed with samtools, and mapped reads are extracted as FASTQ for signal-level processing. QC is computed on mapped BAMs using:
- `samtools flagstat` (alignment summaries),
- `samtools depth -a -m 0` (per-position depth),
- `coverage_plot.py` (coverage figures from depth tables), and
- `create_feather.py`/NanoPlot environment for per-BAM statistics.

Configured defaults are loaded via `nextflow.config` and `conf/modules.config`, notably `params.minimap2_args='-ax splice -uf -k14 --secondary=no'` and `SAMTOOLS_VIEW` filter `-b -F 4`.

### 4. Signal preparation and event alignment
Mapped-read IDs are extracted with seqkit and used to subset pooled FAST5 archives (`fast5_subset`). Multi-read FAST5 files are converted to single-read FAST5 (`multi_to_single_fast5`) for Tombo. In parallel, f5c indexing/event-alignment is performed for each sample/reference pair.

Two eventalign streams are generated:
- `F5C_EVENTALIGN`: read-name eventalign output for Yanocomp and Nanocompore.
- `F5C_EVENTALIGN_XPORE`: read-index eventalign output for xPore compatibility.

Tombo resquiggling is run against target-specific references prior to Tombo differential comparison.

### 5. Modification-calling families
All callers are paired by key `${meta.rrna}_${meta.replicate}` (native vs IVT):

1. Tombo family
- `tombo detect_modifications level_sample_compare` with `--store-p-value`, `--minimum-test-reads 1`, `--statistic-type ks`.
- Tombo stats are converted to full-length per-reference CSV using Tombo Python API (`tombo_stats.LevelStats`).

2. Eventalign-based callers
- Yanocomp: `prep` followed by `gmmtest` (`--fdr-threshold 1.0`, `--min-ks 0.0`).
- Nanocompore: `eventalign_collapse` then `sampcomp` with `--min_coverage 1 --sequence_context 2 --pvalue_thr 1 --logit`.
- xPore: `dataprep` then `diffmod` using generated config (`readcount_min 1`, `readcount_max 1000000`, `pooling false`, `prefiltering false`).

3. BAM-based callers
- ELIGOS2 `pair_diff_mod` with auto-generated whole-reference BED and permissive thresholds (`min_depth 1`, `max_depth 10000`, `pval 1.0`, `oddR 0`, `esb 0`).
- EpiNano-Error mismatch and sum-error analyses from per-site CSVs.
- DiffErr with `-f 1.0` and expression thresholds disabled (`--median-expr-threshold 0 --min-expr-threshold 0`).
- DRUMMER in exome mode with `-p 1.0`, `-f 0`, and odds-ratio floor protection (non-positive values mapped to `1e-6`).
- JACUSA2 `call-2` with all-sites output (`-A`) and permissive thresholds (`-T 0 -c 1 -m 0 -q 0`) after BAM `calmd` preprocessing.

### 6. Downstream benchmarking and metric collation
When `--run_downstream true`, `modules/local/downstream_analysis/main.nf` executes `bin/downstream_analysis/run_analysis.py` on the `modifications/` directory, with optional ground-truth and references CSV injection and run metadata labels (`run_id`, coverage label, quality label). Default downstream analysis settings in this module are `--min-replicates 2` and `--expected-replicates 3`, with optional threshold parameterization.

### 7. Reproducibility and version tracking
Each process emits a `versions.yml`, and the workflow aggregates entries into `pipeline_info/software_versions.yml`. Cross-run consistency was verified by checksum: all primary coverage runs in `covbench_results/`, `results_rasusa/`, `results_rasusa_aln/`, `results_seqkit/`, and `results_seqtk/` share the same manifest hash (`e25779bda6e4a901f13b83abbd81a8a28d9af71e24a84be196731a3c4b8ee1d5`).

The manifest in `results_100x_finalclean/pipeline_info/software_versions.yml` has a distinct hash (`d252b029...`) and contains malformed version fields; it is excluded from primary version reporting.

Audit-style version reporting was applied:
- runtime value from `software_versions.yml` when meaningful,
- fallback/pin from module conda/container spec when runtime value is `unknown` or malformed,
- both shown when runtime and pin diverge.

## Supplementary Methods (Tool-by-Tool)

### 1. MINIMAP2_ALIGN
- Tool/process name: `MINIMAP2_ALIGN`
- Methodological role: map each sample FASTQ to its target-specific reference sequence.
- Inputs/outputs:
  - Biological input: sample-level native or IVT reads.
  - File input: tuple `[meta, reads.fastq]` plus `reference.fa`.
  - File output: `${prefix}.sam`.
- Exact command pattern and key defaults:
```bash
minimap2 \
  ${params.minimap2_args:-'-ax splice -uf -k14 --secondary=no'} \
  -t ${task.cpus} \
  ${reference} \
  ${reads} > ${prefix}.sam
```
- Version reporting line with provenance: executed runtime `2.24-r1122`; module pin/container `bioconda::minimap2=2.24`, `quay.io/biocontainers/minimap2:2.24--h7132678_1`; sources: `results_rasusa/results_100x/pipeline_info/software_versions.yml`, `modules/local/minimap2_align/main.nf`.

### 2. SAMTOOLS_VIEW
- Tool/process name: `SAMTOOLS_VIEW`
- Methodological role: convert alignment output to BAM and filter unmapped reads.
- Inputs/outputs:
  - Biological input: mapped read alignments.
  - File input: tuple `[meta, *.sam]`.
  - File output: `${prefix}.bam`.
- Exact command pattern and key defaults:
```bash
samtools view -b -F 4 -@ ${task.cpus} -o ${prefix}.bam ${sam}
```
- Version reporting line with provenance: executed runtime `1.17`; module pin/container `bioconda::samtools=1.17`, `quay.io/biocontainers/samtools:1.17--h00cdaf9_0`; sources: runtime manifest + `modules/local/samtools_view/main.nf`.

### 3. SAMTOOLS_SORT
- Tool/process name: `SAMTOOLS_SORT`
- Methodological role: sort BAM alignments for indexing and downstream callers.
- Inputs/outputs:
  - Biological input: mapped read alignments.
  - File input: tuple `[meta, bam_or_sam]` from previous step.
  - File output: `${prefix}_sorted.bam`.
- Exact command pattern and key defaults:
```bash
samtools view -S -b -h ${sam} | samtools sort -o ${prefix}_sorted.bam
```
`conf/modules.config` sets `ext.prefix="${meta.id}.sorted"`.
- Version reporting line with provenance: executed runtime `1.17`; module pin/container `samtools=1.17`; sources: runtime manifest + `modules/local/samtools_sort/main.nf` + `conf/modules.config`.

### 4. SAMTOOLS_INDEX
- Tool/process name: `SAMTOOLS_INDEX`
- Methodological role: index sorted BAM files.
- Inputs/outputs:
  - Biological input: mapped/sorted reads.
  - File input: tuple `[meta, *.bam]`.
  - File output: `*.bai`.
- Exact command pattern and key defaults:
```bash
samtools index -@ ${task.cpus} ${bam}
```
- Version reporting line with provenance: executed runtime `1.17`; module pin/container `samtools=1.17`; sources: runtime manifest + `modules/local/samtools_index/main.nf`.

### 5. EXTRACT_MAPPED_READS
- Tool/process name: `EXTRACT_MAPPED_READS`
- Methodological role: export mapped reads to FASTQ for signal-aware downstream tools.
- Inputs/outputs:
  - Biological input: mapped alignments.
  - File input: tuple `[meta, *.sam]`.
  - File output: `${prefix}.fastq`.
- Exact command pattern and key defaults:
```bash
samtools fastq -F 4 --threads ${task.cpus-1} ${sam} > ${prefix}.fastq
```
- Version reporting line with provenance: executed runtime `1.17`; module pin/container `samtools=1.17`; sources: runtime manifest + `modules/local/extract_mapped_reads/main.nf`.

### 6. SAMTOOLS_FLAGSTAT
- Tool/process name: `SAMTOOLS_FLAGSTAT`
- Methodological role: compute mapping summary statistics.
- Inputs/outputs:
  - Biological input: mapped reads.
  - File input: tuple `[meta, *.bam]`.
  - File output: `${prefix}_flagstat.txt`.
- Exact command pattern and key defaults:
```bash
samtools flagstat ${bam} > ${prefix}_flagstat.txt
```
- Version reporting line with provenance: executed runtime `1.17`; module pin/container `samtools=1.17`; sources: runtime manifest + `modules/local/samtools_flagstat/main.nf`.

### 7. SAMTOOLS_DEPTH
- Tool/process name: `SAMTOOLS_DEPTH`
- Methodological role: compute per-position depth across the full reference.
- Inputs/outputs:
  - Biological input: mapped reads.
  - File input: tuple `[meta, *.bam]`.
  - File output: `${prefix}.txt` depth table.
- Exact command pattern and key defaults:
```bash
samtools depth -a -m 0 ${bam} > ${prefix}.txt
```
- Version reporting line with provenance: executed runtime `1.17`; module pin/container `samtools=1.17`; sources: runtime manifest + `modules/local/samtools_depth/main.nf`.

### 8. COVERAGE_PLOT
- Tool/process name: `COVERAGE_PLOT`
- Methodological role: render coverage plots from depth tables.
- Inputs/outputs:
  - Biological input: depth distribution summaries.
  - File input: tuple `[meta, depth.txt]`.
  - File output: `${prefix}.pdf`.
- Exact command pattern and key defaults:
```bash
coverage_plot.py -f ${depth} -t ${meta.id} -o ${prefix}.pdf
```
- Version reporting line with provenance: executed runtime `python 3.9.23; pandas 1.3.4; matplotlib 3.4.3; seaborn 0.11.2`; module pin `python=3.9, pandas=1.3.4, matplotlib=3.4.3, seaborn=0.11.2`; sources: runtime manifest + `modules/local/coverage_plot/main.nf` + `modules/local/coverage_plot/environment.yml`.

### 9. NANOPLOT_BAM
- Tool/process name: `NANOPLOT_BAM`
- Methodological role: generate per-BAM QC feature tables (feather format) in NanoPlot environment.
- Inputs/outputs:
  - Biological input: mapped BAM files.
  - File input: tuple `[meta, bam]`.
  - File output: `${prefix}.feather`.
- Exact command pattern and key defaults:
```bash
create_feather.py --bam ${bam} --output ${prefix}.feather
```
- Version reporting line with provenance: executed runtime `nanoplot 1.41.0`; module pin/container `bioconda::nanoplot=1.41.0`, `quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0`; sources: runtime manifest + `modules/local/nanoplot_bam/main.nf`.

### 10. EXTRACT_READ_IDS
- Tool/process name: `EXTRACT_READ_IDS`
- Methodological role: extract read identifiers from mapped FASTQ files for FAST5 subsetting.
- Inputs/outputs:
  - Biological input: mapped reads (FASTQ).
  - File input: tuple `[meta, fastq]`.
  - File output: `${prefix}.read_ids.txt`.
- Exact command pattern and key defaults:
```bash
seqkit seq -n -i ${fastq} > ${prefix}.read_ids.txt
```
- Version reporting line with provenance: executed runtime `seqkit v2.3.0`; module pin/container `bioconda::seqkit=2.3.1`, `quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0`; sources: runtime manifest + `modules/local/extract_read_ids/main.nf`; note: runtime/pin minor-version mismatch.

### 11. FAST5_SUBSET
- Tool/process name: `FAST5_SUBSET`
- Methodological role: subset pooled FAST5 directories to reads present in mapped FASTQ files.
- Inputs/outputs:
  - Biological input: read IDs corresponding to mapped reads.
  - File input: tuple `[meta, read_ids.txt, fast5_dir]`.
  - File output: `fast5_subset/` directory.
- Exact command pattern and key defaults:
```bash
fast5_subset \
  --input ${fast5_dir} \
  --read_id_list ${read_ids} \
  --save_path fast5_subset \
  --recursive
```
- Version reporting line with provenance: executed runtime `unknown`; module pin/container `bioconda::ont-fast5-api=4.1.0`, `quay.io/biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0`; sources: runtime manifest + `modules/local/fast5_subset/main.nf`.

### 12. MULTI_TO_SINGLE_FAST5
- Tool/process name: `MULTI_TO_SINGLE_FAST5`
- Methodological role: convert multi-read FAST5 files to single-read FAST5 required by Tombo.
- Inputs/outputs:
  - Biological input: subset FAST5 signals.
  - File input: tuple `[meta, fast5_dir]`.
  - File output: `single_fast5_${meta.id}/`.
- Exact command pattern and key defaults:
```bash
multi_to_single_fast5 \
  --input_path ${fast5_dir} \
  --save_path single_fast5_${meta.id} \
  --recursive
```
- Version reporting line with provenance: executed runtime `unknown`; module pin/container `ont-fast5-api=4.1.0`; sources: runtime manifest + `modules/local/multi_to_single_fast5/main.nf`.

### 13. F5C_INDEX
- Tool/process name: `F5C_INDEX`
- Methodological role: create f5c index/readdb artifacts linking FASTQ entries to FAST5 signal records.
- Inputs/outputs:
  - Biological input: mapped reads and corresponding subset FAST5 directories.
  - File input: tuple `[meta, fast5_dir, fastq]` plus optional sequencing summary (`NO_FILE` in this workflow).
  - File output: `*.index`, `*.index.readdb`, `*.index.gzi`, `*.index.fai`.
- Exact command pattern and key defaults:
```bash
f5c index -t ${task.cpus} --iop ${task.cpus} -d ${fast5_dir} ${fastq}
```
- Version reporting line with provenance: executed runtime `unknown`; module pin/container `bioconda::f5c=1.5`, `quay.io/biocontainers/f5c:1.5--hee927d3_2`; sources: runtime manifest + `modules/local/f5c_index/main.nf`.
- Caveat: `errorStrategy='ignore'` in `conf/modules.config`.

### 14. F5C_EVENTALIGN
- Tool/process name: `F5C_EVENTALIGN`
- Methodological role: generate read-name eventalign output for Nanocompore and Yanocomp.
- Inputs/outputs:
  - Biological input: mapped reads + FAST5 signal + reference.
  - File input: tuple `[meta, fastq, index artifacts, bam, bai, fast5, reference]`.
  - File output: `${prefix}_eventalign.txt`.
- Exact command pattern and key defaults:
```bash
f5c eventalign \
  --rna --scale-events --signal-index --print-read-names --samples \
  -t ${task.cpus} -b ${bam} -g ${reference} -r ${fastq} \
  > ${prefix}_eventalign.txt
```
`conf/modules.config` overrides module fallback defaults with `params.f5c_eventalign_args`.
- Version reporting line with provenance: executed runtime `unknown`; module pin/container `f5c=1.5`; sources: runtime manifest + `modules/local/f5c_eventalign/main.nf` + `conf/modules.config` + `nextflow.config`.
- Caveat: `errorStrategy='ignore'`.

### 15. F5C_EVENTALIGN_XPORE
- Tool/process name: `F5C_EVENTALIGN_XPORE`
- Methodological role: generate read-index eventalign output required by xPore dataprep.
- Inputs/outputs:
  - Biological input: mapped reads + FAST5 signal + reference.
  - File input: same structure as `F5C_EVENTALIGN`.
  - File output: `${prefix}_eventalign_xpore.txt`.
- Exact command pattern and key defaults:
```bash
f5c eventalign \
  --rna --scale-events --signal-index \
  -t ${task.cpus} -b ${bam} -g ${reference} -r ${fastq} \
  > ${prefix}_eventalign_xpore.txt
```
- Version reporting line with provenance: executed runtime `unknown`; module pin/container `f5c=1.5`; sources: runtime manifest + `modules/local/f5c_eventalign_xpore/main.nf` + `conf/modules.config`.
- Caveat: `errorStrategy='ignore'`.

### 16. TOMBO_RESQUIGGLE
- Tool/process name: `TOMBO_RESQUIGGLE`
- Methodological role: resquiggle single-read FAST5 signals to target references prior to Tombo comparison.
- Inputs/outputs:
  - Biological input: single-read FAST5 and target reference.
  - File input: tuple `[meta, single_fast5_dir, reference.fa]`.
  - File output: resquiggled FAST5 directory (in place).
- Exact command pattern and key defaults:
```bash
tombo resquiggle \
  --rna --overwrite --num-most-common-errors 5 \
  --processes ${task.cpus} \
  ${fast5} ${reference}
```
- Version reporting line with provenance: executed runtime `1.5.1`; module pin `ont-tombo=1.5.1` (`python=3.7`, `numpy<1.20`); sources: runtime manifest + `modules/local/tombo_resquiggle/main.nf` + `modules/local/tombo_resquiggle/environment.yml`.
- Caveat: `errorStrategy='ignore'`.

### 17. TOMBO_DETECT_MODIFICATIONS
- Tool/process name: `TOMBO_DETECT_MODIFICATIONS`
- Methodological role: pairwise native-vs-IVT differential testing at the signal level.
- Inputs/outputs:
  - Biological input: grouped resquiggled FAST5 directories (`native`, `ivt`) per target+replicate.
  - File input: tuple `[key, types, fast5_dirs]`.
  - File output: `${prefix}.tombo.stats`.
- Exact command pattern and key defaults:
```bash
tombo detect_modifications level_sample_compare \
  --fast5-basedirs ${native_dir} \
  --alternate-fast5-basedirs ${ivt_dir} \
  --statistics-file-basename ${prefix} \
  --store-p-value \
  --minimum-test-reads 1 \
  --statistic-type ks \
  --processes ${task.cpus}
```
- Version reporting line with provenance: executed runtime `1.5.1`; module pin `ont-tombo=1.5.1`; sources: runtime manifest + `modules/local/tombo_detect_modifications/main.nf`.
- Caveat: `errorStrategy='ignore'`.

### 18. TOMBO_TEXT_OUTPUT
- Tool/process name: `TOMBO_TEXT_OUTPUT`
- Methodological role: convert Tombo `.tombo.stats` into full-reference per-position CSV using Tombo Python API.
- Inputs/outputs:
  - Biological input: differential signal statistics and reference sequence.
  - File input: tuple `[key, statistics_file, reference.fa]`.
  - File output: `${prefix}.csv`.
- Exact command pattern and key defaults:
```python
from tombo import tombo_stats
import pandas as pd
sample_level_stats = tombo_stats.LevelStats(statistic_file)
reg_level_stats = sample_level_stats.get_reg_stats(chrom, '+', 1, seq_length)
pd.DataFrame(reg_level_stats).to_csv(f"{prefix}.csv", index=False)
```
- Version reporting line with provenance: executed runtime `tombo 1.5.1; pandas 1.3.4`; module pin `ont-tombo=1.5.1; pandas=1.3.4`; sources: runtime manifest + `modules/local/tombo_text_output/main.nf` + `modules/local/tombo_text_output/environment.yml`.
- Caveat: `errorStrategy='ignore'`.

### 19. YANOCOMP_PREPARE
- Tool/process name: `YANOCOMP_PREPARE`
- Methodological role: prepare eventalign-derived HDF5 feature files for Yanocomp differential testing.
- Inputs/outputs:
  - Biological input: per-sample eventalign tables.
  - File input: tuple `[meta, eventalign.txt]` plus optional summary (`NO_FILE` by default).
  - File output: `${prefix}.hdf5`.
- Exact command pattern and key defaults:
```bash
yanocomp prep -p ${task.cpus} -e ${eventalign} -h ${prefix}.hdf5
```
- Version reporting line with provenance: executed runtime `unknown`; pinned source `yanocomp git tag v0.2` (pip install from `git+https://github.com/bartongroup/yanocomp.git@v0.2`), `python=3.8`; sources: runtime manifest + `modules/local/yanocomp_prepare/main.nf` + `modules/local/yanocomp_prepare/environment.yml` + `modules/local/yanocomp_prepare/Dockerfile`.

### 20. YANOCOMP_ANALYSIS
- Tool/process name: `YANOCOMP_ANALYSIS`
- Methodological role: perform native-vs-IVT GMM-based differential signal testing.
- Inputs/outputs:
  - Biological input: paired native/IVT HDF5 files by target+replicate.
  - File input: tuple `[meta, native.hdf5, ivt.hdf5]`.
  - File output: `${prefix}.bed`, `${prefix}_sm_preds.json`.
- Exact command pattern and key defaults:
```bash
yanocomp gmmtest \
  --fdr-threshold ${params.yanocomp_fdr_threshold} \
  --min-ks ${params.yanocomp_min_ks} \
  -p ${task.cpus} -n 1 \
  -c ${hdf5_native} -t ${hdf5_ivt} \
  -o ${prefix}.bed -s ${prefix}_sm_preds.json
```
With defaults from `nextflow.config`: `yanocomp_fdr_threshold=1.0`, `yanocomp_min_ks=0.0`.
- Version reporting line with provenance: executed runtime `unknown`; pinned source `yanocomp v0.2`; sources: runtime manifest + `modules/local/yanocomp_analysis/main.nf` + `modules/local/yanocomp_analysis/environment.yml`.
- Caveat: `errorStrategy='ignore'`. Representative run contained one ignored Yanocomp failure (`23s_rep1`).

### 21. NANOCOMPORE_EVENTALIGN_COLLAPSE
- Tool/process name: `NANOCOMPORE_EVENTALIGN_COLLAPSE`
- Methodological role: collapse eventalign reads into Nanocompore-ready condition files.
- Inputs/outputs:
  - Biological input: per-sample eventalign tables.
  - File input: tuple `[meta, eventalign.txt]`.
  - File output: `${prefix}/out_eventalign_collapse.tsv` and index files.
- Exact command pattern and key defaults:
```bash
nanocompore eventalign_collapse -i ${eventalign} -o ${prefix} -t ${task.cpus}
```
- Version reporting line with provenance: executed runtime `unknown`; module pin/container `nanocompore=1.0.4`, `quay.io/biocontainers/nanocompore:1.0.4--pyhdfd78af_0`; sources: runtime manifest + `modules/local/nanocompore_eventalign_collapse/main.nf` + `modules/local/nanocompore_eventalign_collapse/environment.yml`.

### 22. NANOCOMPORE_SAMPCOMP
- Tool/process name: `NANOCOMPORE_SAMPCOMP`
- Methodological role: run differential comparison on collapsed native/IVT eventalign data.
- Inputs/outputs:
  - Biological input: paired native and IVT collapsed eventalign files plus reference.
  - File input: tuple `[key, cond1_label, cond2_label, cond1_dirs, cond2_dirs, reference.fa]`.
  - File output: `${prefix}/outSampComp_results.tsv` and associated files.
- Exact command pattern and key defaults:
```bash
nanocompore sampcomp \
  --file_list1 ${cond1_files} \
  --file_list2 ${cond2_files} \
  --label1 ${condition1_label} --label2 ${condition2_label} \
  --fasta ${fasta} --outpath ${prefix} --nthreads ${task.cpus} \
  --min_coverage ${params.nanocompore_min_coverage} \
  --sequence_context 2 --pvalue_thr 1 --logit
```
Default `nanocompore_min_coverage=1`.
- Version reporting line with provenance: executed runtime `unknown`; module pin/container `nanocompore=1.0.4`; sources: runtime manifest + `modules/local/nanocompore_sampcomp/main.nf` + `conf/modules.config`.
- Caveat: `errorStrategy='ignore'`.

### 23. XPORE_DATAPREP
- Tool/process name: `XPORE_DATAPREP`
- Methodological role: convert xpore-compatible eventalign output into dataprep structures.
- Inputs/outputs:
  - Biological input: per-sample eventalign_xpore tables.
  - File input: tuple `[meta, eventalign_xpore.txt]`.
  - File output: `dataprep_${meta.id}/`.
- Exact command pattern and key defaults:
```bash
xpore dataprep --eventalign ${eventalign} --out_dir dataprep_${meta.id} --n_processes ${task.cpus}
```
- Version reporting line with provenance: executed runtime `2.0`; module pin/container `bioconda::xpore=2.1`, `quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0`; sources: runtime manifest + `modules/local/xpore_dataprep/main.nf` + `modules/local/xpore_dataprep/environment.yml`; note: runtime/pin divergence.

### 24. XPORE_DIFFMOD
- Tool/process name: `XPORE_DIFFMOD`
- Methodological role: estimate differential modification between native and IVT conditions.
- Inputs/outputs:
  - Biological input: paired native and IVT xpore dataprep directories.
  - File input: tuple `[key, types, dataprep_dirs]`.
  - File output: `${key}_diffmod/diffmod.table`.
- Exact command pattern and key defaults:
```bash
cat > config.yml <<YAML
out: ${out_dir}
data:
  native:
    rep1: ${native_dir}
  ivt:
    rep1: ${ivt_dir}
criteria:
  readcount_min: 1
  readcount_max: 1000000
method:
  name: gmm
  pooling: false
  prefiltering: false
YAML

xpore diffmod --config config.yml --n_processes ${task.cpus}
```
- Version reporting line with provenance: executed runtime `2.0`; module pin/container `xpore=2.1`; sources: runtime manifest + `modules/local/xpore_diffmod/main.nf` + `modules/local/xpore_diffmod/environment.yml`; note: runtime/pin divergence.
- Caveat: `errorStrategy='ignore'`.

### 25. ELIGOS_PAIR_DIFF_MOD
- Tool/process name: `ELIGOS_PAIR_DIFF_MOD`
- Methodological role: paired native-vs-IVT differential error-signature analysis on BAMs.
- Inputs/outputs:
  - Biological input: paired native/IVT BAMs and target reference.
  - File input: tuple `[meta, native_bam, native_bai, ivt_bam, ivt_bai, reference.fa]`.
  - File output: `${prefix}/` directory with ELIGOS outputs.
- Exact command pattern and key defaults:
```bash
if [ ! -f ${reference}.fai ]; then samtools faidx ${reference}; fi
awk -v OFS='\t' '{print $1, 0, $2, $1, 0, "+"}' ${reference}.fai > regions.bed

eligos2 pair_diff_mod \
  -tbam ${native_bam} -cbam ${ivt_bam} \
  -reg regions.bed -ref ${reference} -t ${task.cpus} \
  --min_depth ${params.eligos_min_depth} \
  --max_depth ${params.eligos_max_depth} \
  --pval ${params.eligos_pval_thr} \
  --oddR ${params.eligos_oddR_thr} \
  --esb ${params.eligos_esb_thr} \
  -o ${prefix}
```
Defaults: `min_depth=1`, `max_depth=10000`, `pval=1.0`, `oddR=0`, `esb=0`.
- Version reporting line with provenance: executed runtime `2.1.0`; module source `pip: eligos2` with container `maestsi/nf-eligos:latest`; sources: runtime manifest + `modules/local/eligos_pair_diff_mod/main.nf` + `modules/local/eligos_pair_diff_mod/environment.yml` + `conf/modules.config`; note: container tag is floating (`latest`).
- Caveat: `errorStrategy='ignore'`.

### 26. EPINANO_ERROR
- Tool/process name: `EPINANO_ERROR`
- Methodological role: mismatch-based and sum-error-based differential comparison between native and IVT.
- Inputs/outputs:
  - Biological input: paired native/IVT BAMs and reference.
  - File input: tuple `[meta, native_bam, native_bai, ivt_bam, ivt_bai, reference.fa]`.
  - File output: `${prefix}/` containing per-site and differential CSV outputs.
- Exact command pattern and key defaults:
```bash
# If epinano_home not provided, clone pinned source
# git checkout eba4700953cc6e6e0ad0a4f846e7e071c43fe51c
python Epinano_Variants.py -r ${reference} -b ${native_bam} -c ${task.cpus} -o native_variants
python Epinano_Variants.py -r ${reference} -b ${ivt_bam} -c ${task.cpus} -o ivt_variants

Rscript Epinano_DiffErr.R \
  -k ivt_fixed.csv -w native_fixed.csv \
  -t ${params.epinano_zscore_threshold} \
  -o ${prefix}/${meta.id}_mismatch \
  -c ${params.epinano_coverage_threshold} \
  -f mis -d ${params.epinano_error_threshold}

python Epinano_sumErr.py --file ${native_var} --out native_sum_err.csv --kmer 0
python Epinano_sumErr.py --file ${ivt_var} --out ivt_sum_err.csv --kmer 0

Rscript Epinano_DiffErr.R \
  -k ivt_fixed.csv -w native_fixed.csv \
  -t ${params.epinano_zscore_threshold} \
  -o ${prefix}/${meta.id}_sumerr \
  -c ${params.epinano_coverage_threshold} \
  -f sum_err -d ${params.epinano_error_threshold}
```
Defaults: `zscore_threshold=0`, `coverage_threshold=0`, `error_threshold=0.1`.
- Version reporting line with provenance: executed runtime `epinano unknown; python 3.12.3; R 4.3.3`; pinned source commit `eba4700953cc6e6e0ad0a4f846e7e071c43fe51c` (cloned at runtime if not preinstalled); sources: runtime manifest + `modules/local/epinano_error/main.nf` + `modules/local/epinano_error/environment.yml` + `conf/modules.config`.
- Caveat: `errorStrategy='ignore'`.

### 27. DIFFERR
- Tool/process name: `DIFFERR`
- Methodological role: differential error-rate testing on paired BAMs.
- Inputs/outputs:
  - Biological input: paired native/IVT BAMs and reference.
  - File input: tuple `[meta, native_bam, native_bai, ivt_bam, ivt_bai, reference.fa]`.
  - File output: `${prefix}.bed` (plus optional `${prefix}.hdf5`).
- Exact command pattern and key defaults:
```bash
differr \
  -b ${native_bam} \
  -a ${ivt_bam} \
  -r ${reference} \
  -o ${prefix}.bed \
  -f ${params.differr_fdr_threshold} \
  --median-expr-threshold 0 --min-expr-threshold 0
```
Default `differr_fdr_threshold=1.0`.
- Version reporting line with provenance: executed runtime `differr unknown; python 3.11.9`; pinned source `git+https://github.com/bartongroup/differr_nanopore_DRS.git@a554d00c894d76cd4e91e2312cf88d611b0af3d2`; sources: runtime manifest + `modules/local/differr/main.nf` + `modules/local/differr/environment.yml` + `conf/modules.config`.
- Caveat: `errorStrategy='ignore'`.

### 28. DRUMMER
- Tool/process name: `DRUMMER`
- Methodological role: paired differential signal testing with odds-ratio and p-value outputs.
- Inputs/outputs:
  - Biological input: paired native/IVT BAMs and reference.
  - File input: tuple `[meta, native_bam, native_bai, ivt_bam, ivt_bai, reference.fa]`.
  - File output: `${prefix}/summary.txt` and related files.
- Exact command pattern and key defaults:
```bash
# If drummer_home not provided, clone pinned source
# git checkout 6683822c6210083e4ab0eecb4b6327e3c55f4c46

python DRUMMER.py \
  -r ${reference} \
  -n ${ref_name} \
  -c ${native_bam} \
  -t ${ivt_bam} \
  -o ${prefix} \
  -a exome \
  -p ${params.drummer_pval_threshold} \
  -z ${effective_odds} \
  -f 0
```
Defaults: `drummer_pval_threshold=1.0`, requested odds ratio `0.0` (sanitized to epsilon if non-positive).
- Version reporting line with provenance: executed runtime `drummer unknown; python 3.10.19`; pinned source commit `6683822c6210083e4ab0eecb4b6327e3c55f4c46`; sources: runtime manifest + `modules/local/drummer/main.nf` + `modules/local/drummer/environment.yml` + `conf/modules.config`.
- Caveat: `errorStrategy='ignore'`.

### 29. JACUSA2
- Tool/process name: `JACUSA2`
- Methodological role: pairwise all-sites comparative calling from BAMs with MD-tag-aware inputs.
- Inputs/outputs:
  - Biological input: paired native/IVT BAMs and reference.
  - File input: tuple `[meta, native_bam, native_bai, ivt_bam, ivt_bai, reference.fa]`.
  - File output: `${prefix}.bed`.
- Exact command pattern and key defaults:
```bash
samtools calmd -b ${native_bam} ${reference} > native_md.bam
samtools index native_md.bam
samtools calmd -b ${ivt_bam} ${reference} > ivt_md.bam
samtools index ivt_md.bam

java -jar ${JACUSA_JAR} call-2 \
  -R ${reference} \
  -p ${task.cpus} \
  -T 0 \
  -A \
  -r ${prefix}.bed \
  -c 1 -m 0 -q 0 \
  native_md.bam ivt_md.bam
```
- Version reporting line with provenance: executed runtime `jacusa2 2.0.4; java 1.8.0_472`; module pin `bioconda::jacusa2=2.0.4`, `samtools=1.17`; sources: runtime manifest + `modules/local/jacusa2/main.nf` + `modules/local/jacusa2/environment.yml` + `conf/modules.config`.
- Caveat: `errorStrategy='ignore'`.

### 30. DOWNSTREAM_ANALYSIS
- Tool/process name: `DOWNSTREAM_ANALYSIS`
- Methodological role: parse tool outputs, harmonize metrics, compare tools/replicates, and produce collation artifacts.
- Inputs/outputs:
  - Biological input: caller outputs under `${outdir}/modifications` plus optional ground-truth positions.
  - File input: `modifications_dir`, `ground_truth`, `references_csv`.
  - File output: `downstream_analysis/` (metadata, by_reference, collation tables, plots, logs).
- Exact command pattern and key defaults:
```bash
python bin/downstream_analysis/run_analysis.py \
  --input-dir ${modifications_dir} \
  --output-dir downstream_analysis \
  --min-replicates 2 \
  --expected-replicates 3 \
  [--threshold <ext.threshold>] \
  [--ground-truth ${ground_truth}] \
  [--references-csv ${references_csv}] \
  [--run-id ${params.downstream_run_id}] \
  [--coverage-label ${params.downstream_coverage_label}] \
  [--quality-label ${params.downstream_quality_label}]
```
- Version reporting line with provenance: executed runtime `python 3.14.3; pandas 3.0.0; sklearn 1.8.0; matplotlib 3.10.8`; module environment constraints `python>=3.10`, `pandas>=1.5.0`, `numpy>=1.24.0`, `scikit-learn>=1.3.0`, `matplotlib>=3.7.0`, `seaborn>=0.12.0`; sources: runtime manifest + `modules/local/downstream_analysis/main.nf` + `modules/local/downstream_analysis/environment.yml`.

## Version Provenance Table

| Process | Tool | Runtime version | Pinned/source version | Version source file(s) | Notes |
|---|---|---|---|---|---|
| MINIMAP2_ALIGN | minimap2 | 2.24-r1122 | bioconda minimap2=2.24; container 2.24--h7132678_1 | `results_rasusa/results_100x/pipeline_info/software_versions.yml`; `modules/local/minimap2_align/main.nf` | Runtime and pin aligned at major/minor level. |
| SAMTOOLS_VIEW | samtools | 1.17 | bioconda samtools=1.17; container 1.17--h00cdaf9_0 | runtime manifest; `modules/local/samtools_view/main.nf` | Aligned. |
| SAMTOOLS_SORT | samtools | 1.17 | bioconda samtools=1.17 | runtime manifest; `modules/local/samtools_sort/main.nf` | Aligned. |
| SAMTOOLS_INDEX | samtools | 1.17 | bioconda samtools=1.17 | runtime manifest; `modules/local/samtools_index/main.nf` | Aligned. |
| EXTRACT_MAPPED_READS | samtools | 1.17 | bioconda samtools=1.17 | runtime manifest; `modules/local/extract_mapped_reads/main.nf` | Aligned. |
| SAMTOOLS_FLAGSTAT | samtools | 1.17 | bioconda samtools=1.17 | runtime manifest; `modules/local/samtools_flagstat/main.nf` | Aligned. |
| SAMTOOLS_DEPTH | samtools | 1.17 | bioconda samtools=1.17 | runtime manifest; `modules/local/samtools_depth/main.nf` | Aligned. |
| COVERAGE_PLOT | python/pandas/matplotlib/seaborn | python 3.9.23; pandas 1.3.4; matplotlib 3.4.3; seaborn 0.11.2 | env pins python=3.9; pandas=1.3.4; matplotlib=3.4.3; seaborn=0.11.2 | runtime manifest; `modules/local/coverage_plot/main.nf`; `modules/local/coverage_plot/environment.yml` | Fully aligned. |
| NANOPLOT_BAM | NanoPlot | 1.41.0 | bioconda nanoplot=1.41.0; container 1.41.0--pyhdfd78af_0 | runtime manifest; `modules/local/nanoplot_bam/main.nf` | Aligned. |
| EXTRACT_READ_IDS | seqkit | seqkit v2.3.0 | bioconda seqkit=2.3.1; container 2.3.1--h9ee0642_0 | runtime manifest; `modules/local/extract_read_ids/main.nf` | Runtime older than pin. |
| FAST5_SUBSET | ont-fast5-api | unknown | bioconda ont-fast5-api=4.1.0; container 4.1.0--pyhdfd78af_0 | runtime manifest; `modules/local/fast5_subset/main.nf` | Runtime probe unresolved. |
| MULTI_TO_SINGLE_FAST5 | ont-fast5-api | unknown | bioconda ont-fast5-api=4.1.0; container 4.1.0--pyhdfd78af_0 | runtime manifest; `modules/local/multi_to_single_fast5/main.nf` | Runtime probe unresolved. |
| F5C_INDEX | f5c | unknown | bioconda f5c=1.5; container 1.5--hee927d3_2 | runtime manifest; `modules/local/f5c_index/main.nf` | Runtime probe unresolved; process configured `errorStrategy=ignore`. |
| F5C_EVENTALIGN | f5c | unknown | bioconda f5c=1.5; container 1.5--hee927d3_2 | runtime manifest; `modules/local/f5c_eventalign/main.nf`; `conf/modules.config` | Runtime probe unresolved; `errorStrategy=ignore`. |
| F5C_EVENTALIGN_XPORE | f5c | unknown | bioconda f5c=1.5; container 1.5--hee927d3_2 | runtime manifest; `modules/local/f5c_eventalign_xpore/main.nf` | Runtime probe unresolved; `errorStrategy=ignore`. |
| TOMBO_RESQUIGGLE | tombo | 1.5.1 | ont-tombo=1.5.1 (python=3.7; numpy<1.20) | runtime manifest; `modules/local/tombo_resquiggle/main.nf`; `modules/local/tombo_resquiggle/environment.yml` | Aligned; `errorStrategy=ignore`. |
| TOMBO_DETECT_MODIFICATIONS | tombo | 1.5.1 | ont-tombo=1.5.1 | runtime manifest; `modules/local/tombo_detect_modifications/main.nf` | Aligned; `errorStrategy=ignore`. |
| TOMBO_TEXT_OUTPUT | tombo/pandas | tombo 1.5.1; pandas 1.3.4 | ont-tombo=1.5.1; pandas=1.3.4 | runtime manifest; `modules/local/tombo_text_output/main.nf`; `modules/local/tombo_text_output/environment.yml` | Aligned; `errorStrategy=ignore`. |
| YANOCOMP_PREPARE | yanocomp | unknown | pip install git+github.com/bartongroup/yanocomp.git@v0.2; python=3.8 | runtime manifest; `modules/local/yanocomp_prepare/main.nf`; `modules/local/yanocomp_prepare/environment.yml`; Dockerfile | Runtime probe unresolved. |
| YANOCOMP_ANALYSIS | yanocomp | unknown | yanocomp source v0.2 | runtime manifest; `modules/local/yanocomp_analysis/main.nf`; `modules/local/yanocomp_analysis/environment.yml` | Runtime probe unresolved; `errorStrategy=ignore`; one ignored task in representative run. |
| NANOCOMPORE_EVENTALIGN_COLLAPSE | nanocompore | unknown | nanocompore=1.0.4; container 1.0.4--pyhdfd78af_0 | runtime manifest; `modules/local/nanocompore_eventalign_collapse/main.nf`; environment.yml | Runtime probe unresolved. |
| NANOCOMPORE_SAMPCOMP | nanocompore | unknown | nanocompore=1.0.4; container 1.0.4--pyhdfd78af_0 | runtime manifest; `modules/local/nanocompore_sampcomp/main.nf`; environment.yml | Runtime probe unresolved; `errorStrategy=ignore`. |
| XPORE_DATAPREP | xpore | 2.0 | xpore=2.1; container 2.1--pyh5e36f6f_0 | runtime manifest; `modules/local/xpore_dataprep/main.nf`; environment.yml | Runtime lower than pin. |
| XPORE_DIFFMOD | xpore | 2.0 | xpore=2.1; container 2.1--pyh5e36f6f_0 | runtime manifest; `modules/local/xpore_diffmod/main.nf`; environment.yml | Runtime lower than pin; `errorStrategy=ignore`. |
| ELIGOS_PAIR_DIFF_MOD | eligos2 | 2.1.0 | env pip eligos2; container `maestsi/nf-eligos:latest` | runtime manifest; `modules/local/eligos_pair_diff_mod/main.nf`; environment.yml | Container uses floating `latest`; `errorStrategy=ignore`. |
| EPINANO_ERROR | EpiNano | unknown; python 3.12.3; R 4.3.3 | source cloned at commit eba4700953cc6e6e0ad0a4f846e7e071c43fe51c; env python>=3.8, samtools=1.17, R>=4.0 | runtime manifest; `modules/local/epinano_error/main.nf`; environment.yml | Runtime package version unresolved; source commit pinned; `errorStrategy=ignore`. |
| DIFFERR | differr | unknown; python 3.11.9 | pip git+https://github.com/bartongroup/differr_nanopore_DRS.git@a554d00c894d76cd4e91e2312cf88d611b0af3d2 | runtime manifest; `modules/local/differr/main.nf`; environment.yml | Runtime CLI version unresolved; source commit pinned; `errorStrategy=ignore`. |
| DRUMMER | DRUMMER | unknown; python 3.10.19 | source cloned at commit 6683822c6210083e4ab0eecb4b6327e3c55f4c46 | runtime manifest; `modules/local/drummer/main.nf`; environment.yml | Runtime package version unresolved; source commit pinned; `errorStrategy=ignore`. |
| JACUSA2 | JACUSA2 | jacusa2 2.0.4; java 1.8.0_472 | bioconda jacusa2=2.0.4; samtools=1.17 | runtime manifest; `modules/local/jacusa2/main.nf`; environment.yml | Aligned; `errorStrategy=ignore`. |
| DOWNSTREAM_ANALYSIS | python/pandas/sklearn/matplotlib | python 3.14.3; pandas 3.0.0; sklearn 1.8.0; matplotlib 3.10.8 | env lower bounds: python>=3.10; pandas>=1.5.0; scikit-learn>=1.3.0; matplotlib>=3.7.0 | runtime manifest; `modules/local/downstream_analysis/main.nf`; environment.yml | Runtime exceeds lower-bound pins. |

## Validation and QA Checks (Applied)

1. Completeness check
- All 30 active workflow processes listed in `workflows/rnamodbench.nf` and subworkflows were documented once in the supplement.

2. Command fidelity check
- Command patterns and defaults were extracted from each module `main.nf` and reconciled with `conf/modules.config` and `nextflow.config` overrides.

3. Version traceability check
- Every version line cites concrete source files (runtime manifest and module pin source).

4. Consistency check
- Software manifest consistency confirmed by shared checksum `e25779bda6e4a901f13b83abbd81a8a28d9af71e24a84be196731a3c4b8ee1d5` across primary run directories.

5. Reproducibility check
- Runtime profile, Nextflow version/build, Wave/Fusion status, and representative command are reported from `results_rasusa/results_100x/pipeline_info/execution_report_2026-02-17_13-26-10.html`.

6. Caveat check (error strategy)
- Modules explicitly configured with `errorStrategy='ignore'` and relevant to this scope: `TOMBO_RESQUIGGLE`, `TOMBO_DETECT_MODIFICATIONS`, `TOMBO_TEXT_OUTPUT`, `F5C_INDEX`, `F5C_EVENTALIGN`, `F5C_EVENTALIGN_XPORE`, `YANOCOMP_ANALYSIS`, `NANOCOMPORE_SAMPCOMP`, `XPORE_DIFFMOD`, `ELIGOS_PAIR_DIFF_MOD`, `EPINANO_ERROR`, `DIFFERR`, `DRUMMER`, `JACUSA2`.

