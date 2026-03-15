# Usage

## Running the pipeline

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --references references.csv \
    --outdir results \
    -profile singularity
```

Use `-resume` to restart from cached results after an interruption. Use `-bg`, `screen`, or `tmux` to run in the background.

## Samplesheet (`--input`)

CSV with the following columns:

| Column | Description |
| --- | --- |
| `sample` | Unique sample ID (no spaces) |
| `fastq` | Path to FASTQ file (`.fq`, `.fastq`, optionally `.gz`) |
| `type` | `native` or `ivt` |
| `replicate` | Replicate identifier (`rep1`, `rep2`, ...) |
| `fast5_dir` | Directory containing FAST5 files |
| `target` | Target label (e.g., `16s`, `23s`) |

Each `target + replicate` must have both a `native` and `ivt` sample.

## References (`--references`)

CSV mapping targets to reference FASTA files:

```csv
target,reference
16s,/refs/k12_16S.fa
23s,/refs/k12_23S.fa
```

## Parameters

### Input/output

| Parameter | Description | Default |
| --- | --- | --- |
| `--input` | Samplesheet CSV | Required |
| `--references` | Target-to-reference CSV | Required |
| `--outdir` | Output directory | `./results` |

### Mapping

| Parameter | Description | Default |
| --- | --- | --- |
| `--minimap2_args` | Additional minimap2 arguments | `-ax splice -uf -k14 --secondary=no` |

### Modification calling thresholds

| Parameter | Description | Default |
| --- | --- | --- |
| `--yanocomp_fdr_threshold` | Yanocomp FDR threshold | `1.0` |
| `--yanocomp_min_ks` | Yanocomp minimum KS statistic | `0.0` |
| `--nanocompore_min_coverage` | Nanocompore minimum coverage | `1` |
| `--eligos_min_depth` | ELIGOS minimum read depth | `1` |
| `--eligos_pval_thr` | ELIGOS p-value threshold | `1.0` |
| `--epinano_error_threshold` | EpiNano error-rate threshold | `0.1` |
| `--differr_fdr_threshold` | DiffErr FDR threshold | `1.0` |
| `--drummer_pval_threshold` | DRUMMER p-value threshold | `1.0` |

> Defaults are set permissively (no filtering) so all tested positions are reported for downstream ROC/PR analysis.

### Downstream analysis

| Parameter | Description | Default |
| --- | --- | --- |
| `--run_downstream` | Enable downstream analysis | `false` |
| `--ground_truth` | Ground truth CSV for benchmarking | None |
| `--downstream_differr_score_field` | DiffErr score: `both`, `g_fdr_neglog10`, `g_stat` | `both` |
| `--downstream_run_id` | Run identifier | None |
| `--downstream_coverage_label` | Coverage label | None |
| `--downstream_quality_label` | Quality label | None |

### Resource limits

| Parameter | Description | Default |
| --- | --- | --- |
| `--max_cpus` | Maximum CPUs | Auto-detected |
| `--max_memory` | Maximum memory | Auto-detected (90%) |
| `--max_time` | Maximum wall time | `240.h` |

## Profiles

| Profile | Description |
| --- | --- |
| `singularity` | Singularity containers |
| `docker` | Docker containers |
| `conda` / `mamba` | Conda environments |
| `apptainer` | Apptainer containers |
| `test` | Minimal test dataset |

Multiple profiles can be combined: `-profile test,docker`

## Custom configuration

Use `-c custom.config` to override resource requests or process settings. See `conf/base.config` for defaults.
