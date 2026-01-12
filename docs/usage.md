# nf-core/rnamodifications: Usage

## Introduction

nf-core/rnamodifications is a bioinformatics pipeline for detecting RNA modifications from Oxford Nanopore direct RNA sequencing data. It compares native RNA samples against in-vitro transcribed (IVT) controls to identify modified positions using multiple detection algorithms.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/rnamodifications \
    --input samplesheet.csv \
    --outdir results \
    --ref_16s /path/to/16s_reference.fa \
    --ref_23s /path/to/23s_reference.fa \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```
work                # Directory containing the nextflow working files
results             # Finished results (configurable, see below)
.nextflow_log       # Log file from Nextflow
```

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

### Samplesheet format

The samplesheet must be a comma-separated file with the following columns:

| Column      | Description                                                  |
| ----------- | ------------------------------------------------------------ |
| `sample`    | Custom sample name. Must be unique and cannot contain spaces |
| `fastq`     | Path to FASTQ file                                           |
| `type`      | Sample type: either `native` or `ivt`                        |
| `replicate` | Replicate identifier (e.g., `rep1`, `rep2`)                  |
| `fast5_dir` | Path to directory containing FAST5 files                     |

An example samplesheet:

```csv
sample,fastq,type,replicate,fast5_dir
native_rep1,/data/native_rep1.fastq.gz,native,rep1,/data/fast5/native_rep1
ivt_rep1,/data/ivt_rep1.fastq.gz,ivt,rep1,/data/fast5/ivt_rep1
native_rep2,/data/native_rep2.fastq.gz,native,rep2,/data/fast5/native_rep2
ivt_rep2,/data/ivt_rep2.fastq.gz,ivt,rep2,/data/fast5/ivt_rep2
```

### Important notes

- Each replicate must have both a `native` and `ivt` sample for modification calling
- The pipeline will warn and skip replicates with missing native/ivt pairs
- FAST5 files are required for signal-level analysis (Tombo, F5C)

## Reference files

The pipeline requires two reference FASTA files:

```bash
--ref_16s /path/to/16s_reference.fa
--ref_23s /path/to/23s_reference.fa
```

These should contain the 16S and 23S rRNA sequences for your organism.

## Core Nextflow arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility.

The pipeline also dynamically loads configurations from [nf-core/configs](https://github.com/nf-core/configs) when it runs. See the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation) for more information.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`.

- `docker`: A generic configuration profile to be used with Docker
- `singularity`: A generic configuration profile to be used with Singularity
- `conda`: A generic configuration profile to be used with Conda
- `test`: A minimal test dataset configuration

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

### `-c`

Specify the path to a specific config file.

### Resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3ber3568c0eb9f9f/conf/base.config#L18), it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline will terminate.

To change the resource requests, please see the [max resources documentation](https://nf-co.re/docs/usage/configuration#max-resources).

## Pipeline parameters

### Input/output options

| Parameter | Description                        | Default     |
| --------- | ---------------------------------- | ----------- |
| `--input` | Path to samplesheet CSV file       | Required    |
| `--outdir`| Output directory for results       | `./results` |

### Reference options

| Parameter   | Description                   | Default  |
| ----------- | ----------------------------- | -------- |
| `--ref_16s` | Path to 16S rRNA reference    | Required |
| `--ref_23s` | Path to 23S rRNA reference    | Required |

### Mapping options

| Parameter         | Description              | Default                          |
| ----------------- | ------------------------ | -------------------------------- |
| `--minimap2_args` | Additional minimap2 args | `-ax splice -uf -k14 --secondary=no` |

### Signal processing options

| Parameter                | Description                 | Default                                        |
| ------------------------ | --------------------------- | ---------------------------------------------- |
| `--tombo_resquiggle_args`| Tombo resquiggle arguments  | `--rna --overwrite --num-most-common-errors 5` |
| `--f5c_eventalign_args`  | F5C eventalign arguments    | `--rna --scale-events --signal-index --print-read-names --samples` |

### Modification calling thresholds

| Parameter                    | Description                     | Default |
| ---------------------------- | ------------------------------- | ------- |
| `--yanocomp_fdr_threshold`   | Yanocomp FDR threshold          | `1.0`   |
| `--xpore_pvalue_threshold`   | Xpore p-value threshold         | `0.05`  |
| `--xpore_diffmod_threshold`  | Xpore diff modification cutoff  | `0.1`   |
| `--nanocompore_min_coverage` | Nanocompore min coverage        | `1`     |
| `--eligos_min_depth`         | ELIGOS minimum read depth       | `5`     |

### Resource limits

| Parameter     | Description       | Default   |
| ------------- | ----------------- | --------- |
| `--max_cpus`  | Maximum CPUs      | Auto-detected |
| `--max_memory`| Maximum memory    | Auto-detected (90% of system) |
| `--max_time`  | Maximum time      | `240.h`   |

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal. You can also use `screen` / `tmux` or similar tools to prevent the connection to your remote computer being lost.

Alternatively, you can use a workflow manager like [Tower](https://tower.nf/).
