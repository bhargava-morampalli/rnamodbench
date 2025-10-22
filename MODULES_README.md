# Module Implementation Guide

This document provides a template and guidelines for implementing the remaining modules following nf-core standards.

## Module Structure Template

Each module should follow this structure:

```groovy
process MODULE_NAME {
    tag "$meta.id"
    label 'process_[single|low|medium|high]'

    conda "bioconda::tool=version"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tool:version' :
        'biocontainers/tool:version' }"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.output"), emit: output
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tool_command \\
        $args \\
        -t $task.cpus \\
        $input_file \\
        > ${prefix}.output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version 2>&1)
    END_VERSIONS
    """
}
```

## Modules to Implement

Below is a list of modules that need to be created based on your current pipeline:

### Completed Modules ✓
- [x] SAMPLESHEET_CHECK
- [x] MINIMAP2_ALIGN
- [x] SAMTOOLS_VIEW
- [x] TOMBO_RESQUIGGLE
- [x] F5C_EVENTALIGN
- [x] YANOCOMP_ANALYSIS

### Modules to Create

#### SAM/BAM Processing
1. **SAMTOOLS_SORT** - Sort BAM files
2. **SAMTOOLS_INDEX** - Index BAM files
3. **SAMTOOLS_FLAGSTAT** - Generate mapping statistics
4. **SAMTOOLS_DEPTH** - Calculate depth per position

#### Read Extraction
5. **EXTRACT_MAPPED_READS** - Extract mapped reads from SAM to FASTQ
6. **EXTRACT_READ_IDS** - Extract read IDs from FASTQ files

#### FAST5 Processing
7. **FAST5_SUBSET** - Subset FAST5 files based on read IDs
8. **MULTI_TO_SINGLE_FAST5** - Convert multi-read to single-read FAST5

#### F5C Processing
9. **F5C_INDEX** - Index FASTQ for f5c

#### Tombo Processing
10. **TOMBO_DETECT_MODIFICATIONS** - Detect modifications comparing native vs IVT
11. **TOMBO_TEXT_OUTPUT** - Extract text output from Tombo stats

#### Yanocomp Processing
12. **YANOCOMP_PREPARE** - Prepare HDF5 files for yanocomp

#### Xpore Processing
13. **XPORE_DATAPREP** - Prepare data for xpore
14. **XPORE_DIFFMOD** - Differential modification analysis

#### QC
15. **NANOPLOT_BAM** - Generate read quality plots
16. **COVERAGE_PLOT** - Create coverage plots

## Example Module Implementations

### SAMTOOLS_SORT

```groovy
process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.sorted.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
```

### SAMTOOLS_INDEX

```groovy
process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam_bai
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools index -@ $task.cpus $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
```

### EXTRACT_MAPPED_READS

```groovy
process EXTRACT_MAPPED_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.fastq"), emit: fastq
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -F 4 $sam | \\
        samtools fastq -@ $task.cpus - > ${prefix}_mapped.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mapped.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
```

### EXTRACT_READ_IDS

```groovy
process EXTRACT_READ_IDS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.txt"), emit: read_ids
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'NR%4==1 {print substr(\$1,2)}' $fastq > ${prefix}_read_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_read_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
```

## Best Practices

1. **Always include a `meta` map** in inputs and outputs for sample tracking
2. **Use `task.ext.args`** for customizable parameters
3. **Include version tracking** in every module
4. **Add resource labels** (process_single, process_low, process_medium, process_high)
5. **Implement stub sections** for pipeline testing
6. **Use `when` directive** for conditional execution
7. **Follow nf-core naming conventions**: TOOL_SUBTOOL format
8. **Document container sources** for both Singularity and Docker

## Container Management

For custom tools not in biocontainers:
- Store Singularity images in a centralized location
- Use absolute paths in the module definition
- Document how to build/obtain the container

Example:
```groovy
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    '/path/to/custom/tool.sif' :
    'dockerhub/custom-tool:latest' }"
```

## Testing

Each module should be tested with:
1. Small test data
2. Expected outputs verified
3. Version tracking confirmed
4. Resource usage monitored
