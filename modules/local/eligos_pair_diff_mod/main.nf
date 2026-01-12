process ELIGOS_PAIR_DIFF_MOD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://maestsi/nf-eligos:latest' :
        'maestsi/nf-eligos:latest' }"

    input:
    tuple val(meta), path(native_bam), path(native_bai), path(ivt_bam), path(ivt_bai), path(reference)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "*.log"                      , emit: log, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_eligos"
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}.log) 2>&1

    echo "=== ELIGOS_PAIR_DIFF_MOD started at \$(date) ==="
    echo "Sample: ${meta.id}"
    echo "Native BAM: $native_bam"
    echo "IVT BAM: $ivt_bam"
    echo "Reference: $reference"
    echo "Arguments: $args"
    echo "CPUs: $task.cpus"
    echo "==========================================="

    # Index reference if not already indexed
    if [ ! -f ${reference}.fai ]; then
        samtools faidx $reference
    fi

    # Generate BED file from reference FASTA index (covers entire reference)
    awk -v OFS='\\t' '{print \$1, 0, \$2, \$1, 0, "+"}' ${reference}.fai > regions.bed

    # Run ELIGOS pair_diff_mod
    eligos2 pair_diff_mod \\
        -tbam $native_bam \\
        -cbam $ivt_bam \\
        -reg regions.bed \\
        -ref $reference \\
        -t $task.cpus \\
        $args \\
        -o $prefix

    echo "=== ELIGOS_PAIR_DIFF_MOD completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eligos2: \$(eligos2 --version 2>&1 | sed 's/^.*eligos2 //; s/ .*\$//' || echo "2.1.0")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_eligos"
    """
    mkdir -p $prefix
    touch ${prefix}/test_paired_diff_mod_result.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eligos2: 2.1.0
    END_VERSIONS
    """
}
