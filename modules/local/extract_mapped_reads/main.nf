process EXTRACT_MAPPED_READS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.fastq"), emit: fastq
    tuple val("${task.process}"), val('samtools'), eval('echo $(samtools --version 2>&1) | sed \'s/^.*samtools //; s/Using.*$//\' || echo unknown'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        fastq \\
        -F 4 \\
        --threads ${task.cpus-1} \\
        $args \\
        $sam \\
        > ${prefix}.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastq
    """
}
