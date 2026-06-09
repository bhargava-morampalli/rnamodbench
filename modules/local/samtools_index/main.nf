process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bai"), emit: bai
    tuple val("${task.process}"), val('samtools'), eval('(samtools --version 2>/dev/null || echo unknown) | head -n 1 | cut -d " " -f 2'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        $args \\
        $bam
    """

    stub:
    """
    touch ${bam}.bai
    """
}
