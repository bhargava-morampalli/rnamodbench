process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.txt"), emit: depth
    tuple val("${task.process}"), val('samtools'), eval('(samtools --version 2>/dev/null || echo unknown) | head -n 1 | cut -d " " -f 2'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools depth -a -m 0 $bam > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
