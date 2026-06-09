process SAMTOOLS_VIEW {
    tag "$meta.id"
    label 'process_low'  // rnamodbench standard: samtools view is lightweight

    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval('(samtools --version 2>/dev/null || echo unknown) | head -n 1 | cut -d " " -f 2'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-b -F 4'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        $sam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    """
}
