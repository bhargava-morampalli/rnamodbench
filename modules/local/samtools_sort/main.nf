process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*_sorted.bam"), emit: bam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -S -b -h $sam | samtools sort -o ${prefix}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
