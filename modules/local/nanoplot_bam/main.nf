process NANOPLOT_BAM {
    tag "$meta.id"
    label 'process_low'  // rnamodbench standard: nanoplot uses process_low

    conda "bioconda::nanoplot=1.41.0"
    container "quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.feather"), emit: stats
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_feather.py --bam $bam --output ${prefix}.feather

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed 's/NanoPlot //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.feather

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed 's/NanoPlot //')
    END_VERSIONS
    """
}
