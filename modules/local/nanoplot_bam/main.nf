process NANOPLOT_BAM {
    tag "$meta.id"
    label 'process_low'  // rnamodbench standard: nanoplot uses process_low

    conda "bioconda::nanoplot=1.41.0"
    container "quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.feather"), emit: stats
    tuple val("${task.process}"), val('nanoplot'), eval('(NanoPlot --version 2>/dev/null || echo unknown) | cut -d " " -f 2'), emit: versions_nanoplot, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_feather.py --bam $bam --output ${prefix}.feather
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.feather
    """
}
