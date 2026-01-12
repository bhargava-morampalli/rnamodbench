process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::minimap2=2.24"
    container "quay.io/biocontainers/minimap2:2.24--h7132678_1"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Alignment mode per benchmarking papers:
    // - For GENOME alignment (with splicing): use '-ax splice -uf -k14 --secondary=no' (default)
    // - For TRANSCRIPTOME alignment: override with ext.args = '-ax map-ont --secondary=no'
    def args = task.ext.args ?: '-ax splice -uf -k14 --secondary=no'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $reference \\
        $reads \\
        > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
