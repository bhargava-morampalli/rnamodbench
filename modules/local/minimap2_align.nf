process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::minimap2=2.24"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.24--h7132678_1' :
        'biocontainers/minimap2:2.24--h7132678_1' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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
