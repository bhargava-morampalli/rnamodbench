process YANOCOMP_PREPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::yanocomp=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/home/bmorampa/containers/yanocomp.sif' :
        'yanocomp:latest' }"

    input:
    tuple val(meta), path(eventalign)

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    yanocomp prep \\
        -p ${task.cpus} \\
        -e $eventalign \\
        -h ${prefix}.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yanocomp: \$(yanocomp --version 2>&1 | sed 's/^.*yanocomp //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hdf5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yanocomp: \$(yanocomp --version 2>&1 | sed 's/^.*yanocomp //; s/ .*\$//')
    END_VERSIONS
    """
}
