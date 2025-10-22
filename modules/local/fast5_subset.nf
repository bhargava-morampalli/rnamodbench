process FAST5_SUBSET {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ont-fast5-api=4.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-fast5-api:4.1.0--pyhdfd78af_0' :
        'biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fast5_dir), path(readids)

    output:
    tuple val(meta), path("${prefix}_fast5s"), emit: fast5
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fast5s
    fast5_subset -i $fast5_dir -s ${prefix}_fast5s -l $readids --recursive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ont-fast5-api: \$(echo \$(fast5_subset --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_fast5s

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ont-fast5-api: \$(echo \$(fast5_subset --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
