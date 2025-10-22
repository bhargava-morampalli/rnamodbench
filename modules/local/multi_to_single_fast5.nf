process MULTI_TO_SINGLE_FAST5 {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ont-fast5-api=4.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-fast5-api:4.1.0--pyhdfd78af_0' :
        'biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fast5_dir)

    output:
    tuple val(meta), path("${prefix}_single"), emit: single_fast5
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_single
    multi_to_single_fast5 --input_path $fast5_dir --save_path ${prefix}_single --recursive

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ont-fast5-api: \$(echo \$(multi_to_single_fast5 --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_single

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ont-fast5-api: \$(echo \$(multi_to_single_fast5 --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
