process MULTI_TO_SINGLE_FAST5 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ont-fast5-api=4.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-fast5-api:4.1.0--pyhdfd78af_0' :
        'biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fast5_dir)

    output:
    tuple val(meta), path("single_fast5"), emit: fast5
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multi_to_single_fast5 \\
        $args \\
        --input_path ${fast5_dir} \\
        --save_path single_fast5 \\
        --recursive

    echo -e '"${task.process}":\n  ont-fast5-api: \$(multi_to_single_fast5 --version | sed -nE "s/multi_to_single_fast5, version (.*)/\\1/p")' > versions.yml
    """
}