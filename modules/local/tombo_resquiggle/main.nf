process TOMBO_RESQUIGGLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // No pre-built container available for ont-tombo - Wave will build from conda environment
    container null

    input:
    tuple val(meta), path(fast5), path(reference)

    output:
    tuple val(meta), path(fast5), emit: resquiggled
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--rna --overwrite --num-most-common-errors 5'
    """
    tombo resquiggle \\
        $args \\
        --processes $task.cpus \\
        $fast5 \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tombo: \$(tombo --version 2>&1 | sed 's/^.*tombo //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${fast5}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tombo: 1.5.1
    END_VERSIONS
    """
}
