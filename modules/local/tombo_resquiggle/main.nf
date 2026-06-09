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
    tuple val("${task.process}"), val('tombo'), eval('tombo --version 2>&1 | grep -oP \'[0-9]+\\.[0-9]+[0-9.]*\' | head -1 || echo unknown'), emit: versions_tombo, topic: versions

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
    """

    stub:
    """
    mkdir -p ${fast5}
    """
}
