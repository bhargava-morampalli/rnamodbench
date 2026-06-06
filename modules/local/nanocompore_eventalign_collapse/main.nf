process NANOCOMPORE_EVENTALIGN_COLLAPSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/nanocompore:1.0.4--pyhdfd78af_0"

    input:
    tuple val(meta), path(eventalign)

    output:
    tuple val(meta), path("${prefix}"), emit: collapsed
    tuple val("${task.process}"), val('nanocompore'), eval('nanocompore --version 2>/dev/null | grep -oP \'[0-9]+\\.[0-9]+[0-9.]*\' | head -1 || echo unknown'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    nanocompore eventalign_collapse \\
        -i $eventalign \\
        -o ${prefix} \\
        -t ${task.cpus}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/eventalign.index
    touch ${prefix}/out_eventalign_collapse.tsv
    touch ${prefix}/out_eventalign_collapse.idx
    """
}
