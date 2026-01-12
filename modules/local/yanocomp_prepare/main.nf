process YANOCOMP_PREPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // No pre-built container available for yanocomp - Wave will build from conda environment
    container null

    input:
    tuple val(meta), path(eventalign)
    path summary  // Optional eventalign summary file per benchmarking papers

    output:
    tuple val(meta), path("*.hdf5"), emit: hdf5
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Add -s flag if summary file is provided (per benchmarking papers)
    def summary_arg = summary.name != 'NO_FILE' ? "-s ${summary}" : ""
    """
    yanocomp prep \\
        -p ${task.cpus} \\
        -e $eventalign \\
        ${summary_arg} \\
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
        yanocomp: 1.0.0
    END_VERSIONS
    """
}
