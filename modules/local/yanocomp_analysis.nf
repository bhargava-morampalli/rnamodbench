process YANOCOMP_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::yanocomp=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/home/bmorampa/containers/yanocomp.sif' :
        'yanocomp:latest' }"

    input:
    tuple val(meta), path(hdf5_native), path(hdf5_ivt)

    output:
    tuple val(meta), path("*.bed")             , emit: bed
    tuple val(meta), path("*_sm_preds.json")   , emit: json
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--fdr-threshold 1 --min-ks 0'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    yanocomp gmmtest \\
        $args \\
        -p $task.cpus \\
        -n 1 \\
        -c $hdf5_native \\
        -t $hdf5_ivt \\
        -o ${prefix}.bed \\
        -s ${prefix}_sm_preds.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yanocomp: \$(yanocomp --version 2>&1 | sed 's/^.*yanocomp //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}_sm_preds.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yanocomp: \$(yanocomp --version 2>&1 | sed 's/^.*yanocomp //; s/ .*\$//')
    END_VERSIONS
    """
}
