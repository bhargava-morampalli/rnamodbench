process YANOCOMP_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // No pre-built container available for yanocomp - Wave will build from conda environment
    container null

    input:
    tuple val(meta), path(hdf5_native), path(hdf5_ivt)

    output:
    tuple val(meta), path("*.bed")             , emit: bed
    tuple val(meta), path("*_sm_preds.json")   , emit: json
    path "*.log"                               , emit: log, optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--fdr-threshold 1 --min-ks 0'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}.log) 2>&1

    echo "=== YANOCOMP_ANALYSIS started at \$(date) ==="
    echo "Sample: ${meta.id}"
    echo "Native HDF5: $hdf5_native"
    echo "IVT HDF5: $hdf5_ivt"
    echo "Arguments: $args"
    echo "CPUs: $task.cpus"
    echo "==========================================="

    yanocomp gmmtest \\
        $args \\
        -p $task.cpus \\
        -n 1 \\
        -c $hdf5_native \\
        -t $hdf5_ivt \\
        -o ${prefix}.bed \\
        -s ${prefix}_sm_preds.json

    echo "=== YANOCOMP_ANALYSIS completed at \$(date) ==="

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
        yanocomp: 1.0.0
    END_VERSIONS
    """
}
