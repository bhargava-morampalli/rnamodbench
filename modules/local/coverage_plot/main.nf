process COVERAGE_PLOT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // No pre-built container available - Wave will build from conda environment
    container null

    input:
    tuple val(meta), path(depth)

    output:
    tuple val(meta), path("*.pdf"), emit: plot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverage_plot.py -f $depth -t ${meta.id} -o ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
        seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
