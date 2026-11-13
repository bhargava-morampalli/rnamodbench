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
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed \'s/Python //\' || echo unknown'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c \'import pandas; print(pandas.__version__)\' 2>/dev/null || echo unknown'), emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval('python -c \'import matplotlib; print(matplotlib.__version__)\' 2>/dev/null || echo unknown'), emit: versions_matplotlib, topic: versions
    tuple val("${task.process}"), val('seaborn'), eval('python -c \'import seaborn; print(seaborn.__version__)\' 2>/dev/null || echo unknown'), emit: versions_seaborn, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverage_plot.py -f $depth -t ${meta.id} -o ${prefix}.pdf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf
    """
}
