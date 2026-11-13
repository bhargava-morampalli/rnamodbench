process DOWNSTREAM_ANALYSIS {
    tag "downstream_analysis"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container null

    input:
    path(modifications_dir)
    val(ground_truth)
    val(references_csv)

    output:
    path("downstream_analysis")         , emit: results
    path "*.log"                        , emit: log, optional: true
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed \'s/Python //\' || echo unknown'), topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c \'import pandas; print(pandas.__version__)\' 2>/dev/null || echo unknown'), topic: versions
    tuple val("${task.process}"), val('sklearn'), eval('python -c \'import sklearn; print(sklearn.__version__)\' 2>/dev/null || echo unknown'), topic: versions
    tuple val("${task.process}"), val('matplotlib'), eval('python -c \'import matplotlib; print(matplotlib.__version__)\' 2>/dev/null || echo unknown'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gt_arg = (ground_truth && ground_truth != 'NO_FILE') ? "--ground-truth ${ground_truth}" : ''
    def ref_arg = (references_csv && references_csv != 'NO_FILE') ? "--references-csv ${references_csv}" : ''
    def run_id_arg = params.downstream_run_id ? "--run-id ${params.downstream_run_id}" : ''
    def coverage_label_arg = params.downstream_coverage_label ? "--coverage-label ${params.downstream_coverage_label}" : ''
    def quality_label_arg = params.downstream_quality_label ? "--quality-label ${params.downstream_quality_label}" : ''
    def differr_score_field_arg = params.downstream_differr_score_field ? "--differr-score-field ${params.downstream_differr_score_field}" : ''
    """
    exec > >(tee -a downstream_analysis.log) 2>&1

    echo "=== DOWNSTREAM_ANALYSIS started at \$(date) ==="
    echo "Modifications dir: ${modifications_dir}"
    echo "Ground truth:      ${ground_truth}"
    echo "References CSV:    ${references_csv}"
    echo "==========================================="

    python ${projectDir}/bin/downstream_analysis/run_analysis.py \\
        --input-dir ${modifications_dir} \\
        --output-dir downstream_analysis \\
        ${gt_arg} \\
        ${ref_arg} \\
        ${run_id_arg} \\
        ${coverage_label_arg} \\
        ${quality_label_arg} \\
        ${differr_score_field_arg} \\
        ${args}

    echo "=== DOWNSTREAM_ANALYSIS completed at \$(date) ==="
    """

    stub:
    """
    mkdir -p downstream_analysis/metrics
    mkdir -p downstream_analysis/tool_comparison
    mkdir -p downstream_analysis/replicate_analysis
    mkdir -p downstream_analysis/visualizations
    mkdir -p downstream_analysis/report

    touch downstream_analysis/metrics/all_metrics.csv
    touch downstream_analysis/tool_comparison/summary.csv
    touch downstream_analysis/report/report.html
    """
}
