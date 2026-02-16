process DOWNSTREAM_ANALYSIS {
    tag "downstream_analysis"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container null
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path(modifications_dir)
    val(ground_truth)
    val(references_csv)

    output:
    path("downstream_analysis")              , emit: results
    path("downstream_analysis/metadata")     , emit: metadata, optional: true
    path("downstream_analysis/by_reference") , emit: by_reference, optional: true
    path("downstream_analysis/collation")    , emit: collation, optional: true
    path "*.log"                             , emit: log, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def threshold = task.ext.threshold ?: ''
    def threshold_arg = threshold ? "--threshold ${threshold}" : ''
    def min_reps = task.ext.min_replicates ?: '2'
    def expected_reps = task.ext.expected_replicates ?: '3'

    def gt_arg = (ground_truth && ground_truth != 'NO_FILE') ? "--ground-truth ${ground_truth}" : ''
    def ref_arg = (references_csv && references_csv != 'NO_FILE') ? "--references-csv ${references_csv}" : ''

    def run_id_arg = params.downstream_run_id ? "--run-id ${params.downstream_run_id}" : ''
    def coverage_label_arg = params.downstream_coverage_label ? "--coverage-label ${params.downstream_coverage_label}" : ''
    def quality_label_arg = params.downstream_quality_label ? "--quality-label ${params.downstream_quality_label}" : ''

    """
    exec > >(tee -a downstream_analysis.log) 2>&1

    echo "=== DOWNSTREAM_ANALYSIS started at \$(date) ==="
    echo "Modifications dir: $modifications_dir"
    echo "Ground truth: $ground_truth"
    echo "References CSV: $references_csv"
    echo "=============================================="

    python ${projectDir}/bin/downstream_analysis/run_analysis.py \\
        --input-dir "$modifications_dir" \\
        --output-dir downstream_analysis \\
        --min-replicates $min_reps \\
        --expected-replicates $expected_reps \\
        $threshold_arg \\
        $gt_arg \\
        $ref_arg \\
        $run_id_arg \\
        $coverage_label_arg \\
        $quality_label_arg \\
        $args

    echo "=== DOWNSTREAM_ANALYSIS completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        sklearn: \$(python -c "import sklearn; print(sklearn.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p downstream_analysis/metadata
    mkdir -p downstream_analysis/by_reference
    mkdir -p downstream_analysis/collation

    touch downstream_analysis/metadata/data_quality_report.csv
    touch downstream_analysis/collation/metrics_long.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10.0
        pandas: 2.0.0
        sklearn: 1.3.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}
