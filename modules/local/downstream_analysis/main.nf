process DOWNSTREAM_ANALYSIS {
    tag "downstream_analysis"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container null

    input:
    path(tombo_outputs)
    path(yanocomp_outputs)
    path(nanocompore_outputs)
    path(xpore_outputs)
    path(eligos_outputs)
    path(epinano_outputs)
    path(differr_outputs)
    path(drummer_outputs)
    path(jacusa2_outputs)
    path(ground_truth)

    output:
    path("downstream_analysis")         , emit: results
    path("downstream_analysis/report")  , emit: report, optional: true
    path("downstream_analysis/metrics") , emit: metrics, optional: true
    path "*.log"                        , emit: log, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gt_arg = ground_truth.name != 'NO_FILE' ? "--ground-truth ${ground_truth}" : ''
    def threshold = task.ext.threshold ?: ''
    def threshold_arg = threshold ? "--threshold ${threshold}" : ''
    def min_reps = task.ext.min_replicates ?: '2'
    def expected_reps = task.ext.expected_replicates ?: '3'

    def gt_arg = (ground_truth && ground_truth != 'NO_FILE') ? "--ground-truth ${ground_truth}" : ''
    def ref_arg = (references_csv && references_csv != 'NO_FILE') ? "--references-csv ${references_csv}" : ''

    def run_id_arg = params.downstream_run_id ? "--run-id ${params.downstream_run_id}" : ''
    def coverage_label_arg = params.downstream_coverage_label ? "--coverage-label ${params.downstream_coverage_label}" : ''
    def quality_label_arg = params.downstream_quality_label ? "--quality-label ${params.downstream_quality_label}" : ''
    def differr_score_field_arg = params.downstream_differr_score_field ? "--differr-score-field ${params.downstream_differr_score_field}" : ''

    """
    # Capture stdout and stderr to log file
    exec > >(tee -a downstream_analysis.log) 2>&1

    echo "=== DOWNSTREAM_ANALYSIS started at \$(date) ==="
    echo "Tombo outputs: $tombo_outputs"
    echo "Yanocomp outputs: $yanocomp_outputs"
    echo "Nanocompore outputs: $nanocompore_outputs"
    echo "xPore outputs: $xpore_outputs"
    echo "ELIGOS outputs: $eligos_outputs"
    echo "EpiNano outputs: $epinano_outputs"
    echo "DiffErr outputs: $differr_outputs"
    echo "DRUMMER outputs: $drummer_outputs"
    echo "JACUSA2 outputs: $jacusa2_outputs"
    echo "Ground truth: $ground_truth"
    echo "==========================================="

    # Build tool arguments
    TOOL_ARGS=""

    # Process each tool's outputs
    if [ -e "$tombo_outputs" ] && [ "$tombo_outputs" != "NO_FILE" ]; then
        for f in $tombo_outputs/*.csv; do
            if [ -f "\$f" ]; then
                TOOL_ARGS="\$TOOL_ARGS --tombo \$f"
                break
            fi
        done
    fi

    if [ -e "$yanocomp_outputs" ] && [ "$yanocomp_outputs" != "NO_FILE" ]; then
        for f in $yanocomp_outputs/*.bed; do
            if [ -f "\$f" ]; then
                TOOL_ARGS="\$TOOL_ARGS --yanocomp \$f"
                break
            fi
        done
    fi

    if [ -e "$nanocompore_outputs" ] && [ "$nanocompore_outputs" != "NO_FILE" ]; then
        for d in $nanocompore_outputs/*_sampcomp; do
            if [ -d "\$d" ]; then
                TOOL_ARGS="\$TOOL_ARGS --nanocompore \$d"
                break
            fi
        done
    fi

    if [ -e "$xpore_outputs" ] && [ "$xpore_outputs" != "NO_FILE" ]; then
        for d in $xpore_outputs/*_diffmod; do
            if [ -d "\$d" ]; then
                TOOL_ARGS="\$TOOL_ARGS --xpore \$d"
                break
            fi
        done
    fi

    if [ -e "$eligos_outputs" ] && [ "$eligos_outputs" != "NO_FILE" ]; then
        for d in $eligos_outputs/*_eligos; do
            if [ -d "\$d" ]; then
                TOOL_ARGS="\$TOOL_ARGS --eligos \$d"
                break
            fi
        done
    fi

    if [ -e "$epinano_outputs" ] && [ "$epinano_outputs" != "NO_FILE" ]; then
        for d in $epinano_outputs/*_epinano*; do
            if [ -d "\$d" ]; then
                TOOL_ARGS="\$TOOL_ARGS --epinano \$d"
                break
            fi
        done
    fi

    if [ -e "$differr_outputs" ] && [ "$differr_outputs" != "NO_FILE" ]; then
        for f in $differr_outputs/*.bed; do
            if [ -f "\$f" ]; then
                TOOL_ARGS="\$TOOL_ARGS --differr \$f"
                break
            fi
        done
    fi

    if [ -e "$drummer_outputs" ] && [ "$drummer_outputs" != "NO_FILE" ]; then
        for d in $drummer_outputs/*_drummer; do
            if [ -d "\$d" ]; then
                TOOL_ARGS="\$TOOL_ARGS --drummer \$d"
                break
            fi
        done
    fi

    if [ -e "$jacusa2_outputs" ] && [ "$jacusa2_outputs" != "NO_FILE" ]; then
        for f in $jacusa2_outputs/*.bed; do
            if [ -f "\$f" ]; then
                TOOL_ARGS="\$TOOL_ARGS --jacusa2 \$f"
                break
            fi
        done
    fi

    echo "Tool arguments: \$TOOL_ARGS"

    # Run downstream analysis
    python ${projectDir}/bin/downstream_analysis/run_analysis.py \\
        \$TOOL_ARGS \\
        $gt_arg \\
        $ref_arg \\
        $run_id_arg \\
        $coverage_label_arg \\
        $quality_label_arg \\
        $differr_score_field_arg \\
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
    mkdir -p downstream_analysis/metrics
    mkdir -p downstream_analysis/tool_comparison
    mkdir -p downstream_analysis/replicate_analysis
    mkdir -p downstream_analysis/visualizations
    mkdir -p downstream_analysis/report

    touch downstream_analysis/metrics/all_metrics.csv
    touch downstream_analysis/tool_comparison/summary.csv
    touch downstream_analysis/report/report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10.0
        pandas: 2.0.0
        sklearn: 1.3.0
        matplotlib: 3.7.0
    END_VERSIONS
    """
}
