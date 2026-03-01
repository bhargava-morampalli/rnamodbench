process GENERATE_ERROR_REPORT {
    tag "${run_dir}"
    label 'process_single'

    conda "${projectDir}/modules/local/downstream_analysis/environment.yml"
    container null

    input:
    val run_dir

    output:
    path "process_status.tsv", emit: process_status
    path "log_events.tsv", emit: log_events
    path "tool_availability_per_run.tsv", emit: availability
    path "error_summary.html", emit: html
    path "error_summary.csv" , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "=== GENERATE_ERROR_REPORT started at \$(date) ==="
    echo "Run directory: $run_dir"

    python ${projectDir}/bin/generate_error_report.py \\
        --run-dir "$run_dir" \\
        --pipeline-info-dir . \\
        --output error_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch process_status.tsv
    touch log_events.tsv
    touch tool_availability_per_run.tsv
    touch error_summary.html
    touch error_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.10.0
        pandas: 2.0.0
    END_VERSIONS
    """
}
