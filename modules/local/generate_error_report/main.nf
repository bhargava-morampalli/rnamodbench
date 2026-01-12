process GENERATE_ERROR_REPORT {
    tag "error_report"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container null

    input:
    path outdir
    path logs

    output:
    path "error_summary.html", emit: html
    path "error_summary.csv" , emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Copy the error report generator script
    cp ${projectDir}/bin/generate_error_report.py .

    # Generate the error summary report
    python generate_error_report.py \\
        --outdir $outdir \\
        --output error_summary

    # Move reports to current directory
    mv ${outdir}/pipeline_info/error_summary.html . 2>/dev/null || touch error_summary.html
    mv ${outdir}/pipeline_info/error_summary.csv . 2>/dev/null || touch error_summary.csv
    """

    stub:
    """
    touch error_summary.html
    touch error_summary.csv
    """
}
