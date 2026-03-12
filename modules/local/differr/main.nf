process DIFFERR {
    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container null

    input:
    tuple val(meta), path(native_bam), path(native_bai), path(ivt_bam), path(ivt_bai), path(reference)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val(meta), path("${prefix}.hdf5"), emit: hdf5, optional: true
    path "*.log"                           , emit: log, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_differr"
    // Set FDR threshold to 1.0 to output ALL positions for ROC curve analysis
    def fdr_threshold = task.ext.fdr ?: '1.0'
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}.log) 2>&1

    echo "=== DIFFERR started at \$(date) ==="
    echo "Sample: ${meta.id}"
    echo "Native BAM: $native_bam"
    echo "IVT BAM: $ivt_bam"
    echo "Reference: $reference"
    echo "FDR threshold: $fdr_threshold"
    echo "Arguments: $args"
    echo "==========================================="

    if ! command -v differr >/dev/null 2>&1; then
        echo "ERROR: 'differr' command not found in environment. Install DiffErr in the module environment/container."
        exit 1
    fi

    # Index reference if not already indexed
    if [ ! -f ${reference}.fai ]; then
        samtools faidx $reference
    fi

    # Run differr comparing native (WT/modified) vs IVT (control/unmodified)
    # -b: native/WT samples (modifications present)
    # -a: IVT samples (modifications absent)
    # Using FDR threshold of 1.0 to output ALL positions for ROC curve analysis
    differr \\
        -b $native_bam \\
        -a $ivt_bam \\
        -r $reference \\
        -o ${prefix}.bed \\
        -f ${fdr_threshold} \\
        $args

    # Also try to get raw output files if differr creates intermediate files
    find . -name "*.bed" -o -name "*.csv" -o -name "*.tsv" | xargs -I {} cp {} . 2>/dev/null || true

    # Create header if no output
    if [ ! -f "${prefix}.bed" ] || [ ! -s "${prefix}.bed" ]; then
        echo -e "chr\\tstart\\tend\\todds_ratio\\tg_stat\\tpval\\tfdr" > ${prefix}.bed
        echo -e "# No positions found - check if BAM files have sufficient coverage" >> ${prefix}.bed
    fi

    echo "=== DIFFERR completed at \$(date) ==="

    # differr does not support --version in some releases; keep probe non-fatal
    differr_version=\$((differr --version 2>/dev/null || true) | head -1 | sed -E 's/^[^0-9]*([0-9]+(\\.[0-9]+)*)?.*\$/\\1/')
    if [ -z "\$differr_version" ]; then
        differr_version="unknown"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        differr: \$differr_version
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_differr"
    """
    # Create stub BED file with expected columns
    echo -e "chr\\tstart\\tend\\todds_ratio\\tg_stat\\tpval\\tfdr" > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        differr: 0.2
        python: 3.8.0
    END_VERSIONS
    """
}
