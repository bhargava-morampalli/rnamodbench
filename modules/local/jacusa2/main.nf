process JACUSA2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container null

    input:
    tuple val(meta), path(native_bam), path(native_bai), path(ivt_bam), path(ivt_bai), path(reference)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path "*.log"                          , emit: log, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_jacusa2"
    // Set minimum score to 0 to output ALL positions for ROC curve analysis
    def min_score = task.ext.min_score ?: '0'
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}.log) 2>&1

    echo "=== JACUSA2 started at \$(date) ==="
    echo "Sample: ${meta.id}"
    echo "Native BAM: $native_bam"
    echo "IVT BAM: $ivt_bam"
    echo "Reference: $reference"
    echo "Min score: $min_score"
    echo "Arguments: $args"
    echo "CPUs: $task.cpus"
    echo "==========================================="

    # Index reference if not already indexed
    if [ ! -f ${reference}.fai ]; then
        samtools faidx $reference
    fi

    # Ensure MD field is populated in BAM files
    samtools calmd -b $native_bam $reference > native_md.bam 2>/dev/null
    samtools index native_md.bam

    samtools calmd -b $ivt_bam $reference > ivt_md.bam 2>/dev/null
    samtools index ivt_md.bam

    # Find JACUSA2 jar location (varies by installation)
    # Search in common locations: CONDA_PREFIX, /opt/conda, /usr
    JACUSA_JAR=""
    for search_dir in "\${CONDA_PREFIX:-/opt/conda}" /opt/conda /usr; do
        if [ -d "\$search_dir" ]; then
            JACUSA_JAR=\$(find "\$search_dir" -name "JACUSA*.jar" 2>/dev/null | head -1)
            if [ -n "\$JACUSA_JAR" ]; then
                break
            fi
        fi
    done

    if [ -z "\$JACUSA_JAR" ]; then
        echo "ERROR: Could not find JACUSA2 jar file" >&2
        exit 1
    fi

    # Run JACUSA2 call-2 for pairwise comparison
    # Condition 1: Native (modified)
    # Condition 2: IVT (unmodified/control)
    # Using -T 0 (min score threshold) and -A to output ALL positions for ROC curve analysis
    # -A: Show all sites including sites without variants
    # Output includes per-position scores useful for ROC curves
    # Note: JACUSA2 uses -R for reference FASTA (not -r which is for result file!)
    java -jar \$JACUSA_JAR call-2 \\
        -R $reference \\
        -p $task.cpus \\
        -T ${min_score} \\
        -A \\
        -r ${prefix}.bed \\
        $args \\
        native_md.bam \\
        ivt_md.bam

    # Clean up temporary files
    rm -f native_md.bam native_md.bam.bai ivt_md.bam ivt_md.bam.bai

    echo "=== JACUSA2 completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jacusa2: \$(java -jar \$JACUSA_JAR 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || echo "2.0.4")
        java: \$(java -version 2>&1 | head -1 | sed 's/.*version "\\([^"]*\\)".*/\\1/')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_jacusa2"
    """
    # Create stub BED file with JACUSA2 output format
    echo -e "##JACUSA2 v2.0.4" > ${prefix}.bed
    echo -e "#contig\\tstart\\tend\\tname\\tscore\\tstrand\\tinfo\\tfilter\\tref" >> ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jacusa2: 2.0.4
        java: 1.8.0
    END_VERSIONS
    """
}
