process EPINANO_ERROR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container null

    input:
    tuple val(meta), path(native_bam), path(native_bai), path(ivt_bam), path(ivt_bai), path(reference)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "*.log"                      , emit: log, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_epinano_error"
    def zscore_threshold = task.ext.zscore ?: '3'
    def coverage_threshold = task.ext.coverage ?: '30'
    def error_threshold = task.ext.error_threshold ?: '0.1'
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}.log) 2>&1

    echo "=== EPINANO_ERROR started at \$(date) ==="
    echo "Sample: ${meta.id}"
    echo "Native BAM: $native_bam"
    echo "IVT BAM: $ivt_bam"
    echo "Reference: $reference"
    echo "Z-score threshold: $zscore_threshold"
    echo "Coverage threshold: $coverage_threshold"
    echo "Error threshold: $error_threshold"
    echo "CPUs: $task.cpus"
    echo "==========================================="

    mkdir -p ${prefix}

    # Clone EpiNano from GitHub
    if [ ! -d "EpiNano" ]; then
        git clone --depth 1 https://github.com/novoalab/EpiNano.git
    fi

    # Index reference if not already indexed
    if [ ! -f ${reference}.fai ]; then
        samtools faidx $reference
    fi

    # Extract variants for native (WT/modified) sample using EpiNano
    # Using -r (lowercase) for reference as per current EpiNano API
    python EpiNano/Epinano_Variants.py \\
        -r $reference \\
        -b $native_bam \\
        -c $task.cpus \\
        -o native_variants || true

    # Extract variants for IVT (KO/unmodified) sample
    python EpiNano/Epinano_Variants.py \\
        -r $reference \\
        -b $ivt_bam \\
        -c $task.cpus \\
        -o ivt_variants || true

    # Find the per-site CSV files - EpiNano outputs files like: sample.fwd.per.site.csv
    # Try multiple patterns to handle different EpiNano versions
    native_var=\$(find native_variants -name "*.per.site.csv" -o -name "*.per_site.*.csv" 2>/dev/null | head -1)
    ivt_var=\$(find ivt_variants -name "*.per.site.csv" -o -name "*.per_site.*.csv" 2>/dev/null | head -1)

    # If not found in output dirs, check current directory (older EpiNano versions)
    if [ -z "\$native_var" ]; then
        native_var=\$(find . -maxdepth 1 -name "*\$(basename $native_bam .bam)*.per.site.csv" 2>/dev/null | head -1)
    fi
    if [ -z "\$ivt_var" ]; then
        ivt_var=\$(find . -maxdepth 1 -name "*\$(basename $ivt_bam .bam)*.per.site.csv" 2>/dev/null | head -1)
    fi

    # Copy per-site CSV files to output - these contain RAW per-position data for ROC curves
    # Columns: #Ref, pos, base, strand, cov, q_mean, q_median, q_std, mis, ins, del
    if [ -f "\$native_var" ]; then
        cp "\$native_var" ${prefix}/native_per_site.csv
    fi
    if [ -f "\$ivt_var" ]; then
        cp "\$ivt_var" ${prefix}/ivt_per_site.csv
    fi

    # Fix column header: R treats '#' as comment, rename '#Ref' to 'X.Ref' (R's expected format)
    if [ -f "\$native_var" ]; then
        sed 's/^#Ref/X.Ref/' "\$native_var" > native_fixed.csv
    fi
    if [ -f "\$ivt_var" ]; then
        sed 's/^#Ref/X.Ref/' "\$ivt_var" > ivt_fixed.csv
    fi

    # Run differential error analysis if both per-site files exist
    # Epinano_DiffErr.R calculates sum_err internally from mis+ins+del columns
    # Output ALL positions with z-scores (set threshold to 0 to get all)
    if [ -f "native_fixed.csv" ] && [ -f "ivt_fixed.csv" ]; then
        Rscript EpiNano/Epinano_DiffErr.R \\
            -k "ivt_fixed.csv" \\
            -w "native_fixed.csv" \\
            -t $zscore_threshold \\
            -o ${prefix}/${meta.id} \\
            -c $coverage_threshold \\
            -f sum_err \\
            -d $error_threshold \\
            $args || true
    fi

    # Move any generated files to output directory
    mv *.csv ${prefix}/ 2>/dev/null || true
    mv *.pdf ${prefix}/ 2>/dev/null || true

    # Create header file if no diff_err results generated
    if [ ! -f "${prefix}/${meta.id}"*.csv ] 2>/dev/null; then
        echo -e "# EpiNano DiffErr analysis" > ${prefix}/${meta.id}_diff_err.csv
        echo -e "# No differential error results could be generated" >> ${prefix}/${meta.id}_diff_err.csv
        echo -e "# Check that per_site.csv files were created and have sufficient coverage" >> ${prefix}/${meta.id}_diff_err.csv
    fi

    echo "=== EPINANO_ERROR completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epinano: \$(cd EpiNano && git describe --tags 2>/dev/null || echo "1.2")
        python: \$(python --version 2>&1 | sed 's/Python //')
        R: \$(R --version 2>&1 | head -1 | sed 's/R version \\([^ ]*\\).*/\\1/')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_epinano_error"
    """
    mkdir -p ${prefix}
    touch ${prefix}/native_per_site.csv
    touch ${prefix}/ivt_per_site.csv
    touch ${prefix}/${meta.id}_diff_err.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epinano: 1.2
        python: 3.8.0
        R: 4.0.0
    END_VERSIONS
    """
}
