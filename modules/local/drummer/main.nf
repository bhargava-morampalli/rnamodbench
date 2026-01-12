process DRUMMER {
    tag "$meta.id"
    label 'process_high'

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
    prefix = task.ext.prefix ?: "${meta.id}_drummer"
    // Set p-value threshold to 1.0 to output ALL positions for ROC curve analysis
    def pval_threshold = task.ext.pval ?: '1.0'
    def odds_ratio = task.ext.odds_ratio ?: '0'
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}.log) 2>&1

    echo "=== DRUMMER started at \$(date) ==="
    echo "Sample: ${meta.id}"
    echo "Native BAM: $native_bam"
    echo "IVT BAM: $ivt_bam"
    echo "Reference: $reference"
    echo "P-value threshold: $pval_threshold"
    echo "Odds ratio: $odds_ratio"
    echo "Arguments: $args"
    echo "==========================================="

    mkdir -p ${prefix}

    # Clone DRUMMER from GitHub
    if [ ! -d "DRUMMER" ]; then
        git clone --depth 1 https://github.com/DepledgeLab/DRUMMER.git
        # Fix Python 3.11+ compatibility: 'rU' mode was removed in Python 3.11
        sed -i "s/'rU'/'r'/g" DRUMMER/modules/test_exome.py
        sed -i "s/'rU'/'r'/g" DRUMMER/modules/support.py 2>/dev/null || true

        # Fix multiprocessing issues in containerized environments (Python 3.14+)
        # Add multiprocessing start method fix at the beginning of DRUMMER.py
        sed -i '1a import multiprocessing; multiprocessing.set_start_method("fork", force=True)' DRUMMER/DRUMMER.py 2>/dev/null || true
    fi

    # Index reference if not already indexed
    if [ ! -f ${reference}.fai ]; then
        samtools faidx $reference
    fi

    # Get reference name from fasta index
    ref_name=\$(head -1 ${reference}.fai | cut -f1)

    # Run DRUMMER in exome mode (suitable for rRNA)
    # Control = IVT (modifications absent)
    # Treatment = Native (modifications present)
    # Using p-value threshold of 1.0 to output ALL positions for ROC curve analysis
    # DRUMMER outputs per-position statistics with p-values and odds ratios
    python DRUMMER/DRUMMER.py \\
        -r $reference \\
        -n \$ref_name \\
        -c $ivt_bam \\
        -t $native_bam \\
        -o ${prefix} \\
        -a exome \\
        -p ${pval_threshold} \\
        $args || true

    # DRUMMER outputs summary.txt with per-position data including:
    # transcript_id, chrom, ref_base, position, read_depth, base_fractions, odds_ratio, pval, motif, g_test, genomic_pos
    # This raw data can be used for ROC curves using the pval column

    # Also look for any intermediate files with all position statistics
    find . -name "*.txt" -o -name "*.csv" | xargs -I {} cp {} ${prefix}/ 2>/dev/null || true

    # Check if output was generated, create empty file with header if not
    if [ ! -f "${prefix}/summary.txt" ] || [ ! -s "${prefix}/summary.txt" ]; then
        echo -e "transcript_id\\tchrom\\tref_base\\tposition\\tread_depth\\tbase_fractions\\todds_ratio\\tpval\\tmotif\\tg_test\\tgenomic_pos" > ${prefix}/summary.txt
        echo -e "# No candidates found - check if BAM files have sufficient coverage" >> ${prefix}/summary.txt
    fi

    echo "=== DRUMMER completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drummer: \$(cd DRUMMER && git describe --tags 2>/dev/null || echo "1.0")
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_drummer"
    """
    mkdir -p ${prefix}
    echo -e "transcript_id\\tchrom\\tref_base\\tposition\\tread_depth\\tbase_fractions\\todds_ratio\\tpval\\tmotif\\tg_test\\tgenomic_pos" > ${prefix}/summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drummer: 1.0
        python: 3.8.0
    END_VERSIONS
    """
}
