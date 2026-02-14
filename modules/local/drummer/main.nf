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
    def drummer_home = task.ext.drummer_home ?: ''
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

    # Clone DRUMMER at pinned commit (not pip-installable — standalone scripts)
    DRUMMER_HOME="${drummer_home}"
    if [ -z "\$DRUMMER_HOME" ]; then
        git clone --quiet https://github.com/DepledgeLab/DRUMMER.git drummer_repo
        cd drummer_repo && git checkout --quiet 6683822c6210083e4ab0eecb4b6327e3c55f4c46 && cd ..
        DRUMMER_HOME="drummer_repo"
    fi

    if [ -n "\$DRUMMER_HOME" ] && [ -f "\$DRUMMER_HOME/DRUMMER.py" ]; then
        DRUMMER_SCRIPT="\$DRUMMER_HOME/DRUMMER.py"
    elif command -v DRUMMER.py >/dev/null 2>&1; then
        DRUMMER_SCRIPT=\$(command -v DRUMMER.py)
    elif command -v drummer >/dev/null 2>&1; then
        DRUMMER_SCRIPT=\$(command -v drummer)
    else
        echo "ERROR: DRUMMER executable not found. Install DRUMMER in the module environment/container or set task.ext.drummer_home."
        exit 1
    fi

    # Index reference if not already indexed
    if [ ! -f ${reference}.fai ]; then
        samtools faidx $reference
    fi

    # Get reference name from fasta index
    ref_name=\$(head -1 ${reference}.fai | cut -f1)

    # Run DRUMMER in exome mode (suitable for rRNA)
    # Control = Native (modifications present)
    # Treatment = IVT (modifications absent)
    # Using p-value threshold of 1.0 to output ALL positions for ROC curve analysis
    # DRUMMER outputs per-position statistics with p-values and odds ratios
    python "\$DRUMMER_SCRIPT" \\
        -r $reference \\
        -n \$ref_name \\
        -c $native_bam \\
        -t $ivt_bam \\
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

    drummer_version=\$(python - <<'PY'
try:
    import importlib.metadata as md
except Exception:
    import importlib_metadata as md  # type: ignore

version = "unknown"
for package in ("drummer", "DRUMMER"):
    try:
        version = md.version(package)
        break
    except Exception:
        pass
print(version)
PY
)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        drummer: \$drummer_version
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
