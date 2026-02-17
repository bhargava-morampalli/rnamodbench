process NANOCOMPORE_SAMPCOMP {
    tag "$key"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/nanocompore:1.0.4--pyhdfd78af_0"

    input:
    tuple val(key), val(condition1_label), val(condition2_label), path(condition1_files), path(condition2_files), path(fasta)

    output:
    tuple val(key), path("${prefix}"), emit: results
    path "*.log"                      , emit: log, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Output ALL positions for ROC curve analysis (pvalue_thr=1), with 5-mer context and logistic regression
    def args = task.ext.args ?: '--sequence_context 2 --pvalue_thr 1 --logit'
    prefix = task.ext.prefix ?: "${key}_sampcomp"
    def cond1_files = condition1_files.collect{ "${it}/out_eventalign_collapse.tsv" }.join(',')
    def cond2_files = condition2_files.collect{ "${it}/out_eventalign_collapse.tsv" }.join(',')
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}.log) 2>&1

    echo "=== NANOCOMPORE_SAMPCOMP started at \$(date) ==="
    echo "Key: ${key}"
    echo "Condition1 (${condition1_label}): ${cond1_files}"
    echo "Condition2 (${condition2_label}): ${cond2_files}"
    echo "Reference: $fasta"
    echo "Arguments: $args"
    echo "CPUs: $task.cpus"
    echo "==========================================="

    nanocompore sampcomp \\
        --file_list1 ${cond1_files} \\
        --file_list2 ${cond2_files} \\
        --label1 ${condition1_label} \\
        --label2 ${condition2_label} \\
        --fasta $fasta \\
        --outpath ${prefix} \\
        --nthreads ${task.cpus} \\
        $args

    echo "=== NANOCOMPORE_SAMPCOMP completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocompore: \$(nanocompore --version 2>/dev/null | grep -oP '[0-9]+\\.[0-9]+[0-9.]*' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${key}_sampcomp"
    """
    mkdir -p ${prefix}
    touch ${prefix}/outSampComp.db
    touch ${prefix}/outSampComp_results.tsv
    touch ${prefix}/outSampComp_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocompore: \$(nanocompore --version 2>/dev/null | grep -oP '[0-9]+\\.[0-9]+[0-9.]*' | head -1 || echo "unknown")
    END_VERSIONS
    """
}
