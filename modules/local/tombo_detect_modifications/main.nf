process TOMBO_DETECT_MODIFICATIONS {
    tag "$key"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // No pre-built container available for ont-tombo - Wave will build from conda environment
    container null

    input:
    tuple val(key), val(types), path(fast5_dirs)

    output:
    tuple val(key), path("*.tombo.stats"), emit: stats
    path "*.log"                         , emit: log, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${key}"
    def native_idx = types.findIndexOf { it == 'native' }
    def ivt_idx = types.findIndexOf { it == 'ivt' }
    def native_dir = fast5_dirs[native_idx]
    def ivt_dir = fast5_dirs[ivt_idx]
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${prefix}_tombo_detect.log) 2>&1

    echo "=== TOMBO_DETECT_MODIFICATIONS started at \$(date) ==="
    echo "Key: ${key}"
    echo "Native dir: ${native_dir}"
    echo "IVT dir: ${ivt_dir}"
    echo "CPUs: $task.cpus"
    echo "==========================================="

    tombo detect_modifications level_sample_compare \\
        --fast5-basedirs ${native_dir} \\
        --alternate-fast5-basedirs ${ivt_dir} \\
        --statistics-file-basename ${prefix} \\
        --store-p-value \\
        --minimum-test-reads 1 \\
        --statistic-type ks \\
        --processes ${task.cpus}

    echo "=== TOMBO_DETECT_MODIFICATIONS completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tombo: \$(tombo --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+[0-9.]*' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${key}"
    """
    touch ${prefix}.tombo.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tombo: 1.5.1
    END_VERSIONS
    """
}
