process TOMBO_DETECT_MODIFICATIONS {
    tag "$key"
    label 'process_high'

    conda "bioconda::ont-tombo=1.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-tombo:1.5.1--py_0' :
        'biocontainers/ont-tombo:1.5.1--py_0' }"

    input:
    tuple val(key), val(types), path(fast5_dirs)

    output:
    tuple val(key), path("*.tombo.stats"), emit: stats
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
    tombo detect_modifications level_sample_compare \\
        --fast5-basedirs ${native_dir} \\
        --alternate-fast5-basedirs ${ivt_dir} \\
        --statistics-file-basename ${prefix} \\
        --store-p-value \\
        --minimum-test-reads 1 \\
        --statistic-type ks \\
        --processes ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tombo: \$(tombo --version 2>&1 | sed 's/^.*tombo version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${key}"
    """
    touch ${prefix}.tombo.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tombo: \$(tombo --version 2>&1 | sed 's/^.*tombo version //; s/ .*\$//')
    END_VERSIONS
    """
}
