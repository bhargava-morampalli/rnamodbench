process XPORE_DIFFMOD {
    tag "$key"
    label 'process_high'

    conda "bioconda::xpore=2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh7cba7a3_0' :
        'biocontainers/xpore:2.1--pyh7cba7a3_0' }"

    input:
    tuple val(key), val(types), path(dataprep_dirs)

    output:
    tuple val(key), path("diffmod"), emit: diffmod
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def native_idx = types.findIndexOf { it == 'native' }
    def ivt_idx = types.findIndexOf { it == 'ivt' }
    def native_dir = dataprep_dirs[native_idx]
    def ivt_dir = dataprep_dirs[ivt_idx]
    """
    # Create xpore config YAML
    cat > config.yml << YAML
out_dir: diffmod
data:
  native:
    rep1: ${native_dir}
  ivt:
    rep1: ${ivt_dir}
YAML

    xpore diffmod --config config.yml --n_processes ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed 's/^.*xpore //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p diffmod
    touch diffmod/diffmod.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed 's/^.*xpore //; s/ .*\$//')
    END_VERSIONS
    """
}
