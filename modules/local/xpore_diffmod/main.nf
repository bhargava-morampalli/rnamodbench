process XPORE_DIFFMOD {
    tag "$key"
    label 'process_medium'  // rnamodbench nanoseq uses process_medium for xpore

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0"

    input:
    tuple val(key), val(types), path(dataprep_dirs)

    output:
    tuple val(key), path("${key}_diffmod"), emit: diffmod
    path "*.log"                          , emit: log, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def native_idx = types.findIndexOf { it == 'native' }
    def ivt_idx = types.findIndexOf { it == 'ivt' }
    def native_dir = dataprep_dirs[native_idx]
    def ivt_dir = dataprep_dirs[ivt_idx]
    def args = task.ext.args ?: ''
    def out_dir = "${key}_diffmod"
    """
    # Capture stdout and stderr to log file
    exec > >(tee -a ${key}_xpore.log) 2>&1

    echo "=== XPORE_DIFFMOD started at \$(date) ==="
    echo "Key: ${key}"
    echo "Native dir: ${native_dir}"
    echo "IVT dir: ${ivt_dir}"
    echo "Arguments: $args"
    echo "CPUs: $task.cpus"
    echo "==========================================="

    # Create xpore config YAML
    cat > config.yml << YAML
out: ${out_dir}
data:
  native:
    rep1: ${native_dir}
  ivt:
    rep1: ${ivt_dir}
YAML

    xpore diffmod --config config.yml --n_processes ${task.cpus} $args

    echo "=== XPORE_DIFFMOD completed at \$(date) ==="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed 's/^.*xpore //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def out_dir = "${key}_diffmod"
    """
    mkdir -p ${out_dir}
    touch ${out_dir}/diffmod.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: 2.1
    END_VERSIONS
    """
}
