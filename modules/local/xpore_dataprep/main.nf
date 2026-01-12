process XPORE_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0"

    input:
    tuple val(meta), path(eventalign)

    output:
    tuple val(meta), path("dataprep_*"), emit: dataprep
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def out_dir = "dataprep_${meta.id}"
    """
    xpore dataprep \\
        --eventalign $eventalign \\
        --out_dir $out_dir \\
        --n_processes ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed 's/^.*xpore //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def out_dir = "dataprep_${meta.id}"
    """
    mkdir -p $out_dir
    touch $out_dir/eventalign.index
    touch $out_dir/data.index
    touch $out_dir/data.json
    touch $out_dir/data.readcount

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: 2.1
    END_VERSIONS
    """
}
