process XPORE_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/xpore:2.1--pyh5e36f6f_0"

    input:
    tuple val(meta), path(eventalign)

    output:
    tuple val(meta), path("dataprep_*"), emit: dataprep
    tuple val("${task.process}"), val('xpore'), eval('xpore --version 2>&1 | sed \'s/^.*xpore //; s/ .*$//\' || echo unknown'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def out_dir = "dataprep_${meta.id}"
    """
    xpore dataprep \\
        --eventalign $eventalign \\
        --out_dir $out_dir \\
        --n_processes ${task.cpus}
    """

    stub:
    def out_dir = "dataprep_${meta.id}"
    """
    mkdir -p $out_dir
    touch $out_dir/eventalign.index
    touch $out_dir/data.index
    touch $out_dir/data.json
    touch $out_dir/data.readcount
    """
}
