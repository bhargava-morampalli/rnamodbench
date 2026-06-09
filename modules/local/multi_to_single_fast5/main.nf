process MULTI_TO_SINGLE_FAST5 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ont-fast5-api=4.1.0"
    container "quay.io/biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(fast5_dir)

    output:
    tuple val(meta), path("single_fast5_*"), emit: fast5
    tuple val("${task.process}"), val('ont-fast5-api'), eval('multi_to_single_fast5 --version 2>&1 | grep -oP \'[0-9]+\\.[0-9]+[0-9.]*\' | head -1 || echo unknown'), emit: versions_ont_fast5_api, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def out_dir = "single_fast5_${meta.id}"
    """
    multi_to_single_fast5 \\
        $args \\
        --input_path ${fast5_dir} \\
        --save_path $out_dir \\
        --recursive
    """

    stub:
    def out_dir = "single_fast5_${meta.id}"
    """
    mkdir -p $out_dir
    """
}