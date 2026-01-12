process FAST5_SUBSET {
    tag "$meta.id"
    label 'process_low'  // Single-threaded I/O-bound tool; extra CPUs won't help

    conda "bioconda::ont-fast5-api=4.1.0"
    container "quay.io/biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(read_ids), path(fast5_dir)

    output:
    tuple val(meta), path("fast5_subset"), emit: fast5
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fast5_subset \\
        $args \\
        --input ${fast5_dir} \\
        --read_id_list ${read_ids} \\
        --save_path fast5_subset \\
        --recursive

    echo -e '"${task.process}":\n  ont-fast5-api: \$(fast5_subset --version | sed -nE "s/fast5_subset, version (.*)/\\1/p")' > versions.yml
    """

    stub:
    """
    mkdir -p fast5_subset

    echo -e '"${task.process}":\n  ont-fast5-api: \$(fast5_subset --version | sed -nE "s/fast5_subset, version (.*)/\\1/p")' > versions.yml
    """
}