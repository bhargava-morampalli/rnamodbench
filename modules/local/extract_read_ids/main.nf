process EXTRACT_READ_IDS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.read_ids.txt"), emit: read_ids
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit seq -n -i $fastq > ${prefix}.read_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.read_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit version //')
    END_VERSIONS
    """
}