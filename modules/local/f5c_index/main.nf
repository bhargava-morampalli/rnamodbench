process F5C_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::f5c=1.5"
    container "quay.io/biocontainers/f5c:1.5--hee927d3_2"

    input:
    tuple val(meta), path(fast5_dir), path(fastq)
    path sequencing_summary  // Optional sequencing summary for faster indexing per benchmarking papers

    output:
    tuple val(meta), path(fastq), path("*.index"), path("*.index.readdb"), path("*.index.gzi"), path("*.index.fai"), path(fast5_dir), emit: indexed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Add -s for sequencing summary if provided, and --iop for I/O threads per benchmarking papers
    def summary_arg = sequencing_summary.name != 'NO_FILE' ? "-s ${sequencing_summary}" : ""
    """
    f5c index -t ${task.cpus} --iop ${task.cpus} ${summary_arg} -d $fast5_dir $fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        f5c: \$(f5c --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+[0-9.]*' | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${fastq}.index
    touch ${fastq}.index.readdb
    touch ${fastq}.index.gzi
    touch ${fastq}.index.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        f5c: 1.1
    END_VERSIONS
    """
}
