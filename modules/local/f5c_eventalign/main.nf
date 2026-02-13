process F5C_EVENTALIGN {
    tag "$meta.id"
    label 'process_medium'  // rnamodbench nanopolish equivalent uses process_medium

    conda "bioconda::f5c=1.5"
    container "quay.io/biocontainers/f5c:1.5--hee927d3_2"

    input:
    tuple val(meta), path(fastq), path(fastq_idx), path(readdb), path(gzi), path(fai), path(bam), path(bai), path(fast5), path(reference)

    output:
    tuple val(meta), path("*.txt"), emit: eventalign
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Added --secondary=yes --min-mapq 0 to match Nanopolish defaults per benchmarking papers
    def args = task.ext.args ?: '--rna --scale-events --signal-index --print-read-names --secondary=yes --min-mapq 0'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    f5c eventalign \\
        $args \\
        -t $task.cpus \\
        -b $bam \\
        -g $reference \\
        -r $fastq \\
        > ${prefix}_eventalign.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        f5c: \$(f5c --version 2>&1 | sed 's/^.*f5c //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_eventalign.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        f5c: 1.1
    END_VERSIONS
    """
}
