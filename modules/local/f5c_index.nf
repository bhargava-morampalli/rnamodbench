process F5C_INDEX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::f5c=1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/f5c:1.1--h0326b38_1' :
        'biocontainers/f5c:1.1--h0326b38_1' }"

    input:
    tuple val(meta), path(fast5_dir), path(fastq)

    output:
    tuple val(meta), path(fastq), path("*.index"), path("*.index.readdb"), path("*.index.gzi"), path("*.index.fai"), path(fast5_dir), emit: indexed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    f5c index -t ${task.cpus} -d $fast5_dir $fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        f5c: \$(f5c --version 2>&1 | sed 's/^.*f5c //; s/ .*\$//')
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
        f5c: \$(f5c --version 2>&1 | sed 's/^.*f5c //; s/ .*\$//')
    END_VERSIONS
    """
}
