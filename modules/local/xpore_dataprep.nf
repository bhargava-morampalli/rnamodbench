process XPORE_DATAPREP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::xpore=2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xpore:2.1--pyh7cba7a3_0' :
        'biocontainers/xpore:2.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(eventalign)

    output:
    tuple val(meta), path("dataprep"), emit: dataprep
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    xpore dataprep \\
        --eventalign $eventalign \\
        --out_dir dataprep \\
        --n_processes ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed 's/^.*xpore //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dataprep
    touch dataprep/eventalign.index
    touch dataprep/data.index
    touch dataprep/data.json
    touch dataprep/data.readcount

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xpore: \$(xpore --version 2>&1 | sed 's/^.*xpore //; s/ .*\$//')
    END_VERSIONS
    """
}
