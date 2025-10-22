process TOMBO_TEXT_OUTPUT {
    tag "$key"
    label 'process_low'

    conda "bioconda::ont-tombo=1.5.1 conda-forge::pandas=1.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ont-tombo:1.5.1--py_0' :
        'biocontainers/ont-tombo:1.5.1--py_0' }"

    input:
    tuple val(key), path(statistic)

    output:
    tuple val(key), path("*.csv"), emit: bed
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${key}"
    def rrna_type = key.split('_')[0]
    // Define chromosome names and coordinates based on rRNA type
    def chrom = rrna_type == '16s' ? '16s_88_rrsE' : '23s_78_rrlB'
    def start = 1
    def end = rrna_type == '16s' ? 1813 : 3163
    """
    #!/usr/bin/env python
    from tombo import tombo_stats
    import pandas as pd

    # Load Tombo statistics
    sample_level_stats = tombo_stats.LevelStats("$statistic")
    
    # Get regional statistics
    reg_level_stats = sample_level_stats.get_reg_stats('${chrom}', '+', ${start}, ${end})
    
    # Convert to DataFrame and save as CSV
    pd.DataFrame(reg_level_stats).to_csv("${prefix}.csv", index=False)

    # Create versions file
    import tombo
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    tombo: {tombo.__version__}\\n')
        f.write(f'    pandas: {pd.__version__}\\n')
    """

    stub:
    def prefix = task.ext.prefix ?: "${key}"
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tombo: 1.5.1
        pandas: 1.3.4
    END_VERSIONS
    """
}
