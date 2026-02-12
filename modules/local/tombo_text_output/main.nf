process TOMBO_TEXT_OUTPUT {
    tag "$key"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // No pre-built container available for ont-tombo - Wave will build from conda environment
    container null

    input:
    tuple val(key), path(statistic), path(reference)

    output:
    tuple val(key), path("*.csv"), emit: bed
    path "*.log"                 , emit: log, optional: true
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${key}"
    """
    #!/usr/bin/env python
    import sys
    import logging
    from datetime import datetime

    # Set up logging to file and stdout
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('${prefix}_tombo_text.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    logger.info("=== TOMBO_TEXT_OUTPUT started ===")
    logger.info(f"Key: ${key}")
    logger.info(f"Statistic file: $statistic")
    logger.info(f"Reference file: ${reference}")

    # Parse reference FASTA to get chromosome name and sequence length
    chrom = None
    seq_length = 0
    with open("${reference}") as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                chrom = line[1:].split()[0]
            else:
                seq_length += len(line)

    if not chrom:
        raise ValueError(f"No sequence found in reference FASTA: ${reference}")

    logger.info(f"Chromosome: {chrom}")
    logger.info(f"Region: 1-{seq_length}")

    from tombo import tombo_stats
    import pandas as pd

    try:
        # Load Tombo statistics
        logger.info("Loading Tombo statistics...")
        sample_level_stats = tombo_stats.LevelStats("$statistic")

        # Get regional statistics for the full reference
        logger.info("Extracting regional statistics...")
        reg_level_stats = sample_level_stats.get_reg_stats(chrom, '+', 1, seq_length)

        # Convert to DataFrame and save as CSV
        logger.info(f"Saving results to ${prefix}.csv")
        pd.DataFrame(reg_level_stats).to_csv("${prefix}.csv", index=False)
        logger.info("=== TOMBO_TEXT_OUTPUT completed successfully ===")
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise

    # Create versions file
    import subprocess
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        tombo_version = subprocess.check_output(['tombo', '--version'], stderr=subprocess.STDOUT).decode().strip().split()[-1]
        f.write('    tombo: ' + tombo_version + '\\n')
        f.write('    pandas: ' + pd.__version__ + '\\n')
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
