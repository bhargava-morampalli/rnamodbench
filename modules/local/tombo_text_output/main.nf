process TOMBO_TEXT_OUTPUT {
    tag "$key"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // No pre-built container available for ont-tombo - Wave will build from conda environment
    container null

    input:
    tuple val(key), path(statistic), path(reference)

    output:
    tuple val(key), path("*.csv"), emit: csv
    path "*.log"                 , emit: log, optional: true
    tuple val("${task.process}"), val('tombo'), eval('tombo --version 2>&1 | grep -oP \'[0-9]+\\.[0-9]+[0-9.]*\' | head -1 || echo unknown'), topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c \'import pandas; print(pandas.__version__)\' 2>/dev/null || echo unknown'), topic: versions

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

    """

    stub:
    def prefix = task.ext.prefix ?: "${key}"
    """
    touch ${prefix}.csv

    """
}
