#!/bin/bash -euo pipefail
samtools \
    fastq \
    -F 4 \
    --threads 1 \
     \
    test.sam \
    > test.fastq

cat <<-END_VERSIONS > versions.yml
"EXTRACT_MAPPED_READS":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
