#!/bin/bash -euo pipefail
touch test_stub.fastq

cat <<-END_VERSIONS > versions.yml
"EXTRACT_MAPPED_READS":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
