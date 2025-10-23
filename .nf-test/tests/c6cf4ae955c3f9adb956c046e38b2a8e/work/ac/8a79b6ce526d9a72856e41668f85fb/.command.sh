#!/bin/bash -euo pipefail
samtools \
    index \
    -@ 2 \
     \
    test.bam

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_INDEX":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
