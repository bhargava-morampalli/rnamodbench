#!/bin/bash -euo pipefail
samtools view \
    -b -F 4 \
    -@ 6 \
    -o test_meta.bam \
    test.sam

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_VIEW":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
