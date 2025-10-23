#!/bin/bash -euo pipefail
samtools view -S -b -h test.sam | samtools sort -o test.sorted_sorted.bam

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_SORT":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
