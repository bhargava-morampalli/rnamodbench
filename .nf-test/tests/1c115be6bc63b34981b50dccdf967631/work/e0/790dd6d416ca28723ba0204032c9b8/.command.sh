#!/bin/bash -euo pipefail
samtools depth -a -m 0 test.bam > test.txt

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_DEPTH":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
