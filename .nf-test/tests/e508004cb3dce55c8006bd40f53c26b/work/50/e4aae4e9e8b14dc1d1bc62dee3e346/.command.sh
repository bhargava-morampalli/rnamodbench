#!/bin/bash -euo pipefail
touch test.bam.bai

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_INDEX":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
