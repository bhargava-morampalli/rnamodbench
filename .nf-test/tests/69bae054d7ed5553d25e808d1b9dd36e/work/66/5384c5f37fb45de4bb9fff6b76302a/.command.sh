#!/bin/bash -euo pipefail
touch test_stub.bam

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_VIEW":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
