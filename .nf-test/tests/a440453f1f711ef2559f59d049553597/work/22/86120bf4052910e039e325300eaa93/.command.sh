#!/bin/bash -euo pipefail
touch test_stub.sorted_sorted.bam

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_SORT":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
