#!/bin/bash -euo pipefail
touch test_stub.txt

cat <<-END_VERSIONS > versions.yml
"SAMTOOLS_DEPTH":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
