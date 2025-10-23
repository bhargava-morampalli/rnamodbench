#!/bin/bash -euo pipefail
minimap2 \
    -ax splice -uf -k14 --secondary=no \
    -t 12 \
    reference.fasta \
    test.fastq.gz \
    > test.sam

cat <<-END_VERSIONS > versions.yml
"MINIMAP2_ALIGN":
    minimap2: $(minimap2 --version 2>&1)
END_VERSIONS
