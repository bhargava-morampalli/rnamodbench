#!/bin/bash -euo pipefail
seqkit seq -n -i test_nanopore.fastq.gz > test_nanopore_gz.read_ids.txt

cat <<-END_VERSIONS > versions.yml
"EXTRACT_READ_IDS":
    seqkit: $(seqkit version | sed 's/seqkit version //')
END_VERSIONS
