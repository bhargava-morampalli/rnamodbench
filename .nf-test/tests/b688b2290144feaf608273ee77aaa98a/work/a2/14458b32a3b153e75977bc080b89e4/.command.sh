#!/bin/bash -euo pipefail
seqkit seq -n -i empty.fastq > empty_test.read_ids.txt

cat <<-END_VERSIONS > versions.yml
"EXTRACT_READ_IDS":
    seqkit: $(seqkit version | sed 's/seqkit version //')
END_VERSIONS
