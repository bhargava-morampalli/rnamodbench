#!/bin/bash -euo pipefail
touch test_stub.read_ids.txt

cat <<-END_VERSIONS > versions.yml
"EXTRACT_READ_IDS":
    seqkit: $(seqkit version | sed 's/seqkit version //')
END_VERSIONS
