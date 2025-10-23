#!/bin/bash -euo pipefail
touch test_stub.sam

cat <<-END_VERSIONS > versions.yml
"MINIMAP2_ALIGN":
    minimap2: $(minimap2 --version 2>&1)
END_VERSIONS
