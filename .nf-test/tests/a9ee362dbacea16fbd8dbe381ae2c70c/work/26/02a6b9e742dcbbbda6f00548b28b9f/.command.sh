#!/bin/bash -euo pipefail
mkdir -p single_fast5

  echo -e '"MULTI_TO_SINGLE_FAST5":
ont-fast5-api: $(multi_to_single_fast5 --version | sed -nE "s/multi_to_single_fast5, version (.*)/\1/p")' > versions.yml
