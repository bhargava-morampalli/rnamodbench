#!/bin/bash -euo pipefail
mkdir -p fast5_subset

  echo -e '"FAST5_SUBSET":
ont-fast5-api: $(fast5_subset --version | sed -nE "s/fast5_subset, version (.*)/\1/p")' > versions.yml
