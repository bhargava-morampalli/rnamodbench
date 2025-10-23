#!/bin/bash -euo pipefail
fast5_subset \
       \
      --input /home/bmorampa/new_improved/modules/local/fast5_subset/tests/data/fast5_input \
      --read_id_list read_ids.txt \
      --save_path fast5_subset \
      --recursive

  echo -e '"FAST5_SUBSET":
ont-fast5-api: $(fast5_subset --version | sed -nE "s/fast5_subset, version (.*)/\1/p")' > versions.yml
