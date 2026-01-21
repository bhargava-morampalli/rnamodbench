#!/bin/bash
set -euo pipefail

# Your current approach (splice mode)
minimap2 -ax splice -uf -k14 --secondary=no MG1655.fa dRNA_k12_7aug_bc_1.fastq > alignment_test/splice_uf_k14/dRNA_k12_7aug_bc_1_spliced.sam
minimap2 -ax splice -uf -k14 --secondary=no MG1655.fa dRNA_k12_7aug_bc_2.fastq > alignment_test/splice_uf_k14/dRNA_k12_7aug_bc_2_spliced.sam
minimap2 -ax splice -uf -k14 --secondary=no MG1655.fa dRNA_k12_7aug_bc_3.fastq > alignment_test/splice_uf_k14/dRNA_k12_7aug_bc_3_spliced.sam

# Bacterial-standard approach (map-ont)  
minimap2 -ax map-ont --secondary=no MG1655.fa dRNA_k12_7aug_bc_1.fastq > alignment_test/map_ont_default/dRNA_k12_7aug_bc_1_map_ont_default.sam
minimap2 -ax map-ont --secondary=no MG1655.fa dRNA_k12_7aug_bc_2.fastq > alignment_test/map_ont_default/dRNA_k12_7aug_bc_2_map_ont_default.sam
minimap2 -ax map-ont --secondary=no MG1655.fa dRNA_k12_7aug_bc_3.fastq > alignment_test/map_ont_default/dRNA_k12_7aug_bc_3_map_ont_default.sam


# bacterial forward strand only
minimap2 -ax map-ont --for-only -k14 --secondary=no MG1655.fa dRNA_k12_7aug_bc_1.fastq > alignment_test/map_ont_foronly_k14/dRNA_k12_7aug_bc_1_map_ont_foronly_k14.sam
minimap2 -ax map-ont --for-only -k14 --secondary=no MG1655.fa dRNA_k12_7aug_bc_2.fastq > alignment_test/map_ont_foronly_k14/dRNA_k12_7aug_bc_2_map_ont_foronly_k14.sam
minimap2 -ax map-ont --for-only -k14 --secondary=no MG1655.fa dRNA_k12_7aug_bc_3.fastq > alignment_test/map_ont_foronly_k14/dRNA_k12_7aug_bc_3_map_ont_foronly_k14.sam

# For your extended rRNA references (16S and 23S)
# Splice mode
minimap2 -ax splice -uf -k14 --secondary=no k12_16S.fa dRNA_k12_7aug_bc_1.fastq > alignment_comparison/splice_uf_k14/16S.sam
minimap2 -ax splice -uf -k14 --secondary=no k12_23S.fa dRNA_k12_7aug_bc_1.fastq > alignment_comparison/splice_uf_k14/23S.sam

# Map-ont mode
minimap2 -ax map-ont -uf --secondary=no k12_16S.fa dRNA_k12_7aug_bc_1.fastq > alignment_comparison/map_ont/16S.sam
minimap2 -ax map-ont -uf --secondary=no k12_23S.fa dRNA_k12_7aug_bc_1.fastq > alignment_comparison/map_ont/23S.sam