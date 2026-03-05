#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FAST5_SUBSET } from './modules/local/fast5_subset/main.nf'

workflow {
    
    // Native inputs
    ch_native_ids = Channel.fromPath("${projectDir}/tests/data/ids_native.txt")
    ch_native_dir = Channel.fromPath('/absolute/path/to/k12_native_fast5', type: 'dir')
    ch_native_meta = Channel.value([ id: 'native_rep1' ])

    // IVT inputs
    ch_ivt_ids = Channel.fromPath("${projectDir}/tests/data/ids_ivt.txt")
    ch_ivt_dir = Channel.fromPath('/absolute/path/to/k12_ivt_fast5', type: 'dir')
    ch_ivt_meta = Channel.value([ id: 'ivt_rep1' ])

    // Run subsetting
    FAST5_SUBSET_NATIVE(ch_native_meta, ch_native_ids, ch_native_dir)
    FAST5_SUBSET_IVT(ch_ivt_meta, ch_ivt_ids, ch_ivt_dir)
}

process FAST5_SUBSET_NATIVE {
    tag "$meta.id"
    publishDir "${projectDir}/tests/data/fast5", mode: 'copy'
    
    conda "bioconda::ont-fast5-api=4.1.0"
    container "quay.io/biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0"

    input:
    val meta
    path read_ids
    path fast5_dir

    output:
    path "native_rep1"

    script:
    """
    fast5_subset \\
        --input ${fast5_dir} \\
        --read_id_list ${read_ids} \\
        --save_path native_rep1 \\
        --recursive
    """
}

process FAST5_SUBSET_IVT {
    tag "$meta.id"
    publishDir "${projectDir}/tests/data/fast5", mode: 'copy'
    
    conda "bioconda::ont-fast5-api=4.1.0"
    container "quay.io/biocontainers/ont-fast5-api:4.1.0--pyhdfd78af_0"

    input:
    val meta
    path read_ids
    path fast5_dir

    output:
    path "ivt_rep1"

    script:
    """
    fast5_subset \\
        --input ${fast5_dir} \\
        --read_id_list ${read_ids} \\
        --save_path ivt_rep1 \\
        --recursive
    """
}
