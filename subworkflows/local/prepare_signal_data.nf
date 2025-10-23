#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: PREPARE_SIGNAL_DATA
 * Purpose: Prepare signal-level data for downstream analysis (Tombo, f5c, Yanocomp, Xpore).
 * Inputs:
 *   - mapped_fastq_16s_native: [ val(meta), path(fastq) ]
 *   - mapped_fastq_16s_ivt:    [ val(meta), path(fastq) ]
 *   - mapped_fastq_23s_native: [ val(meta), path(fastq) ]
 *   - mapped_fastq_23s_ivt:    [ val(meta), path(fastq) ]
 *   - mapped_bam_16s_native:   [ val(meta), path(bam) ]
 *   - mapped_bam_16s_ivt:      [ val(meta), path(bam) ]
 *   - mapped_bam_23s_native:   [ val(meta), path(bam) ]
 *   - mapped_bam_23s_ivt:      [ val(meta), path(bam) ]
 *   - ref_16s:                 path(reference)
 *   - ref_23s:                 path(reference)
 * Outputs:
 *   - tombo_resquiggled_16s:   [ val(meta), path(fast5) ]
 *   - tombo_resquiggled_23s:   [ val(meta), path(fast5) ]
 *   - f5c_ready_16s_native:    [ val(meta), path(fastq), path(fast5) ]
 *   - f5c_ready_16s_ivt:       [ val(meta), path(fastq), path(fast5) ]
 *   - f5c_ready_23s_native:    [ val(meta), path(fastq), path(fast5) ]
 *   - f5c_ready_23s_ivt:       [ val(meta), path(fastq), path(fast5) ]
 *   - versions:                [ path(versions.yml) ]
 */

include { EXTRACT_READ_IDS      } from '../../modules/local/extract_read_ids'
include { FAST5_SUBSET          } from '../../modules/local/fast5_subset'
include { MULTI_TO_SINGLE_FAST5 } from '../../modules/local/multi_to_single_fast5'
include { TOMBO_RESQUIGGLE      } from '../../modules/local/tombo_resquiggle'

workflow PREPARE_SIGNAL_DATA {
    take:
        mapped_fastq_16s_native  // [ val(meta), path(fastq) ]
        mapped_fastq_16s_ivt     // [ val(meta), path(fastq) ]
        mapped_fastq_23s_native  // [ val(meta), path(fastq) ]
        mapped_fastq_23s_ivt     // [ val(meta), path(fastq) ]
        mapped_bam_16s_native    // [ val(meta), path(bam) ]
        mapped_bam_16s_ivt       // [ val(meta), path(bam) ]
        mapped_bam_23s_native    // [ val(meta), path(bam) ]
        mapped_bam_23s_ivt       // [ val(meta), path(bam) ]
        ref_16s                  // path(reference)
        ref_23s                  // path(reference)

    main:
        // Combine all mapped fastq files into a single channel
        ch_mapped_fastq_all = mapped_fastq_16s_native
            .mix(mapped_fastq_16s_ivt)
            .mix(mapped_fastq_23s_native)
            .mix(mapped_fastq_23s_ivt)

        // 1. Extract read IDs from mapped FASTQ files
        EXTRACT_READ_IDS ( ch_mapped_fastq_all )

        // 2. Subset the main FAST5 directory to get relevant files
        FAST5_SUBSET ( EXTRACT_READ_IDS.out.read_ids )

        // 3. Convert multi-read FAST5 to single-read FAST5 for Tombo
        MULTI_TO_SINGLE_FAST5 ( FAST5_SUBSET.out.fast5 )

        // Collect all version files from modules
        versions_ch = Channel.empty()
        versions_ch = versions_ch
            .mix(EXTRACT_READ_IDS.out.versions)
            .mix(FAST5_SUBSET.out.versions)
            .mix(MULTI_TO_SINGLE_FAST5.out.versions)
            .collect()

        // Prepare input for Tombo resquiggle: [ meta, fastq, single_fast5_dir, reference ]
        ch_tombo_input = ch_mapped_fastq_all
            .join(MULTI_TO_SINGLE_FAST5.out.fast5)
            .map { meta, fastq, single_fast5 ->
                def reference = (meta.rrna == '16s') ? ref_16s : ref_23s
                [ meta, fastq, single_fast5, reference ]
            }

        // 4. Run Tombo resquiggle
        TOMBO_RESQUIGGLE ( ch_tombo_input )
        versions_ch = versions_ch.mix(TOMBO_RESQUIGGLE.out.versions).collect()

        // Group resquiggled files by rRNA type for modification calling
        ch_tombo_resquiggled_grouped = TOMBO_RESQUIGGLE.out.resquiggled.branch { meta, fast5 ->
            is_16s: meta.rrna == '16s'
            is_23s: meta.rrna == '23s'
        }

        // Prepare input for f5c: [ meta, fastq, multi_fast5_dir ]
        // This is for f5c index and eventalign, which can handle multi-read FAST5
        ch_f5c_ready = ch_mapped_fastq_all
            .join(FAST5_SUBSET.out.fast5)

        // Branch f5c-ready data by rRNA type and sample type
        ch_f5c_ready_branched = ch_f5c_ready.branch { meta, fastq, fast5 ->
            f5c_16s_native: meta.rrna == '16s' && meta.type == 'native'
            f5c_16s_ivt:    meta.rrna == '16s' && meta.type == 'ivt'
            f5c_23s_native: meta.rrna == '23s' && meta.type == 'native'
            f5c_23s_ivt:    meta.rrna == '23s' && meta.type == 'ivt'
        }

    emit:
        tombo_resquiggled_16s = ch_tombo_resquiggled_grouped.is_16s
        tombo_resquiggled_23s = ch_tombo_resquiggled_grouped.is_23s

        f5c_ready_16s_native = ch_f5c_ready_branched.f5c_16s_native
        f5c_ready_16s_ivt    = ch_f5c_ready_branched.f5c_16s_ivt
        f5c_ready_23s_native = ch_f5c_ready_branched.f5c_23s_native
        f5c_ready_23s_ivt    = ch_f5c_ready_branched.f5c_23s_ivt

        versions = versions_ch
}