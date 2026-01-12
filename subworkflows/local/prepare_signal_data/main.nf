#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: PREPARE_SIGNAL_DATA
 * Purpose: Prepare signal-level data for downstream analysis (Tombo, f5c, Yanocomp, Xpore).
 *
 * Inputs:
 *   - mapped_fastqs:   [ val(meta), path(fastq) ] - unified channel for all mapped FASTQs
 *   - mapped_bams:     [ val(meta), path(bam), path(bai) ] - unified channel for all mapped BAMs
 *   - native_fast5_dir: path(directory) - FAST5 directory for native samples
 *   - ivt_fast5_dir:    path(directory) - FAST5 directory for IVT samples
 *
 * Outputs:
 *   - single_fast5: [ val(meta), path(fast5) ] - single-read FAST5 files for Tombo
 *   - f5c_ready:    [ val(meta), path(fastq), path(bam), path(bai), path(fast5) ] - for f5c processing
 *   - versions:     [ path(versions.yml) ]
 *
 * Note: Data is processed dynamically based on meta.rrna and meta.type.
 *       No hardcoded target types - supports any rRNA target.
 */

include { EXTRACT_READ_IDS      } from '../../../modules/local/extract_read_ids'
include { FAST5_SUBSET          } from '../../../modules/local/fast5_subset'
include { MULTI_TO_SINGLE_FAST5 } from '../../../modules/local/multi_to_single_fast5'

workflow PREPARE_SIGNAL_DATA {
    take:
        mapped_fastqs     // [ val(meta), path(fastq) ] - unified channel
        mapped_bams       // [ val(meta), path(bam), path(bai) ] - unified channel
        native_fast5_dir  // path(directory)
        ivt_fast5_dir     // path(directory)

    main:
        // Join mapped FASTQs with BAMs by sample ID
        // Result: [ meta, fastq, bam, bai ]
        ch_mapped_all = mapped_fastqs.join(mapped_bams)

        // 1. Extract read IDs from mapped FASTQ files
        ch_fastqs_for_ids = ch_mapped_all.map { meta, fastq, bam, bai -> [ meta, fastq ] }
        EXTRACT_READ_IDS ( ch_fastqs_for_ids )

        // 2. Subset the FAST5 directory to get relevant files
        // Select correct FAST5 directory based on sample type (native/ivt)
        // Wrap directory paths in channels for combine operation
        ch_native_fast5 = Channel.value(native_fast5_dir)
        ch_ivt_fast5 = Channel.value(ivt_fast5_dir)

        ch_fast5_subset_input = EXTRACT_READ_IDS.out.read_ids
            .combine(ch_native_fast5)
            .combine(ch_ivt_fast5)
            .map { meta, read_ids, native_dir, ivt_dir ->
                def fast5_dir = (meta.type == 'native') ? native_dir : ivt_dir
                [ meta, read_ids, fast5_dir ]
            }

        FAST5_SUBSET ( ch_fast5_subset_input )

        // 3. Convert multi-read FAST5 to single-read FAST5 for Tombo
        MULTI_TO_SINGLE_FAST5 ( FAST5_SUBSET.out.fast5 )

        // Collect all version files from modules
        versions_ch = Channel.empty()
        versions_ch = versions_ch
            .mix(EXTRACT_READ_IDS.out.versions)
            .mix(FAST5_SUBSET.out.versions)
            .mix(MULTI_TO_SINGLE_FAST5.out.versions)
            .collect()

        // Prepare input for f5c: [ meta, fastq, bam, bai, fast5 ]
        // Join original inputs (fastq, bam, bai) with subsetted fast5s
        ch_f5c_ready = ch_mapped_all
            .join(FAST5_SUBSET.out.fast5)
            .map { meta, fastq, bam, bai, fast5 ->
                [ meta, fastq, bam, bai, fast5 ]
            }

    emit:
        single_fast5 = MULTI_TO_SINGLE_FAST5.out.fast5  // [ val(meta), path(fast5) ]
        f5c_ready    = ch_f5c_ready                      // [ val(meta), path(fastq), path(bam), path(bai), path(fast5) ]
        versions     = versions_ch
}
