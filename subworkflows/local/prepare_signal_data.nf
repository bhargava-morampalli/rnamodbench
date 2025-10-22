//
// Prepare signal-level data: extract FAST5 files and perform tombo resquiggle
//

include { EXTRACT_READ_IDS       } from '../../modules/local/extract_read_ids'
include { FAST5_SUBSET           } from '../../modules/local/fast5_subset'
include { MULTI_TO_SINGLE_FAST5  } from '../../modules/local/multi_to_single_fast5'
include { TOMBO_RESQUIGGLE       } from '../../modules/local/tombo_resquiggle'
include { F5C_INDEX              } from '../../modules/local/f5c_index'

workflow PREPARE_SIGNAL_DATA {
    take:
    mapped_fastq_16s_native  // channel: [ val(meta), path(fastq) ]
    mapped_fastq_16s_ivt     // channel: [ val(meta), path(fastq) ]
    mapped_fastq_23s_native  // channel: [ val(meta), path(fastq) ]
    mapped_fastq_23s_ivt     // channel: [ val(meta), path(fastq) ]
    mapped_bam_16s_native    // channel: [ val(meta), path(bam) ]
    mapped_bam_16s_ivt       // channel: [ val(meta), path(bam) ]
    mapped_bam_23s_native    // channel: [ val(meta), path(bam) ]
    mapped_bam_23s_ivt       // channel: [ val(meta), path(bam) ]
    native_fast5             // channel: path(fast5_dir)
    ivt_fast5                // channel: path(fast5_dir)
    ref_16s                  // channel: path(reference)
    ref_23s                  // channel: path(reference)

    main:
    ch_versions = Channel.empty()

    //
    // Combine all mapped FASTQs
    //
    ch_all_fastqs = mapped_fastq_16s_native
        .mix(mapped_fastq_16s_ivt)
        .mix(mapped_fastq_23s_native)
        .mix(mapped_fastq_23s_ivt)

    //
    // MODULE: Extract read IDs from mapped FASTQs
    //
    EXTRACT_READ_IDS ( ch_all_fastqs )
    ch_versions = ch_versions.mix(EXTRACT_READ_IDS.out.versions.first())

    //
    // Separate read IDs by type
    //
    EXTRACT_READ_IDS.out.read_ids
        .branch {
            native_16s: it[0].type == 'native' && it[0].rrna == '16s'
            native_23s: it[0].type == 'native' && it[0].rrna == '23s'
            ivt_16s: it[0].type == 'ivt' && it[0].rrna == '16s'
            ivt_23s: it[0].type == 'ivt' && it[0].rrna == '23s'
        }
        .set { read_ids }

    //
    // MODULE: Extract relevant FAST5 files using read IDs
    //
    FAST5_SUBSET (
        read_ids.native_16s.map { meta, ids -> [ meta, native_fast5, ids ] }
            .mix(read_ids.native_23s.map { meta, ids -> [ meta, native_fast5, ids ] })
            .mix(read_ids.ivt_16s.map { meta, ids -> [ meta, ivt_fast5, ids ] })
            .mix(read_ids.ivt_23s.map { meta, ids -> [ meta, ivt_fast5, ids ] })
    )
    ch_versions = ch_versions.mix(FAST5_SUBSET.out.versions.first())

    //
    // MODULE: Convert multi-read FAST5 to single-read FAST5
    //
    MULTI_TO_SINGLE_FAST5 ( FAST5_SUBSET.out.fast5 )
    ch_versions = ch_versions.mix(MULTI_TO_SINGLE_FAST5.out.versions.first())

    //
    // Separate single FAST5s by type and add reference
    //
    MULTI_TO_SINGLE_FAST5.out.single_fast5
        .branch {
            native_16s: it[0].type == 'native' && it[0].rrna == '16s'
            native_23s: it[0].type == 'native' && it[0].rrna == '23s'
            ivt_16s: it[0].type == 'ivt' && it[0].rrna == '16s'
            ivt_23s: it[0].type == 'ivt' && it[0].rrna == '23s'
        }
        .set { single_fast5s }

    //
    // MODULE: Tombo resquiggle for each sample
    //
    TOMBO_RESQUIGGLE (
        single_fast5s.native_16s.map { meta, fast5 -> [ meta, fast5, ref_16s ] }
            .mix(single_fast5s.native_23s.map { meta, fast5 -> [ meta, fast5, ref_23s ] })
            .mix(single_fast5s.ivt_16s.map { meta, fast5 -> [ meta, fast5, ref_16s ] })
            .mix(single_fast5s.ivt_23s.map { meta, fast5 -> [ meta, fast5, ref_23s ] })
    )
    ch_versions = ch_versions.mix(TOMBO_RESQUIGGLE.out.versions.first())

    //
    // Group tombo results by rRNA type for comparison
    //
    TOMBO_RESQUIGGLE.out.resquiggled
        .map { meta, fast5 ->
            def key = "${meta.rrna}_${meta.replicate}"
            [ key, meta.type, fast5 ]
        }
        .groupTuple()
        .set { ch_tombo_grouped }

    //
    // MODULE: F5C index for f5c eventalign
    //
    F5C_INDEX (
        FAST5_SUBSET.out.fast5.join(ch_all_fastqs, by: 0)
    )
    ch_versions = ch_versions.mix(F5C_INDEX.out.versions.first())

    //
    // Combine indexed FASTQ with corresponding BAM
    //
    ch_all_bams = mapped_bam_16s_native
        .mix(mapped_bam_16s_ivt)
        .mix(mapped_bam_23s_native)
        .mix(mapped_bam_23s_ivt)

    F5C_INDEX.out.indexed
        .join(ch_all_bams, by: 0)
        .map { meta, fastq, fastq_idx, readdb, gzi, fai, fast5, bam, bai ->
            [ meta, fastq, fastq_idx, readdb, gzi, fai, bam, bai, fast5 ]
        }
        .set { ch_f5c_ready }

    //
    // Separate f5c ready data by type
    //
    ch_f5c_ready
        .branch {
            native_16s: it[0].type == 'native' && it[0].rrna == '16s'
            native_23s: it[0].type == 'native' && it[0].rrna == '23s'
            ivt_16s: it[0].type == 'ivt' && it[0].rrna == '16s'
            ivt_23s: it[0].type == 'ivt' && it[0].rrna == '23s'
        }
        .set { f5c_ready }

    emit:
    // Tombo resquiggled FAST5s grouped by rRNA and replicate
    tombo_resquiggled_16s = ch_tombo_grouped.filter { it[0].contains('16s') }  // channel: [ key, [types], [fast5s] ]
    tombo_resquiggled_23s = ch_tombo_grouped.filter { it[0].contains('23s') }  // channel: [ key, [types], [fast5s] ]

    // F5C ready data
    f5c_ready_16s_native = f5c_ready.native_16s  // channel: [ val(meta), paths... ]
    f5c_ready_16s_ivt    = f5c_ready.ivt_16s     // channel: [ val(meta), paths... ]
    f5c_ready_23s_native = f5c_ready.native_23s  // channel: [ val(meta), paths... ]
    f5c_ready_23s_ivt    = f5c_ready.ivt_23s     // channel: [ val(meta), paths... ]

    versions = ch_versions                       // channel: [ versions.yml ]
}
