//
// Signal processing: f5c eventalign
//

include { F5C_EVENTALIGN } from '../../modules/local/f5c_eventalign'

workflow SIGNAL_PROCESSING {
    take:
    f5c_ready_16s_native  // channel: [ val(meta), paths... ]
    f5c_ready_16s_ivt     // channel: [ val(meta), paths... ]
    f5c_ready_23s_native  // channel: [ val(meta), paths... ]
    f5c_ready_23s_ivt     // channel: [ val(meta), paths... ]
    ref_16s               // channel: path(reference)
    ref_23s               // channel: path(reference)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: F5C eventalign
    //
    F5C_EVENTALIGN (
        f5c_ready_16s_native.map { meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5 ->
            [ meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5, ref_16s ]
        }
        .mix(f5c_ready_16s_ivt.map { meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5 ->
            [ meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5, ref_16s ]
        })
        .mix(f5c_ready_23s_native.map { meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5 ->
            [ meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5, ref_23s ]
        })
        .mix(f5c_ready_23s_ivt.map { meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5 ->
            [ meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5, ref_23s ]
        })
    )
    ch_versions = ch_versions.mix(F5C_EVENTALIGN.out.versions.first())

    //
    // Separate eventalign outputs by type
    //
    F5C_EVENTALIGN.out.eventalign
        .branch {
            native_16s: it[0].type == 'native' && it[0].rrna == '16s'
            native_23s: it[0].type == 'native' && it[0].rrna == '23s'
            ivt_16s: it[0].type == 'ivt' && it[0].rrna == '16s'
            ivt_23s: it[0].type == 'ivt' && it[0].rrna == '23s'
        }
        .set { eventalign_results }

    emit:
    eventalign_16s_native = eventalign_results.native_16s  // channel: [ val(meta), path(eventalign) ]
    eventalign_16s_ivt    = eventalign_results.ivt_16s     // channel: [ val(meta), path(eventalign) ]
    eventalign_23s_native = eventalign_results.native_23s  // channel: [ val(meta), path(eventalign) ]
    eventalign_23s_ivt    = eventalign_results.ivt_23s     // channel: [ val(meta), path(eventalign) ]
    versions              = ch_versions                    // channel: [ versions.yml ]
}
