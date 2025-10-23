#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: signal_processing
 * Purpose: Run f5c eventalign on signal-level data for all sample types and rRNA references.
 * Inputs:
 *   - f5c_ready_16s_native: [ val(meta), path(fastq), path(idx), path(readdb), path(gzi), path(fai), path(bam), path(bai), path(fast5) ]
 *   - f5c_ready_16s_ivt:    [ val(meta), path(fastq), path(idx), path(readdb), path(gzi), path(fai), path(bam), path(bai), path(fast5) ]
 *   - f5c_ready_23s_native: [ val(meta), path(fastq), path(idx), path(readdb), path(gzi), path(fai), path(bam), path(bai), path(fast5) ]
 *   - f5c_ready_23s_ivt:    [ val(meta), path(fastq), path(idx), path(readdb), path(gzi), path(fai), path(bam), path(bai), path(fast5) ]
 *   - ref_16s:              path(reference)
 *   - ref_23s:              path(reference)
 * Outputs:
 *   - eventalign_16s_native: [ val(meta), path(eventalign) ]
 *   - eventalign_16s_ivt:    [ val(meta), path(eventalign) ]
 *   - eventalign_23s_native: [ val(meta), path(eventalign) ]
 *   - eventalign_23s_ivt:    [ val(meta), path(eventalign) ]
 *   - versions:              [ path(versions.yml) ]
 */

include { F5C_EVENTALIGN } from '../../modules/local/f5c_eventalign'

workflow signal_processing {
    take:
        f5c_ready_16s_native  // [ val(meta), path(fastq), ... ]
        f5c_ready_16s_ivt     // [ val(meta), path(fastq), ... ]
        f5c_ready_23s_native  // [ val(meta), path(fastq), ... ]
        f5c_ready_23s_ivt     // [ val(meta), path(fastq), ... ]
        ref_16s               // path(reference)
        ref_23s               // path(reference)

    main:
        // Prepare eventalign input for each rRNA/sample type
        ch_eventalign_input = f5c_ready_16s_native.map { meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5 ->
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

        // Run f5c eventalign
        F5C_EVENTALIGN ( ch_eventalign_input )

        // Collect version info
        versions_ch = Channel.empty().mix(F5C_EVENTALIGN.out.versions)

        // Branch eventalign outputs by sample type and rRNA
        F5C_EVENTALIGN.out.eventalign
            .branch {
                native_16s: it[0].type == 'native' && it[0].rrna == '16s'
                ivt_16s:    it[0].type == 'ivt'    && it[0].rrna == '16s'
                native_23s: it[0].type == 'native' && it[0].rrna == '23s'
                ivt_23s:    it[0].type == 'ivt'    && it[0].rrna == '23s'
            }
            .set { eventalign_results }

    emit:
        eventalign_16s_native = eventalign_results.native_16s  // [ val(meta), path(eventalign) ]
        eventalign_16s_ivt    = eventalign_results.ivt_16s     // [ val(meta), path(eventalign) ]
        eventalign_23s_native = eventalign_results.native_23s  // [ val(meta), path(eventalign) ]
        eventalign_23s_ivt    = eventalign_results.ivt_23s     // [ val(meta), path(eventalign) ]
        versions              = versions_ch                    // [ path(versions.yml) ]
}
