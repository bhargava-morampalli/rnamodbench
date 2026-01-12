#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: SIGNAL_PROCESSING
 * Purpose: Run f5c eventalign and Tombo resquiggle on signal-level data.
 *
 * Inputs:
 *   - single_fast5: [ val(meta), path(fast5) ] - single-read FAST5 files
 *   - f5c_ready:    [ val(meta), path(fastq), path(bam), path(bai), path(fast5) ] - for f5c
 *   - ref_map:      val(map) - map of target -> reference file (e.g., ['16s': /path/to/16s.fa])
 *
 * Outputs:
 *   - eventalign:       [ val(meta), path(eventalign) ] - for nanocompore/yanocomp (read_name column)
 *   - eventalign_xpore: [ val(meta), path(eventalign) ] - for xpore (read_index column)
 *   - tombo_resquiggled: [ val(meta), path(fast5) ] - resquiggled FAST5 files
 *   - versions:         [ path(versions.yml) ]
 *
 * Note: Data is processed dynamically based on meta.rrna.
 *       Reference selection is done using the ref_map parameter.
 */

include { F5C_EVENTALIGN       } from '../../../modules/local/f5c_eventalign'
include { F5C_EVENTALIGN_XPORE } from '../../../modules/local/f5c_eventalign_xpore'
include { F5C_INDEX            } from '../../../modules/local/f5c_index'
include { TOMBO_RESQUIGGLE     } from '../../../modules/local/tombo_resquiggle'

workflow SIGNAL_PROCESSING {
    take:
        single_fast5  // [ val(meta), path(fast5) ]
        f5c_ready     // [ val(meta), path(fastq), path(bam), path(bai), path(fast5) ]
        ref_map       // val(map) - { '16s': ref_file, '23s': ref_file, ... }

    main:
        // Combine f5c_ready data with appropriate reference based on meta.rrna
        ch_inputs = f5c_ready
            .combine(ref_map)
            .map { meta, fastq, bam, bai, fast5, refs ->
                def target = meta.rrna.toLowerCase()
                def reference = refs[target]
                if (!reference) {
                    error "No reference found for target '${target}'. Available: ${refs.keySet()}"
                }
                [ meta, fastq, bam, bai, fast5, reference ]
            }

        // F5C Index
        // Input: [ meta, fast5_dir, fastq ], sequencing_summary (optional)
        ch_index_input = ch_inputs.map { meta, fastq, bam, bai, fast5, ref -> [ meta, fast5, fastq ] }
        ch_no_summary = file('NO_FILE')  // Placeholder for optional sequencing_summary
        F5C_INDEX ( ch_index_input, ch_no_summary )

        // F5C Eventalign
        // Input: [ meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5, ref ]
        ch_eventalign_input = F5C_INDEX.out.indexed
            .join(ch_inputs.map { meta, fastq, bam, bai, fast5, ref -> [ meta, bam, bai, ref ] })
            .map { meta, fastq, idx, readdb, gzi, fai, fast5, bam, bai, ref ->
                [ meta, fastq, idx, readdb, gzi, fai, bam, bai, fast5, ref ]
            }

        F5C_EVENTALIGN ( ch_eventalign_input )

        // F5C Eventalign for xpore (without --print-read-names, produces read_index column)
        F5C_EVENTALIGN_XPORE ( ch_eventalign_input )

        // Tombo Processing
        // Input: [ meta, single_fast5, reference ]
        ch_tombo_input = single_fast5
            .combine(ref_map)
            .map { meta, fast5, refs ->
                def target = meta.rrna.toLowerCase()
                def reference = refs[target]
                if (!reference) {
                    error "No reference found for target '${target}'. Available: ${refs.keySet()}"
                }
                [ meta, fast5, reference ]
            }

        TOMBO_RESQUIGGLE ( ch_tombo_input )

        // Collect version info
        versions_ch = Channel.empty()
            .mix(F5C_EVENTALIGN.out.versions)
            .mix(F5C_EVENTALIGN_XPORE.out.versions)
            .mix(F5C_INDEX.out.versions)
            .mix(TOMBO_RESQUIGGLE.out.versions)
            .collect()

    emit:
        // Unified eventalign outputs (consumers can filter by meta.rrna and meta.type)
        eventalign       = F5C_EVENTALIGN.out.eventalign       // [ val(meta), path(eventalign) ]
        eventalign_xpore = F5C_EVENTALIGN_XPORE.out.eventalign // [ val(meta), path(eventalign) ]
        tombo_resquiggled = TOMBO_RESQUIGGLE.out.resquiggled   // [ val(meta), path(fast5) ]
        versions         = versions_ch                          // [ path(versions.yml) ]
}
