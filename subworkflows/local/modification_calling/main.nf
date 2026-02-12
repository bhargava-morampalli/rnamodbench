#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: MODIFICATION_CALLING
 * Purpose: Perform RNA modification calling using multiple tools:
 *          Tombo, Yanocomp, Nanocompore, Xpore, ELIGOS, EpiNano-Error, DiffErr, DRUMMER, JACUSA2.
 *
 * Inputs:
 *   - tombo_resquiggled: [ val(meta), path(fast5) ] - resquiggled FAST5 files
 *   - eventalign:        [ val(meta), path(eventalign) ] - for nanocompore/yanocomp (read_name column)
 *   - eventalign_xpore:  [ val(meta), path(eventalign) ] - for xpore (read_index column)
 *   - bams:              [ val(meta), path(bam), path(bai) ] - mapped BAM files
 *   - ref_map:           val(map) - { 'target': ref_file, ... }
 *
 * Outputs:
 *   - tombo_bed:           [ val(meta), path(bed) ]
 *   - yanocomp_bed:        [ val(meta), path(bed) ]
 *   - nanocompore_results: [ val(key), path(results) ]
 *   - xpore_table:         [ val(key), path(table) ]
 *   - eligos_results:      [ val(meta), path(results) ]
 *   - epinano_results:     [ val(meta), path(results) ]
 *   - differr_bed:         [ val(meta), path(bed) ]
 *   - drummer_results:     [ val(meta), path(results) ]
 *   - jacusa2_bed:         [ val(meta), path(bed) ]
 *   - versions:            [ path(versions.yml) ]
 *
 * Note: All tools group samples by ${meta.rrna}_${meta.replicate} and pair native/ivt.
 *       This works dynamically for any target type, not just 16s/23s.
 */

include { TOMBO_DETECT_MODIFICATIONS      } from '../../../modules/local/tombo_detect_modifications'
include { TOMBO_TEXT_OUTPUT               } from '../../../modules/local/tombo_text_output'
include { YANOCOMP_PREPARE                } from '../../../modules/local/yanocomp_prepare'
include { YANOCOMP_ANALYSIS               } from '../../../modules/local/yanocomp_analysis'
include { NANOCOMPORE_EVENTALIGN_COLLAPSE } from '../../../modules/local/nanocompore_eventalign_collapse'
include { NANOCOMPORE_SAMPCOMP            } from '../../../modules/local/nanocompore_sampcomp'
include { XPORE_DATAPREP                  } from '../../../modules/local/xpore_dataprep'
include { XPORE_DIFFMOD                   } from '../../../modules/local/xpore_diffmod'
include { ELIGOS_PAIR_DIFF_MOD            } from '../../../modules/local/eligos_pair_diff_mod'
include { EPINANO_ERROR                   } from '../../../modules/local/epinano_error'
include { DIFFERR                         } from '../../../modules/local/differr'
include { DRUMMER                         } from '../../../modules/local/drummer'
include { JACUSA2                         } from '../../../modules/local/jacusa2'
// Note: NANORMS is not included as it requires EpiNano per-site CSV files which are not currently generated

workflow MODIFICATION_CALLING {
    take:
        tombo_resquiggled  // [ val(meta), path(fast5) ]
        eventalign         // [ val(meta), path(eventalign) ] - for nanocompore/yanocomp
        eventalign_xpore   // [ val(meta), path(eventalign) ] - for xpore
        bams               // [ val(meta), path(bam), path(bai) ]
        ref_map            // val(map) - { 'target': ref_file, ... }

    main:
        // Collect all version info
        versions_ch = Channel.empty()

        // =====================================================================
        // TOMBO: Detect modifications by comparing native vs IVT
        // =====================================================================
        // Group by target+replicate, collect native and IVT samples
        ch_tombo_ready = tombo_resquiggled
            .map { meta, fast5 ->
                def key = "${meta.rrna}_${meta.replicate}"
                [ key, meta.type, fast5 ]
            }
            .groupTuple()

        TOMBO_DETECT_MODIFICATIONS ( ch_tombo_ready )
        versions_ch = versions_ch.mix(TOMBO_DETECT_MODIFICATIONS.out.versions)

        // TOMBO: Extract text output from stats files
        // Wire reference FASTA from ref_map using target extracted from the key
        ch_tombo_stats_with_ref = TOMBO_DETECT_MODIFICATIONS.out.stats
            .combine(ref_map)
            .map { key, stats, refs ->
                def rrna = key.split('_')[0]
                def reference = refs[rrna.toLowerCase()]
                if (!reference) {
                    error "No reference found for target '${rrna}' in ref_map. Available: ${refs.keySet()}"
                }
                [ key, stats, reference ]
            }

        TOMBO_TEXT_OUTPUT ( ch_tombo_stats_with_ref )
        versions_ch = versions_ch.mix(TOMBO_TEXT_OUTPUT.out.versions)

        // =====================================================================
        // YANOCOMP: Prepare HDF5 files from f5c eventalign
        // =====================================================================
        ch_no_summary = file('NO_FILE')  // Placeholder for optional summary file
        YANOCOMP_PREPARE ( eventalign, ch_no_summary )
        versions_ch = versions_ch.mix(YANOCOMP_PREPARE.out.versions)

        // Join native and IVT HDF5 files for comparison
        ch_yanocomp_ready = YANOCOMP_PREPARE.out.hdf5
            .map { meta, hdf5 ->
                def key = "${meta.rrna}_${meta.replicate}"
                [ key, meta.type, hdf5 ]
            }
            .groupTuple()
            .filter { key, types, hdf5s ->
                def has_native = types.contains('native')
                def has_ivt = types.contains('ivt')
                if (!has_native || !has_ivt) {
                    log.warn "Skipping yanocomp for ${key}: missing ${!has_native ? 'native' : 'ivt'} sample"
                }
                has_native && has_ivt
            }
            .map { key, types, hdf5s ->
                def meta = [:]
                meta.id = key
                meta.rrna = key.split('_')[0]
                meta.replicate = key.split('_').drop(1).join('_')  // Handle replicate names with underscores
                def native_idx = types.findIndexOf { it == 'native' }
                def ivt_idx = types.findIndexOf { it == 'ivt' }
                [ meta, hdf5s[native_idx], hdf5s[ivt_idx] ]
            }

        YANOCOMP_ANALYSIS ( ch_yanocomp_ready )
        versions_ch = versions_ch.mix(YANOCOMP_ANALYSIS.out.versions)

        // =====================================================================
        // NANOCOMPORE: Collapse eventalign files
        // =====================================================================
        NANOCOMPORE_EVENTALIGN_COLLAPSE ( eventalign )
        versions_ch = versions_ch.mix(NANOCOMPORE_EVENTALIGN_COLLAPSE.out.versions)

        // Group and pair native/IVT for nanocompore
        ch_nanocompore_grouped = NANOCOMPORE_EVENTALIGN_COLLAPSE.out.collapsed
            .map { meta, collapsed_dir ->
                def key = "${meta.rrna}_${meta.replicate}"
                [ key, meta.type, collapsed_dir, meta.rrna ]
            }
            .groupTuple()
            .filter { key, types, collapsed_dirs, rrnas ->
                def has_native = types.contains('native')
                def has_ivt = types.contains('ivt')
                if (!has_native || !has_ivt) {
                    log.warn "Skipping nanocompore for ${key}: missing ${!has_native ? 'native' : 'ivt'} sample"
                }
                has_native && has_ivt
            }
            .map { key, types, collapsed_dirs, rrnas ->
                def rrna = rrnas[0]
                def native_idx = types.findIndexOf { it == 'native' }
                def ivt_idx = types.findIndexOf { it == 'ivt' }
                def native_dir = collapsed_dirs[native_idx]
                def ivt_dir = collapsed_dirs[ivt_idx]
                [ key, 'native', 'ivt', native_dir, ivt_dir, rrna ]
            }

        // Combine with appropriate reference based on rRNA type
        ch_nanocompore_ready = ch_nanocompore_grouped
            .combine(ref_map)
            .map { key, native_label, ivt_label, native_dir, ivt_dir, rrna, refs ->
                def reference = refs[rrna.toLowerCase()]
                if (!reference) {
                    error "No reference found for target '${rrna}'. Available: ${refs.keySet()}"
                }
                [ key, native_label, ivt_label, native_dir, ivt_dir, reference ]
            }

        NANOCOMPORE_SAMPCOMP ( ch_nanocompore_ready )
        versions_ch = versions_ch.mix(NANOCOMPORE_SAMPCOMP.out.versions)

        // =====================================================================
        // XPORE: Data preparation using xpore-compatible eventalign
        // =====================================================================
        XPORE_DATAPREP ( eventalign_xpore )
        versions_ch = versions_ch.mix(XPORE_DATAPREP.out.versions)

        // Join native and IVT xpore data for comparison
        ch_xpore_ready = XPORE_DATAPREP.out.dataprep
            .map { meta, dataprep_dir ->
                def key = "${meta.rrna}_${meta.replicate}"
                [ key, meta.type, dataprep_dir ]
            }
            .groupTuple()

        XPORE_DIFFMOD ( ch_xpore_ready )
        versions_ch = versions_ch.mix(XPORE_DIFFMOD.out.versions)

        // =====================================================================
        // BAM-BASED TOOLS: ELIGOS, EPINANO, DIFFERR, DRUMMER, JACUSA2
        // =====================================================================
        // Group BAMs by target+replicate, pair native/IVT, combine with reference
        ch_bam_paired = bams
            .map { meta, bam, bai ->
                def key = "${meta.rrna}_${meta.replicate}"
                [ key, meta, bam, bai, meta.type ]
            }
            .groupTuple()
            .filter { key, metas, bams_list, bais_list, types ->
                def has_native = types.contains('native')
                def has_ivt = types.contains('ivt')
                if (!has_native || !has_ivt) {
                    log.warn "Skipping BAM-based tools for ${key}: missing ${!has_native ? 'native' : 'ivt'} sample"
                }
                has_native && has_ivt
            }
            .map { key, metas, bams_list, bais_list, types ->
                def native_idx = types.findIndexOf { it == 'native' }
                def ivt_idx = types.findIndexOf { it == 'ivt' }
                def meta = [:]
                meta.id = key
                meta.rrna = metas[0].rrna
                meta.replicate = metas[0].replicate
                [ meta, bams_list[native_idx], bais_list[native_idx], bams_list[ivt_idx], bais_list[ivt_idx] ]
            }
            .combine(ref_map)
            .map { meta, native_bam, native_bai, ivt_bam, ivt_bai, refs ->
                def reference = refs[meta.rrna.toLowerCase()]
                if (!reference) {
                    error "No reference found for target '${meta.rrna}'. Available: ${refs.keySet()}"
                }
                [ meta, native_bam, native_bai, ivt_bam, ivt_bai, reference ]
            }

        // Run all BAM-based tools
        ELIGOS_PAIR_DIFF_MOD ( ch_bam_paired )
        versions_ch = versions_ch.mix(ELIGOS_PAIR_DIFF_MOD.out.versions)

        EPINANO_ERROR ( ch_bam_paired )
        versions_ch = versions_ch.mix(EPINANO_ERROR.out.versions)

        DIFFERR ( ch_bam_paired )
        versions_ch = versions_ch.mix(DIFFERR.out.versions)

        DRUMMER ( ch_bam_paired )
        versions_ch = versions_ch.mix(DRUMMER.out.versions)

        JACUSA2 ( ch_bam_paired )
        versions_ch = versions_ch.mix(JACUSA2.out.versions)

        // NANORMS: Currently disabled as it requires EpiNano per-site CSV files
        ch_nanorms_results = Channel.empty()

    emit:
        tombo_bed           = TOMBO_TEXT_OUTPUT.out.bed
        yanocomp_bed        = YANOCOMP_ANALYSIS.out.bed
        nanocompore_results = NANOCOMPORE_SAMPCOMP.out.results
        xpore_table         = XPORE_DIFFMOD.out.diffmod
        eligos_results      = ELIGOS_PAIR_DIFF_MOD.out.results
        epinano_results     = EPINANO_ERROR.out.results
        differr_bed         = DIFFERR.out.bed
        drummer_results     = DRUMMER.out.results
        jacusa2_bed         = JACUSA2.out.bed
        nanorms_results     = ch_nanorms_results
        versions            = versions_ch.collect()
}
