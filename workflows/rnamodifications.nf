#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * RNAMODIFICATIONS Workflow
 *
 * Dynamic RNA modification detection pipeline supporting arbitrary target rRNA types.
 * Samples are mapped directly to their designated references based on the 'target' column
 * in the samplesheet (no dual-mapping overhead).
 */

include { MAPPING_RRNA         } from '../subworkflows/local/mapping_rrna'
include { QC_STATS             } from '../subworkflows/local/qc_stats'
include { PREPARE_SIGNAL_DATA  } from '../subworkflows/local/prepare_signal_data'
include { SIGNAL_PROCESSING    } from '../subworkflows/local/signal_processing'
include { MODIFICATION_CALLING } from '../subworkflows/local/modification_calling'

// Import nf-schema functions
include { samplesheetToList } from 'plugin/nf-schema'

workflow RNAMODIFICATIONS {
    ch_versions = Channel.empty()

    // =========================================================================
    // PARSE REFERENCES CSV
    // =========================================================================
    // Format: target,reference (e.g., 16s,/path/to/16s.fa)
    // Creates channel: [ val(target), path(reference) ]

    if (params.references) {
        ch_references = Channel
            .fromPath(params.references, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def target = row.target.toLowerCase()
                def ref_file = file(row.reference, checkIfExists: true)
                [ target, ref_file ]
            }

        // Also create a reference map for downstream use
        ch_ref_map = ch_references.collect().map { ref_list ->
            def ref_map = [:]
            ref_list.each { target, ref_file ->
                ref_map[target] = ref_file
            }
            ref_map
        }
    } else if (params.ref_16s && params.ref_23s) {
        // DEPRECATED: Support legacy ref_16s/ref_23s parameters for backwards compatibility
        log.warn "Using deprecated ref_16s/ref_23s parameters. Please migrate to --references CSV format."
        ch_references = Channel.of(
            ['16s', file(params.ref_16s, checkIfExists: true)],
            ['23s', file(params.ref_23s, checkIfExists: true)]
        )
        ch_ref_map = Channel.value([
            '16s': file(params.ref_16s),
            '23s': file(params.ref_23s)
        ])
    } else {
        error "ERROR: Please provide --references parameter (CSV file mapping targets to reference FASTA files)"
    }

    // =========================================================================
    // PARSE SAMPLESHEET
    // =========================================================================
    // Uses nf-schema plugin to validate against schema
    // Creates channel: [ val(meta), path(fastq) ]
    // meta contains: id, single_end, type, replicate, fast5_dir, rrna

    reads = Channel
        .fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.single_end = true
            meta.type = row.type
            meta.replicate = row.replicate
            meta.fast5_dir = file(row.fast5_dir)
            meta.rrna = row.target.toLowerCase()  // Target rRNA type (e.g., '16s', '23s')
            [ meta, file(row.fastq) ]
        }

    // Get unique targets from samplesheet for validation
    ch_targets = reads.map { meta, fastq -> meta.rrna }.unique()

    // Validate that all targets in samplesheet have corresponding references
    ch_targets
        .combine(ch_ref_map)
        .map { target, ref_map ->
            if (!ref_map.containsKey(target)) {
                error "ERROR: Target '${target}' in samplesheet has no corresponding reference in references CSV. Available: ${ref_map.keySet()}"
            }
        }
        .collect()

    // Validate that each target+replicate has both native AND IVT samples
    // All downstream modification calling tools require native vs IVT comparison
    reads
        .map { meta, fastq ->
            def key = "${meta.rrna}_${meta.replicate}"
            [ key, meta.type ]
        }
        .groupTuple()
        .map { key, types ->
            def has_native = types.contains('native')
            def has_ivt = types.contains('ivt')
            if (!has_native || !has_ivt) {
                def missing = []
                if (!has_native) missing << 'native'
                if (!has_ivt) missing << 'ivt'
                error "ERROR: Target+replicate '${key}' is missing ${missing.join(' and ')} sample(s). " +
                      "All modification calling tools require both native and IVT samples for comparison."
            }
        }
        .collect()

    // Split reads by type for downstream FAST5 handling
    native_reads = reads.filter { it[0].type == 'native' }
    ivt_reads = reads.filter { it[0].type == 'ivt' }

    // =========================================================================
    // MAPPING
    // =========================================================================
    // Maps each sample ONCE to its designated reference (based on meta.rrna)

    MAPPING_RRNA(reads, ch_references)
    ch_versions = ch_versions.mix(MAPPING_RRNA.out.versions)

    // =========================================================================
    // QC STATS
    // =========================================================================
    // Unified BAM channel - QC_STATS processes all BAMs regardless of target

    QC_STATS(MAPPING_RRNA.out.mapped_bams)
    ch_versions = ch_versions.mix(QC_STATS.out.versions)

    // =========================================================================
    // PREPARE SIGNAL DATA
    // =========================================================================
    // Get FAST5 directories from reads
    native_fast5 = native_reads.map { it[0].fast5_dir }.first()
    ivt_fast5 = ivt_reads.map { it[0].fast5_dir }.first()

    PREPARE_SIGNAL_DATA(
        MAPPING_RRNA.out.mapped_fastqs,
        MAPPING_RRNA.out.mapped_bams,
        native_fast5,
        ivt_fast5
    )
    ch_versions = ch_versions.mix(PREPARE_SIGNAL_DATA.out.versions)

    // =========================================================================
    // SIGNAL PROCESSING
    // =========================================================================

    SIGNAL_PROCESSING(
        PREPARE_SIGNAL_DATA.out.single_fast5,
        PREPARE_SIGNAL_DATA.out.f5c_ready,
        ch_ref_map
    )
    ch_versions = ch_versions.mix(SIGNAL_PROCESSING.out.versions)

    // =========================================================================
    // MODIFICATION CALLING
    // =========================================================================

    MODIFICATION_CALLING(
        SIGNAL_PROCESSING.out.tombo_resquiggled,
        SIGNAL_PROCESSING.out.eventalign,
        SIGNAL_PROCESSING.out.eventalign_xpore,
        MAPPING_RRNA.out.mapped_bams,
        ch_ref_map
    )
    ch_versions = ch_versions.mix(MODIFICATION_CALLING.out.versions)

    // =========================================================================
    // COLLECT VERSIONS
    // =========================================================================

    ch_versions.unique().map { it.text }.collectFile(
        storeDir: "${params.outdir}/pipeline_info",
        name: 'software_versions.yml',
        sort: true,
        newLine: true
    )
}
