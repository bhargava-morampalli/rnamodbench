#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * RNAMODBENCH Workflow
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
include { DOWNSTREAM_ANALYSIS  } from '../modules/local/downstream_analysis'

workflow RNAMODBENCH {
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
                def target = row instanceof Map ? row.target : row[0]
                def reference = row instanceof Map ? row.reference : row[1]
                def ref_file = file(reference, checkIfExists: true)
                target = target.toString().toLowerCase()
                [ target, ref_file ]
            }

        // Also create a reference map for downstream use
        ch_ref_map = ch_references.collect().map { ref_list ->
            def ref_map = [:]
            if (ref_list && ref_list[0] instanceof List) {
                ref_list.each { item ->
                    def target = item[0]
                    def ref_file = item[1]
                    ref_map[target.toString()] = ref_file
                }
            } else {
                ref_list.collate(2).each { item ->
                    if (item.size() == 2) {
                        def target = item[0]
                        def ref_file = item[1]
                        ref_map[target.toString()] = ref_file
                    }
                }
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
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def target = row.target ?: row.rrna
            if (!target) {
                error "ERROR: Samplesheet row for sample '${row.sample}' is missing required target column."
            }

            def meta = [:]
            meta.id = row.sample
            meta.single_end = true
            meta.type = row.type
            meta.replicate = row.replicate
            meta.fast5_dir = file(row.fast5_dir, checkIfExists: true)
            meta.rrna = target.toString().toLowerCase()  // Target rRNA type (e.g., '16s', '23s')
            [ meta, file(row.fastq, checkIfExists: true) ]
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
    // Pooled FAST5 mode: require exactly one canonical FAST5 directory per type.
    // Canonicalization uses real paths so symlink/path variants collapse.
    native_fast5 = native_reads
        .map { meta, fastq ->
            try {
                meta.fast5_dir.toRealPath().toString()
            } catch (Exception e) {
                error "ERROR: Could not resolve FAST5 directory for sample '${meta.id}' (${meta.type}): ${meta.fast5_dir}"
            }
        }
        .unique()
        .collect()
        .map { dirs ->
            def sorted_dirs = dirs.collect { it }.sort()
            if (sorted_dirs.size() != 1) {
                error "ERROR: Pooled FAST5 mode requires exactly one FAST5 directory for type 'native'. Found ${sorted_dirs.size()}: ${sorted_dirs.join(', ')}"
            }
            file(sorted_dirs[0], checkIfExists: true)
        }

    ivt_fast5 = ivt_reads
        .map { meta, fastq ->
            try {
                meta.fast5_dir.toRealPath().toString()
            } catch (Exception e) {
                error "ERROR: Could not resolve FAST5 directory for sample '${meta.id}' (${meta.type}): ${meta.fast5_dir}"
            }
        }
        .unique()
        .collect()
        .map { dirs ->
            def sorted_dirs = dirs.collect { it }.sort()
            if (sorted_dirs.size() != 1) {
                error "ERROR: Pooled FAST5 mode requires exactly one FAST5 directory for type 'ivt'. Found ${sorted_dirs.size()}: ${sorted_dirs.join(', ')}"
            }
            file(sorted_dirs[0], checkIfExists: true)
        }

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
    // OPTIONAL DOWNSTREAM ANALYSIS (default: disabled)
    // =========================================================================

    if (params.run_downstream) {
        def gt_input = params.ground_truth ? file(params.ground_truth, checkIfExists: true).toString() : 'NO_FILE'
        def refs_input = params.references ? file(params.references, checkIfExists: true).toString() : 'NO_FILE'

        ch_modifications_dir = MODIFICATION_CALLING.out.versions
            .map { _ -> file("${params.outdir}/modifications") }

        DOWNSTREAM_ANALYSIS(
            ch_modifications_dir,
            Channel.value(gt_input),
            Channel.value(refs_input)
        )
        ch_versions = ch_versions.mix(DOWNSTREAM_ANALYSIS.out.versions)
    }

    // =========================================================================
    // COLLECT VERSIONS
    // =========================================================================

    ch_versions
        // Subworkflows emit either a single versions.yml path or a collected list of paths.
        .flatMap { item -> item instanceof Collection ? item : [item] }
        .map { version_file ->
            version_file
                .text
                .readLines()
                .findAll { line -> line.trim() != 'END_VERSIONS' }
                .join('\n')
                .trim()
        }
        .filter { entry -> entry }
        .unique()
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true
        )
}
