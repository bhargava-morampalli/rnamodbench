#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: MAPPING_RRNA
 * Purpose: Map reads to their designated rRNA references based on target specified in samplesheet.
 *          Each sample is mapped ONCE to its designated reference (no dual-mapping).
 *
 * Inputs:
 *   - reads:      [ val(meta), path(fastq) ] - meta must contain: id, type, rrna, replicate, fast5_dir
 *   - references: [ val(target), path(reference) ] - mapping of target names to reference files
 *
 * Outputs:
 *   - mapped_bams:   [ val(meta), path(bam), path(bai) ] - sorted, indexed BAM files
 *   - mapped_fastqs: [ val(meta), path(fastq) ] - extracted mapped reads as FASTQ
 *   - versions:      [ path(versions.yml) ]
 *
 * Note: Downstream branching should be done by consumers using meta.rrna and meta.type
 */

include { MINIMAP2_ALIGN     } from '../../../modules/local/minimap2_align'
include { SAMTOOLS_VIEW      } from '../../../modules/local/samtools_view'
include { SAMTOOLS_SORT      } from '../../../modules/local/samtools_sort'
include { SAMTOOLS_INDEX     } from '../../../modules/local/samtools_index'
include { EXTRACT_MAPPED_READS } from '../../../modules/local/extract_mapped_reads'

workflow MAPPING_RRNA {
    take:
        reads       // [ val(meta), path(fastq) ] - meta.rrna contains target (e.g., '16s', '23s')
        references  // [ val(target), path(reference) ] - e.g., ['16s', /path/to/16s.fa]

    main:
        // Normalize references channel for joining
        // Input: [ val(target), path(reference) ] -> [ val(target_lowercase), path(reference) ]
        ch_refs_normalized = references.map { target, ref_file ->
            [ target.toLowerCase(), ref_file ]
        }

        // Transform reads to have target as join key
        // Input: [ val(meta), path(fastq) ] -> [ val(target), val(meta), path(fastq) ]
        ch_reads_keyed = reads.map { meta, fastq ->
            [ meta.rrna.toLowerCase(), meta, fastq ]
        }

        // Join reads with their matching reference on target key
        // Result: [ target, meta, fastq, ref_file ]
        ch_reads_with_ref = ch_reads_keyed
            .combine(ch_refs_normalized)
            .filter { read_target, meta, fastq, ref_target, ref_file ->
                read_target == ref_target
            }
            .map { read_target, meta, fastq, ref_target, ref_file ->
                [ meta, fastq, ref_file ]
            }

        // Map reads to their designated reference using minimap2
        // Input: [ meta, fastq, reference ]
        MINIMAP2_ALIGN (
            ch_reads_with_ref.map { meta, fastq, ref -> [ meta, fastq ] },
            ch_reads_with_ref.map { meta, fastq, ref -> ref }
        )

        // Convert SAM to BAM
        SAMTOOLS_VIEW ( MINIMAP2_ALIGN.out.sam )

        // Sort BAM files
        SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )

        // Index BAM files
        SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

        // Extract mapped reads from SAM files as FASTQ
        EXTRACT_MAPPED_READS ( MINIMAP2_ALIGN.out.sam )

        // Join BAM and BAI files
        ch_bam_bai = SAMTOOLS_SORT.out.bam
            .join(SAMTOOLS_INDEX.out.bai)

        // Collect all version files from modules
        versions_ch = Channel.empty()
        versions_ch = versions_ch
            .mix(MINIMAP2_ALIGN.out.versions)
            .mix(SAMTOOLS_VIEW.out.versions)
            .mix(SAMTOOLS_SORT.out.versions)
            .mix(SAMTOOLS_INDEX.out.versions)
            .mix(EXTRACT_MAPPED_READS.out.versions)
            .collect()

    emit:
        mapped_bams   = ch_bam_bai                    // [ val(meta), path(bam), path(bai) ]
        mapped_fastqs = EXTRACT_MAPPED_READS.out.fastq // [ val(meta), path(fastq) ]
        versions      = versions_ch
}
