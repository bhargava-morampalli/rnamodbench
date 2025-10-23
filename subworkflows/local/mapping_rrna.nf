#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: MAPPING_RRNA
 * Purpose: Map reads to 16S and 23S rRNA references, extract mapped reads, and produce sorted/indexed BAMs for downstream analysis.
 * Inputs:
 *   - native_reads: [ val(meta), path(fastq) ]
 *   - ivt_reads:    [ val(meta), path(fastq) ]
 *   - ref_16s:      path(reference)
 *   - ref_23s:      path(reference)
 * Outputs:
 *   - bams_16s_native: [ val(meta), path(bam) ]
 *   - bams_16s_ivt:    [ val(meta), path(bam) ]
 *   - bams_23s_native: [ val(meta), path(bam) ]
 *   - bams_23s_ivt:    [ val(meta), path(bam) ]
 *   - mapped_fastq_16s_native: [ val(meta), path(fastq) ]
 *   - mapped_fastq_16s_ivt:    [ val(meta), path(fastq) ]
 *   - mapped_fastq_23s_native: [ val(meta), path(fastq) ]
 *   - mapped_fastq_23s_ivt:    [ val(meta), path(fastq) ]
 *   - mapped_bam_16s_native:   [ val(meta), path(bam) ]
 *   - mapped_bam_16s_ivt:      [ val(meta), path(bam) ]
 *   - mapped_bam_23s_native:   [ val(meta), path(bam) ]
 *   - mapped_bam_23s_ivt:      [ val(meta), path(bam) ]
 *   - versions:                [ path(versions.yml) ]
 */

include { MINIMAP2_ALIGN as MINIMAP2_16S_NATIVE  } from '../../modules/local/minimap2_align'
include { MINIMAP2_ALIGN as MINIMAP2_16S_IVT     } from '../../modules/local/minimap2_align'
include { MINIMAP2_ALIGN as MINIMAP2_23S_NATIVE  } from '../../modules/local/minimap2_align'
include { MINIMAP2_ALIGN as MINIMAP2_23S_IVT     } from '../../modules/local/minimap2_align'
include { SAMTOOLS_VIEW                           } from '../../modules/local/samtools_view'
include { SAMTOOLS_SORT                           } from '../../modules/local/samtools_sort'
include { SAMTOOLS_INDEX                          } from '../../modules/local/samtools_index'
include { EXTRACT_MAPPED_READS                    } from '../../modules/local/extract_mapped_reads'
include { MINIMAP2_ALIGN as REMAP_16S_NATIVE     } from '../../modules/local/minimap2_align'
include { MINIMAP2_ALIGN as REMAP_16S_IVT        } from '../../modules/local/minimap2_align'
include { MINIMAP2_ALIGN as REMAP_23S_NATIVE     } from '../../modules/local/minimap2_align'
include { MINIMAP2_ALIGN as REMAP_23S_IVT        } from '../../modules/local/minimap2_align'

workflow MAPPING_RRNA {
    take:
        native_reads    // [ val(meta), path(fastq) ]
        ivt_reads       // [ val(meta), path(fastq) ]
        ref_16s         // path(reference)
        ref_23s         // path(reference)

    main:
        // Add rRNA type to metadata
        native_reads_16s = native_reads.map { meta, fastq ->
            def new_meta = meta.clone()
            new_meta.rrna = '16s'
            [ new_meta, fastq ]
        }
        native_reads_23s = native_reads.map { meta, fastq ->
            def new_meta = meta.clone()
            new_meta.rrna = '23s'
            [ new_meta, fastq ]
        }
        ivt_reads_16s = ivt_reads.map { meta, fastq ->
            def new_meta = meta.clone()
            new_meta.rrna = '16s'
            [ new_meta, fastq ]
        }
        ivt_reads_23s = ivt_reads.map { meta, fastq ->
            def new_meta = meta.clone()
            new_meta.rrna = '23s'
            [ new_meta, fastq ]
        }

        // Map reads to references using minimap2
        MINIMAP2_16S_NATIVE ( native_reads_16s, ref_16s )
        MINIMAP2_16S_IVT    ( ivt_reads_16s, ref_16s )
        MINIMAP2_23S_NATIVE ( native_reads_23s, ref_23s )
        MINIMAP2_23S_IVT    ( ivt_reads_23s, ref_23s )

        // Combine all SAMs
        ch_all_sams = MINIMAP2_16S_NATIVE.out.sam
            .mix(MINIMAP2_16S_IVT.out.sam)
            .mix(MINIMAP2_23S_NATIVE.out.sam)
            .mix(MINIMAP2_23S_IVT.out.sam)

        // Convert SAM to BAM, sort and index
        SAMTOOLS_VIEW ( ch_all_sams )
        SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
        SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

        // Separate sorted BAMs by type and rRNA
        SAMTOOLS_SORT.out.bam
            .branch {
                native_16s: it[0].type == 'native' && it[0].rrna == '16s'
                native_23s: it[0].type == 'native' && it[0].rrna == '23s'
                ivt_16s: it[0].type == 'ivt' && it[0].rrna == '16s'
                ivt_23s: it[0].type == 'ivt' && it[0].rrna == '23s'
            }
            .set { sorted_bams }

        // Extract mapped reads from SAM files
        EXTRACT_MAPPED_READS ( ch_all_sams )

        // Separate mapped FASTQs by type and rRNA
        EXTRACT_MAPPED_READS.out.fastq
            .branch {
                native_16s: it[0].type == 'native' && it[0].rrna == '16s'
                native_23s: it[0].type == 'native' && it[0].rrna == '23s'
                ivt_16s: it[0].type == 'ivt' && it[0].rrna == '16s'
                ivt_23s: it[0].type == 'ivt' && it[0].rrna == '23s'
            }
            .set { mapped_fastqs }

        // Remap extracted reads for cleaner BAM files
        REMAP_16S_NATIVE ( mapped_fastqs.native_16s, ref_16s )
        REMAP_16S_IVT    ( mapped_fastqs.ivt_16s, ref_16s )
        REMAP_23S_NATIVE ( mapped_fastqs.native_23s, ref_23s )
        REMAP_23S_IVT    ( mapped_fastqs.ivt_23s, ref_23s )

        // Combine remapped SAMs
        ch_remapped_sams = REMAP_16S_NATIVE.out.sam
            .mix(REMAP_16S_IVT.out.sam)
            .mix(REMAP_23S_NATIVE.out.sam)
            .mix(REMAP_23S_IVT.out.sam)

        // Convert remapped SAMs to sorted BAMs
        SAMTOOLS_VIEW ( ch_remapped_sams )
        SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
        SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

        // Separate remapped BAMs
        SAMTOOLS_SORT.out.bam
            .branch {
                native_16s: it[0].type == 'native' && it[0].rrna == '16s'
                native_23s: it[0].type == 'native' && it[0].rrna == '23s'
                ivt_16s: it[0].type == 'ivt' && it[0].rrna == '16s'
                ivt_23s: it[0].type == 'ivt' && it[0].rrna == '23s'
            }
            .set { remapped_bams }

        // Collect all version files from modules
        versions_ch = Channel.empty()
        versions_ch = versions_ch
            .mix(MINIMAP2_16S_NATIVE.out.versions)
            .mix(MINIMAP2_16S_IVT.out.versions)
            .mix(MINIMAP2_23S_NATIVE.out.versions)
            .mix(MINIMAP2_23S_IVT.out.versions)
            .mix(SAMTOOLS_VIEW.out.versions)
            .mix(SAMTOOLS_SORT.out.versions)
            .mix(SAMTOOLS_INDEX.out.versions)
            .mix(EXTRACT_MAPPED_READS.out.versions)
            .mix(REMAP_16S_NATIVE.out.versions)
            .mix(REMAP_16S_IVT.out.versions)
            .mix(REMAP_23S_NATIVE.out.versions)
            .mix(REMAP_23S_IVT.out.versions)
            .collect()

    emit:
        bams_16s_native = sorted_bams.native_16s
        bams_16s_ivt    = sorted_bams.ivt_16s
        bams_23s_native = sorted_bams.native_23s
        bams_23s_ivt    = sorted_bams.ivt_23s

        mapped_fastq_16s_native = mapped_fastqs.native_16s
        mapped_fastq_16s_ivt    = mapped_fastqs.ivt_16s
        mapped_fastq_23s_native = mapped_fastqs.native_23s
        mapped_fastq_23s_ivt    = mapped_fastqs.ivt_23s

        mapped_bam_16s_native = remapped_bams.native_16s
        mapped_bam_16s_ivt    = remapped_bams.ivt_16s
        mapped_bam_23s_native = remapped_bams.native_23s
        mapped_bam_23s_ivt    = remapped_bams.ivt_23s

        versions = versions_ch
}
