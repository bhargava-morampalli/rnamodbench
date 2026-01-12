#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: QC_STATS
 * Purpose: Generate mapping statistics, depth, coverage plots, and NanoPlot QC for BAM files.
 *
 * Inputs:
 *   - bams: [ val(meta), path(bam), path(bai) ] - unified channel for all BAM files
 *           meta should contain: id, type, rrna, replicate
 *
 * Outputs:
 *   - flagstat: [ val(meta), path(flagstat) ]
 *   - depth:    [ val(meta), path(depth) ]
 *   - plots:    [ val(meta), path(plot) ]
 *   - nanoplot: [ val(meta), path(stats) ]
 *   - versions: [ path(versions.yml) ]
 *
 * Note: This subworkflow processes all BAM files regardless of target type.
 *       Results are organized by meta.rrna and meta.type in output directories.
 */

include { SAMTOOLS_FLAGSTAT } from '../../../modules/local/samtools_flagstat'
include { SAMTOOLS_DEPTH    } from '../../../modules/local/samtools_depth'
include { NANOPLOT_BAM      } from '../../../modules/local/nanoplot_bam'
include { COVERAGE_PLOT     } from '../../../modules/local/coverage_plot'

workflow QC_STATS {
    take:
        bams  // [ val(meta), path(bam), path(bai) ] - unified channel

    main:
        // Extract just BAM files (without BAI) for tools that don't need index
        ch_bams = bams.map { meta, bam, bai -> [ meta, bam ] }

        // Mapping statistics with samtools flagstat
        SAMTOOLS_FLAGSTAT ( ch_bams )

        // Calculate depth per position
        SAMTOOLS_DEPTH ( ch_bams )

        // Generate coverage plots
        COVERAGE_PLOT ( SAMTOOLS_DEPTH.out.depth )

        // Generate read statistics with NanoPlot
        NANOPLOT_BAM ( ch_bams )

        // Collect all version files from modules
        versions_ch = Channel.empty()
        versions_ch = versions_ch
            .mix(SAMTOOLS_FLAGSTAT.out.versions)
            .mix(SAMTOOLS_DEPTH.out.versions)
            .mix(COVERAGE_PLOT.out.versions)
            .mix(NANOPLOT_BAM.out.versions)
            .collect()

    emit:
        flagstat = SAMTOOLS_FLAGSTAT.out.flagstat  // [ val(meta), path(flagstat) ]
        depth    = SAMTOOLS_DEPTH.out.depth        // [ val(meta), path(depth) ]
        plots    = COVERAGE_PLOT.out.plot          // [ val(meta), path(plot) ]
        nanoplot = NANOPLOT_BAM.out.stats          // [ val(meta), path(stats) ]
        versions = versions_ch                     // [ path(versions.yml) ]
}
