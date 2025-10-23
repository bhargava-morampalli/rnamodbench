#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Subworkflow: QC_STATS
 * Purpose: Generate mapping statistics, depth, coverage plots, and NanoPlot QC for BAM files.
 * Inputs:
 *   - bams_16s_native: [ val(meta), path(bam) ]
 *   - bams_16s_ivt:    [ val(meta), path(bam) ]
 *   - bams_23s_native: [ val(meta), path(bam) ]
 *   - bams_23s_ivt:    [ val(meta), path(bam) ]
 * Outputs:
 *   - flagstat: [ val(meta), path(flagstat) ]
 *   - depth:    [ val(meta), path(depth) ]
 *   - plots:    [ val(meta), path(plot) ]
 *   - nanoplot: [ val(meta), path(stats) ]
 *   - versions: [ path(versions.yml) ]
 */

include { SAMTOOLS_FLAGSTAT } from '../../modules/local/samtools_flagstat'
include { SAMTOOLS_DEPTH    } from '../../modules/local/samtools_depth'
include { NANOPLOT_BAM      } from '../../modules/local/nanoplot_bam'
include { COVERAGE_PLOT     } from '../../modules/local/coverage_plot'

workflow QC_STATS {
    take:
        bams_16s_native  // [ val(meta), path(bam) ]
        bams_16s_ivt     // [ val(meta), path(bam) ]
        bams_23s_native  // [ val(meta), path(bam) ]
        bams_23s_ivt     // [ val(meta), path(bam) ]

    main:
        // Combine all BAMs for QC
        ch_all_bams = bams_16s_native
            .mix(bams_16s_ivt)
            .mix(bams_23s_native)
            .mix(bams_23s_ivt)

        // Mapping statistics with samtools flagstat
        SAMTOOLS_FLAGSTAT ( ch_all_bams )

        // Calculate depth per position
        SAMTOOLS_DEPTH ( ch_all_bams )

        // Generate coverage plots
        COVERAGE_PLOT ( SAMTOOLS_DEPTH.out.depth )

        // Generate read statistics with NanoPlot
        NANOPLOT_BAM ( ch_all_bams )

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
