//
// Generate QC statistics and coverage plots
//

include { SAMTOOLS_FLAGSTAT } from '../../modules/local/samtools_flagstat'
include { SAMTOOLS_DEPTH    } from '../../modules/local/samtools_depth'
include { NANOPLOT_BAM      } from '../../modules/local/nanoplot_bam'
include { COVERAGE_PLOT     } from '../../modules/local/coverage_plot'

workflow QC_STATS {
    take:
    bams_16s_native  // channel: [ val(meta), path(bam) ]
    bams_16s_ivt     // channel: [ val(meta), path(bam) ]
    bams_23s_native  // channel: [ val(meta), path(bam) ]
    bams_23s_ivt     // channel: [ val(meta), path(bam) ]

    main:
    ch_versions = Channel.empty()

    //
    // Combine all BAMs for QC
    //
    ch_all_bams = bams_16s_native
        .mix(bams_16s_ivt)
        .mix(bams_23s_native)
        .mix(bams_23s_ivt)

    //
    // MODULE: Mapping statistics with samtools flagstat
    //
    SAMTOOLS_FLAGSTAT ( ch_all_bams )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    //
    // MODULE: Calculate depth per position
    //
    SAMTOOLS_DEPTH ( ch_all_bams )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    //
    // MODULE: Generate coverage plots
    //
    COVERAGE_PLOT ( SAMTOOLS_DEPTH.out.depth )

    //
    // MODULE: Generate read statistics with NanoPlot
    //
    NANOPLOT_BAM ( ch_all_bams )
    ch_versions = ch_versions.mix(NANOPLOT_BAM.out.versions.first())

    emit:
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat  // channel: [ val(meta), path(flagstat) ]
    depth    = SAMTOOLS_DEPTH.out.depth        // channel: [ val(meta), path(depth) ]
    plots    = COVERAGE_PLOT.out.plot          // channel: [ val(meta), path(plot) ]
    nanoplot = NANOPLOT_BAM.out.stats          // channel: [ val(meta), path(stats) ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
