//
// RNA modification calling using multiple tools
//

include { TOMBO_DETECT_MODIFICATIONS } from '../../modules/local/tombo_detect_modifications'
include { TOMBO_TEXT_OUTPUT          } from '../../modules/local/tombo_text_output'
include { YANOCOMP_PREPARE           } from '../../modules/local/yanocomp_prepare'
include { YANOCOMP_ANALYSIS          } from '../../modules/local/yanocomp_analysis'
include { XPORE_DATAPREP             } from '../../modules/local/xpore_dataprep'
include { XPORE_DIFFMOD              } from '../../modules/local/xpore_diffmod'

workflow MODIFICATION_CALLING {
    take:
    tombo_resquiggled_16s    // channel: [ key, [types], [fast5s] ]
    tombo_resquiggled_23s    // channel: [ key, [types], [fast5s] ]
    eventalign_16s_native    // channel: [ val(meta), path(eventalign) ]
    eventalign_16s_ivt       // channel: [ val(meta), path(eventalign) ]
    eventalign_23s_native    // channel: [ val(meta), path(eventalign) ]
    eventalign_23s_ivt       // channel: [ val(meta), path(eventalign) ]

    main:
    ch_versions = Channel.empty()

    //
    // TOMBO: Detect modifications by comparing native vs IVT
    //
    TOMBO_DETECT_MODIFICATIONS (
        tombo_resquiggled_16s.mix(tombo_resquiggled_23s)
    )
    ch_versions = ch_versions.mix(TOMBO_DETECT_MODIFICATIONS.out.versions.first())

    //
    // TOMBO: Extract text output from wiggle files
    //
    TOMBO_TEXT_OUTPUT ( TOMBO_DETECT_MODIFICATIONS.out.stats )
    ch_versions = ch_versions.mix(TOMBO_TEXT_OUTPUT.out.versions.first())

    //
    // YANOCOMP: Prepare HDF5 files from f5c eventalign
    //
    ch_eventalign_all = eventalign_16s_native
        .mix(eventalign_16s_ivt)
        .mix(eventalign_23s_native)
        .mix(eventalign_23s_ivt)

    YANOCOMP_PREPARE ( ch_eventalign_all )
    ch_versions = ch_versions.mix(YANOCOMP_PREPARE.out.versions.first())

    //
    // Join native and IVT HDF5 files for comparison
    //
    YANOCOMP_PREPARE.out.hdf5
        .map { meta, hdf5 ->
            def key = "${meta.rrna}_${meta.replicate}"
            [ key, meta.type, hdf5 ]
        }
        .groupTuple()
        .map { key, types, hdf5s ->
            def meta = [:]
            meta.id = key
            meta.rrna = key.split('_')[0]
            meta.replicate = key.split('_')[1]

            def native_idx = types.findIndexOf { it == 'native' }
            def ivt_idx = types.findIndexOf { it == 'ivt' }

            [ meta, hdf5s[native_idx], hdf5s[ivt_idx] ]
        }
        .set { ch_yanocomp_ready }

    //
    // YANOCOMP: Detect modifications
    //
    YANOCOMP_ANALYSIS ( ch_yanocomp_ready )
    ch_versions = ch_versions.mix(YANOCOMP_ANALYSIS.out.versions.first())

    //
    // XPORE: Prepare data
    //
    XPORE_DATAPREP ( ch_eventalign_all )
    ch_versions = ch_versions.mix(XPORE_DATAPREP.out.versions.first())

    //
    // Join native and IVT xpore data for comparison
    //
    XPORE_DATAPREP.out.dataprep
        .map { meta, dataprep_dir ->
            def key = "${meta.rrna}_${meta.replicate}"
            [ key, meta.type, dataprep_dir ]
        }
        .groupTuple()
        .set { ch_xpore_ready }

    //
    // XPORE: Differential modification analysis
    //
    XPORE_DIFFMOD ( ch_xpore_ready )
    ch_versions = ch_versions.mix(XPORE_DIFFMOD.out.versions.first())

    emit:
    tombo_bed    = TOMBO_TEXT_OUTPUT.out.bed         // channel: [ val(meta), path(bed) ]
    yanocomp_bed = YANOCOMP_ANALYSIS.out.bed         // channel: [ val(meta), path(bed) ]
    xpore_table  = XPORE_DIFFMOD.out.diffmod         // channel: [ val(meta), path(table) ]
    versions     = ch_versions                       // channel: [ versions.yml ]
}
