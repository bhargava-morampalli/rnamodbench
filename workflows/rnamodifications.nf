//
// RNA Modifications workflow - main analysis pipeline
//

include { MAPPING_RRNA } from '../subworkflows/local/mapping_rrna'
include { QC_STATS } from '../subworkflows/local/qc_stats'
include { PREPARE_SIGNAL_DATA } from '../subworkflows/local/prepare_signal_data'
include { SIGNAL_PROCESSING } from '../subworkflows/local/signal_processing'
include { MODIFICATION_CALLING } from '../subworkflows/local/modification_calling'

workflow RNAMODIFICATIONS {

    ch_versions = Channel.empty()

    // Read and parse samplesheet
    ch_samplesheet = Channel.fromPath(params.input)

    ch_reads_all = ch_samplesheet
        .splitCsv(header: true, sep: ',')
        .map { row ->
            [
                [ id: row.sample, single_end: true, type: row.type, replicate: row.replicate, fast5_dir: file(row.fast5_dir) ],
                file(row.fastq)
            ]
        }

    // Branch reads by type
    ch_branched_reads = ch_reads_all.branch { meta, fastq ->
        native: meta.type == 'native'
        ivt: meta.type == 'ivt'
    }

    // Extract read channels
    ch_native_reads = ch_branched_reads.native
    ch_ivt_reads = ch_branched_reads.ivt

    // The fast5_dir is now carried within the meta map of ch_native_reads and ch_ivt_reads
    // and subsequently within the meta map of MAPPING_RRNA.out.mapped_fastq_... channels.
    // PREPARE_SIGNAL_DATA will extract it from there.

    // Create reference channels
    ch_ref_16s = Channel.value(file(params.ref_16s, checkIfExists: true))
    ch_ref_23s = Channel.value(file(params.ref_23s, checkIfExists: true))

    // SUBWORKFLOW: Map reads to rRNA references
    MAPPING_RRNA (
        ch_native_reads,
        ch_ivt_reads,
        ch_ref_16s,
        ch_ref_23s
    )
    ch_versions = ch_versions.mix(MAPPING_RRNA.out.versions)

    // SUBWORKFLOW: QC statistics
    QC_STATS (
        MAPPING_RRNA.out.bams_16s_native,
        MAPPING_RRNA.out.bams_16s_ivt,
        MAPPING_RRNA.out.bams_23s_native,
        MAPPING_RRNA.out.bams_23s_ivt
    )
    ch_versions = ch_versions.mix(QC_STATS.out.versions)

    // SUBWORKFLOW: Prepare signal data
    PREPARE_SIGNAL_DATA (
        MAPPING_RRNA.out.mapped_fastq_16s_native,
        MAPPING_RRNA.out.mapped_fastq_16s_ivt,
        MAPPING_RRNA.out.mapped_fastq_23s_native,
        MAPPING_RRNA.out.mapped_fastq_23s_ivt,
        MAPPING_RRNA.out.mapped_bam_16s_native,
        MAPPING_RRNA.out.mapped_bam_16s_ivt,
        MAPPING_RRNA.out.mapped_bam_23s_native,
        MAPPING_RRNA.out.mapped_bam_23s_ivt,
        ch_ref_16s,
        ch_ref_23s
    )
    ch_versions = ch_versions.mix(PREPARE_SIGNAL_DATA.out.versions)

    // SUBWORKFLOW: Signal processing
    SIGNAL_PROCESSING (
        PREPARE_SIGNAL_DATA.out.f5c_ready_16s_native,
        PREPARE_SIGNAL_DATA.out.f5c_ready_16s_ivt,
        PREPARE_SIGNAL_DATA.out.f5c_ready_23s_native,
        PREPARE_SIGNAL_DATA.out.f5c_ready_23s_ivt,
        ch_ref_16s,
        ch_ref_23s
    )
    ch_versions = ch_versions.mix(SIGNAL_PROCESSING.out.versions)

    // SUBWORKFLOW: Modification calling
    MODIFICATION_CALLING (
        PREPARE_SIGNAL_DATA.out.tombo_resquiggled_16s,
        PREPARE_SIGNAL_DATA.out.tombo_resquiggled_23s,
        SIGNAL_PROCESSING.out.eventalign_16s_native,
        SIGNAL_PROCESSING.out.eventalign_16s_ivt,
        SIGNAL_PROCESSING.out.eventalign_23s_native,
        SIGNAL_PROCESSING.out.eventalign_23s_ivt
    )
    ch_versions = ch_versions.mix(MODIFICATION_CALLING.out.versions)

    // Collect software versions
    ch_versions
        .unique()
        .map { it.text }
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'software_versions.yml',
            sort: true,
            newLine: true
        )
}
