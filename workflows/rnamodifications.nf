/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnamodifications.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK          } from '../subworkflows/local/input_check'
include { MAPPING_RRNA         } from '../subworkflows/local/mapping_rrna'
include { QC_STATS             } from '../subworkflows/local/qc_stats'
include { PREPARE_SIGNAL_DATA  } from '../subworkflows/local/prepare_signal_data'
include { SIGNAL_PROCESSING    } from '../subworkflows/local/signal_processing'
include { MODIFICATION_CALLING } from '../subworkflows/local/modification_calling'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Function to generate software versions from channel
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
            .unique()
            .map { it.text }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAMODIFICATIONS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        params.input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // Create reference channels
    //
    ch_ref_16s = Channel.value(file(params.ref_16s, checkIfExists: true))
    ch_ref_23s = Channel.value(file(params.ref_23s, checkIfExists: true))

    //
    // SUBWORKFLOW: Map reads to 16S and 23S rRNA references
    //
    MAPPING_RRNA (
        INPUT_CHECK.out.native_reads,
        INPUT_CHECK.out.ivt_reads,
        ch_ref_16s,
        ch_ref_23s
    )
    ch_versions = ch_versions.mix(MAPPING_RRNA.out.versions)

    //
    // SUBWORKFLOW: Generate QC statistics and coverage plots
    //
    QC_STATS (
        MAPPING_RRNA.out.bams_16s_native,
        MAPPING_RRNA.out.bams_16s_ivt,
        MAPPING_RRNA.out.bams_23s_native,
        MAPPING_RRNA.out.bams_23s_ivt
    )
    ch_versions = ch_versions.mix(QC_STATS.out.versions)

    //
    // SUBWORKFLOW: Prepare signal-level data (FAST5 extraction, tombo resquiggle, f5c)
    //
    PREPARE_SIGNAL_DATA (
        MAPPING_RRNA.out.mapped_fastq_16s_native,
        MAPPING_RRNA.out.mapped_fastq_16s_ivt,
        MAPPING_RRNA.out.mapped_fastq_23s_native,
        MAPPING_RRNA.out.mapped_fastq_23s_ivt,
        MAPPING_RRNA.out.mapped_bam_16s_native,
        MAPPING_RRNA.out.mapped_bam_16s_ivt,
        MAPPING_RRNA.out.mapped_bam_23s_native,
        MAPPING_RRNA.out.mapped_bam_23s_ivt,
        INPUT_CHECK.out.native_fast5,
        INPUT_CHECK.out.ivt_fast5,
        ch_ref_16s,
        ch_ref_23s
    )
    ch_versions = ch_versions.mix(PREPARE_SIGNAL_DATA.out.versions)

    //
    // SUBWORKFLOW: Signal processing and eventalign
    //
    SIGNAL_PROCESSING (
        PREPARE_SIGNAL_DATA.out.f5c_ready_16s_native,
        PREPARE_SIGNAL_DATA.out.f5c_ready_16s_ivt,
        PREPARE_SIGNAL_DATA.out.f5c_ready_23s_native,
        PREPARE_SIGNAL_DATA.out.f5c_ready_23s_ivt,
        ch_ref_16s,
        ch_ref_23s
    )
    ch_versions = ch_versions.mix(SIGNAL_PROCESSING.out.versions)

    //
    // SUBWORKFLOW: RNA modification calling with multiple tools
    //
    MODIFICATION_CALLING (
        PREPARE_SIGNAL_DATA.out.tombo_resquiggled_16s,
        PREPARE_SIGNAL_DATA.out.tombo_resquiggled_23s,
        SIGNAL_PROCESSING.out.eventalign_16s_native,
        SIGNAL_PROCESSING.out.eventalign_16s_ivt,
        SIGNAL_PROCESSING.out.eventalign_23s_native,
        SIGNAL_PROCESSING.out.eventalign_23s_ivt
    )
    ch_versions = ch_versions.mix(MODIFICATION_CALLING.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_versions.yml',
            sort: true,
            newLine: true
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
