//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    // Separate native and IVT reads
    reads
        .branch {
            native: it[2] == 'native'
            ivt: it[2] == 'ivt'
        }
        .set { branched_reads }

    // Extract FAST5 paths
    ch_native_fast5 = reads
        .filter { it[2] == 'native' }
        .map { it[4] }
        .first()

    ch_ivt_fast5 = reads
        .filter { it[2] == 'ivt' }
        .map { it[4] }
        .first()

    emit:
    native_reads = branched_reads.native.map { meta, fastq, type, replicate, fast5 -> [ meta, fastq ] } // channel: [ val(meta), path(fastq) ]
    ivt_reads    = branched_reads.ivt.map { meta, fastq, type, replicate, fast5 -> [ meta, fastq ] }    // channel: [ val(meta), path(fastq) ]
    native_fast5 = ch_native_fast5                                                                       // channel: path(fast5_dir)
    ivt_fast5    = ch_ivt_fast5                                                                          // channel: path(fast5_dir)
    versions     = ch_versions                                                                           // channel: [ versions.yml ]
}

// Function to get list of [ meta, fastq, type, replicate, fast5_dir ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = true // Direct RNA-seq is always single-end
    meta.type         = row.type // native or ivt
    meta.replicate    = row.replicate

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq}"
    }
    fastq_meta = [ meta, file(row.fastq), row.type, row.replicate, file(row.fast5_dir) ]
    return fastq_meta
}
