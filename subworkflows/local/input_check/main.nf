//
// Check input samplesheet and get read channels
//

// No modules needed - reading CSV directly

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    ch_versions = Channel.empty()

    reads = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.single_end = true
            meta.type = row.type
            meta.replicate = row.replicate
            [ meta, file(row.fastq), row.type, row.replicate, file(row.fast5_dir) ]
        }

    branched_reads = reads
        .branch {
            native: it[2] == 'native'
            ivt: it[2] == 'ivt'
        }

    ch_native_fast5 = reads
        .filter { it[2] == 'native' }
        .map { it[4] }
        .first()

    ch_ivt_fast5 = reads
        .filter { it[2] == 'ivt' }
        .map { it[4] }
        .first()

    emit:
    native_reads = branched_reads.native.map { meta, fastq, t, r, f -> [ meta, fastq ] }
    ivt_reads    = branched_reads.ivt.map { meta, fastq, t, r, f -> [ meta, fastq ] }
    native_fast5 = ch_native_fast5
    ivt_fast5    = ch_ivt_fast5
    versions     = ch_versions
}
