workflow PROCESS_SAMPLE_SHEET {

    take:
        sample_sheet //path sample sheet
        interleaved //value that should be boolean
        sequencing_type_param //string value that should be one of 'amplicon', 'mips', 'hybrid_capture'

    main:

    // returns tuple with [metadata, [reads]]
    // depending on input, may be of form [metadata, 'reads1, reads2'] or
    // [metadata, [reads_interleaved]]
    samples = parse_sample_sheet(sample_sheet, interleaved, sequencing_type_param)

    emit:
        samples
}

def parse_sample_sheet(sample_sheet, interleaved_param, sequencing_technique_param) {

    sample_sheet
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        // declare meta hash and check for optional boolean columns
        def meta = [:]
        // check if column exists, and use param value if not
        def interleaved = row.interleaved?.toBoolean() ?: interleaved_param

        // validate required columns
        if (!row.containsKey('sample') || !row.sample?.trim()) {
            error "Samplesheet missing required 'sample' column or empty value in row: ${row}"
        }
        if (!row.containsKey('fastq1') || !row.fastq1.trim()) {
            error "Samplesheet missing required 'fastq1' column or empty value in row: ${row}"
        }

        def sample_id = row.sample.trim()
        def fastq1 = file(row.fastq1.trim())

        meta.sample = row.sample
        meta.interleaved = interleaved

        // handle optional columns, depending on whether input is interleaved
        def replicate_id = row.containsKey('replicate_id') && row.replicate_id?.trim() ?
            row.replicate_id.trim() : sample_id

        def sequencing_technique = row.containsKey('sequencing_tech') && row.sequencing_tech?.trim() ?
            row.sequencing_tech.trim() : sequencing_technique_param

        meta.replicate_id = replicate_id
        meta.sequencing_tech = sequencing_technique

        // return depending on number of input files (1 if interleaved, 2 otherwise)
        if(interleaved) {
            return tuple(meta, [fastq1])
        }
        else {
            if (!row.containsKey('fastq2') || !row.fastq2.trim()) {
                error "Samplesheet missing required 'fastq2' column or empty value in row: ${row}"
            }

            def fastq2 = file(row.fastq2.trim())

            return tuple(meta, [fastq1, fastq2])
        }
    }

}