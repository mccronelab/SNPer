workflow PROCESS_SAMPLE_SHEET {

    take:
        sample_sheet //path sample sheet
        primer_csv //path to primer CSV
        interleaved //value that should be boolean
        sequencing_type_param //string value that should be one of 'amplicon', 'mips', 'hybrid_capture'
        primer_id_default // string that should match a primer ID in primer_csv

    main:

    // returns tuple with [metadata, [reads]]
    // depending on input, may be of form [metadata, 'reads1, reads2'] or
    // [metadata, [reads_interleaved]]
    samples = parse_sample_sheet(sample_sheet, interleaved, sequencing_type_param, primer_id_default).view()
    // returns Channel where each row looks like: [primer_id, primer_file]
    primer_paths = parse_primer_csv(primer_csv)
    
    // prepare samples for joining with primers
    samples_with_primer = samples
        .map { meta, fastqs -> tuple(meta.primer_id, meta, fastqs) }
        .combine(primer_paths, by: 0)
        // incorporate primer info into sample metadata map
        .map{ primer_id, meta, fastqs, primer_bed -> 
            meta.primer_id = primer_id
            meta.primer_bedfile = primer_bed
            
            return [meta, fastqs] }

    emit:
        samples_with_primer
}

def parse_sample_sheet(sample_sheet, interleaved_param, sequencing_technique_param,
    primer_id_default) {

    sample_sheet
    .splitCsv(header: true, sep: ',', strip: true)
    .map { row -> 
        // declare meta map and check for optional boolean columns
        def meta = [:]
        // check if column exists, and use param value if not
        def interleaved = row.interleaved?.toBoolean() ?: interleaved_param
        def primer_id = row.primer_id?.toString() ?: primer_id_default

        // validate required columns
        if (!row.containsKey('sample') || !row.sample?.trim()) {
            error "Samplesheet missing required 'sample' column or empty value in row: ${row}"
        }
        if (!row.containsKey('fastq1') || !row.fastq1.trim()) {
            error "Samplesheet missing required 'fastq1' column or empty value in row: ${row}"
        }

        def sample_id = row.sample
        def fastq1 = file(row.fastq1)

        meta.sample = sample_id
        meta.interleaved = interleaved
        meta.primer_id = primer_id

        // handle optional columns, depending on whether input is interleaved
        def replicate_id = row.containsKey('replicate_id') && row.replicate_id ?
            row.replicate_id : sample_id

        def sequencing_technique = row.containsKey('sequencing_tech') && row.sequencing_tech ?
            row.sequencing_tech : sequencing_technique_param

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

def parse_primer_csv(primer_csv) {
    // go through CSV file row by row, adding each primer ID as a key to the map
    // with a file path as the value
    primer_csv
    .splitCsv(header: true, sep: ',', strip: true)
    .map { row ->
        // check that primer paths, primer IDs are present
        if (!row.containsKey('primer_bedfile') || !row.primer_bedfile?.trim()) {
            error "Primer CSV missing required 'primer_bedfile' column or empty value in row: ${row}"
        }
        if (!row.containsKey('primer_id') || !row.primer_id?.trim()) {
            error "Primer CSV missing required 'primer_id' column or empty value in row: ${row}"
        }

        // groovy supports implicit returns
        [row.primer_id, file(row.primer_bedfile)]
    }
}