workflow PROCESS_SAMPLE_SHEET {

    take:
        sample_sheet //path sample sheet
        primer_csv //path to primer CSV
        interleaved //value that should be boolean
        read_deduplication_param //string value that should be one of 'amplicon', 'umi', 'positional'
        primer_id_default // string that should match a primer ID in primer_csv

    main:

    // returns tuple with [metadata, [reads]]
    // depending on input, may be of form [metadata, 'reads1, reads2'] or
    // [metadata, [reads_interleaved]]
    samples = parse_sample_sheet(sample_sheet, interleaved, read_deduplication_param, primer_id_default)
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

// Resolve a path read out of a CSV. A relative path resolves against the CSV's
// own parent directory rather than launchDir, so a sheet is portable: it can be
// moved, or run from any directory, and `nextflow run <github-url>` works — there
// launchDir is wherever the user happened to be standing, while projectDir is the
// pulled asset clone, and neither is where the sheet's data lives. Absolute paths
// and remote URIs (s3://, https://, …) are already unambiguous and pass through
// untouched.
def resolve_against_csv(csv_path, raw_path) {
    def path_str = raw_path.toString().trim()

    if (path_str.startsWith('/') || path_str ==~ /^[a-zA-Z][a-zA-Z0-9+.-]*:\/\/.*$/) {
        return file(path_str)
    }

    return file("${file(csv_path).parent}/${path_str}")
}

def parse_sample_sheet(sample_sheet, interleaved_param, read_deduplication_param,
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
        def fastq1 = resolve_against_csv(params.sample_sheet, row.fastq1)

        meta.sample = sample_id
        meta.interleaved = interleaved
        meta.primer_id = primer_id

        // handle optional columns, depending on whether input is interleaved
        def replicate_id = row.containsKey('replicate_id') && row.replicate_id ?
            row.replicate_id : sample_id

        // 'sequencing_tech' was renamed to 'read_deduplication'. Fail loudly rather
        // than let an old sheet fall through to the param default.
        if (row.containsKey('sequencing_tech')) {
            error "Samplesheet column 'sequencing_tech' has been renamed to 'read_deduplication'; rename it in: ${params.sample_sheet}"
        }

        def read_deduplication = row.containsKey('read_deduplication') && row.read_deduplication ?
            row.read_deduplication : read_deduplication_param

        meta.replicate_id = replicate_id
        meta.read_deduplication = read_deduplication

        // verify we handle the read deduplication method
        def renamed_methods = ["mips": "umi", "hybrid-capture": "positional"]
        if (read_deduplication.toLowerCase() in renamed_methods.keySet()) {
            error "read_deduplication value '${read_deduplication}' has been renamed to '${renamed_methods[read_deduplication.toLowerCase()]}'"
        }

        def supported_methods = ["amplicon", "umi", "positional"]
        if (!(read_deduplication.toLowerCase() in supported_methods)) {
            throw new IllegalArgumentException("Value '${read_deduplication}' not in ${supported_methods}")
}

        // return depending on number of input files (1 if interleaved, 2 otherwise)
        if(interleaved) {
            return tuple(meta, [fastq1])
        }
        else {
            if (!row.containsKey('fastq2') || !row.fastq2.trim()) {
                error "Samplesheet missing required 'fastq2' column or empty value in row: ${row}"
            }

            def fastq2 = resolve_against_csv(params.sample_sheet, row.fastq2)

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

        // groovy supports implicit returns.
        // Anchored on the primer CSV's own directory for the same reason as the
        // sample sheet. This is a no-op for the default `${projectDir}/primers.csv`
        // — its parent *is* projectDir — but it stops a user-supplied primer CSV
        // elsewhere on disk from having projectDir glued onto its paths.
        [row.primer_id, resolve_against_csv(params.primer_csv, row.primer_bedfile)]
    }
}