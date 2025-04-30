workflow PROCESS_SAMPLE_SHEET {

    take:
        sample_sheet //path sample sheet
    main:

    samples = sample_sheet
        .splitCsv(header:true,sep:",")
        .map { make_fastq_ch(it) }

    emit:
        samples
}

def make_fastq_ch(LinkedHashMap row){
    def meta = [:]
    meta.sample = row.sample.trim()
    meta.replicate_id = row.replicate_id.trim()
    
    meta_fastq = [ meta, [file(row.fastq1.trim()), file(row.fastq2.trim())] ]

    return meta_fastq

}