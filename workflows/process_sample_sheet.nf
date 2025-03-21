workflow PROCESS_SAMPLE_SHEET {

    take:
        sample_sheet //path sample sheet
    main:

     samples = sample_sheet
        .splitCsv(header:true,sep:",")
        .map{row -> tuple(row.sample.trim(), row.replicate_id.trim(), [file(row.fastq1.trim()), file(row.fastq2.trim())])}

    emit:
        samples
}