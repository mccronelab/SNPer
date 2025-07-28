include { FASTQC } from "../modules/fastqc"
include { TRIMMOMATIC } from "../modules/trimmomatic.nf"
include { NGMERGE } from "../modules/ngmerge.nf"

workflow READS_QC {
    take:
        samples // metadata, tuple ([samples,replicate_id:r],[fastq1,fastq2])

    main:
        // each fastqc channel tuple contains a meta, replicate, [path(read1), path(read2)]
        FASTQC(samples)

        trimmed_samples = TRIMMOMATIC(samples, file(params.trimmomatic_jarfile))
        merged_samples = NGMERGE(trimmed_samples)

    emit:
        processed_samples = merged_samples
}