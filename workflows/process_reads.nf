include { FASTQC } from "../modules/fastqc"
include { FASTP } from "../modules/fastp.nf"

workflow READS_QC {
    take:
        samples // metadata, tuple ([samples,replicate_id:r],[fastq1,fastq2])

    main:
        // each fastqc channel tuple contains a meta, replicate, [path(read1), path(read2)]
        // fastp reports the same read-level statistics in its own JSON/HTML, so FastQC
        // is independently skippable without also giving up trimming -- which is what
        // skip_qc costs, since it bypasses this whole subworkflow.
        if (!params.skip_fastqc) {
            FASTQC(samples)
        }

        fastp_out = FASTP(samples)

        // fastp always writes split R1/R2, including when the input was interleaved,
        // so everything downstream must see a genuine pair -- BWA_MEM branches on
        // meta.interleaved. Build a fresh map here rather than mutating the shared
        // one, which would alias across every item in the channel.
        processed_samples = fastp_out.reads
            .map { meta, reads -> [meta + [interleaved: false], reads] }

    emit:
        processed_samples
}
