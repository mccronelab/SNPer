process NGMERGE {
    input:
        tuple val(meta), path(paired_reads)

    output:
        tuple val(meta), path("${meta.replicate_id}.adapter_removed_[1|2].fastq.gz")

    script:
    def read1 = paired_reads[0]
    def read2 = paired_reads[1]
        """
        NGmerge -a -1 ${read1} -2 ${read2} -o ${meta.replicate_id}.adapter_removed
        """
}