process NGMERGE {
    input:
        tuple val(meta), path(paired_reads)

    output:
        tuple val(meta), path("*.merged.fastq.gz")

    script:
    def read1 = paired_reads[0]
    def read2 = paired_reads[1]
    def base = paired_reads[0].baseName
        """
        NGmerge -1 ${read1} -2 ${read2} -o ${base}.merged.fastq.gz
        """
}