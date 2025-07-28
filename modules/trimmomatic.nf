process TRIMMOMATIC {
    input:
        tuple val(meta), path(paired_reads)
        val trimmomatic_jarfile

    output:
        tuple val(meta), path("${paired_reads[0].baseName}.r[1|2].fastq.gz")

    script:
    def basename = "${paired_reads[0].baseName}"
    """
        java -jar ${trimmomatic_jarfile} PE \
        ${paired_reads} \
        ${basename}.r1.fastq.gz ${basename}.r1_unpaired.fastq.gz \
        ${basename}.r2.fastq.gz ${basename}.r2_unpaired.fastq.gz \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 MINLEN:36
    """
}