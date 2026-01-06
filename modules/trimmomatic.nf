process TRIMMOMATIC {
    tag "${meta.replicate_id}"

    cpus 1
    memory '2G'
    time '4h'

    input:
        tuple val(meta), path(paired_reads)
        val trimmomatic_jarfile

    output:
        tuple val(meta), path("${meta.replicate_id}_paired.r[1|2].fastq.gz")

    script:
    def rep_id = meta.replicate_id
    """
        java -jar ${trimmomatic_jarfile} PE \
        ${paired_reads} \
        ${rep_id}_paired.r1.fastq.gz ${rep_id}.r1_unpaired.fastq.gz \
        ${rep_id}_paired.r2.fastq.gz ${rep_id}.r2_unpaired.fastq.gz \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 MINLEN:36
    """
}