process BWA_MEM {
    label 'process_high'
    errorStrategy 'retry'
    maxRetries 3

    cpus {1 * task.attempt}
    // need to account for potentially increasing CPU allocation
    memory { 2G * task.attempt * task.cpus}
    time { 4.h * task.attempt }


    input:
        tuple val(meta),  path(paired_reads), path(reference)

    output:
        tuple val(meta), path("*.bam")

    script:
    if( meta.interleaved == true)
        """
        bwa index ${reference}
        bwa mem -p ${reference} ${paired_reads} | \
        samtools view -b - > ${meta.replicate_id}.bam
        """
    else
        """
        bwa index ${reference}
        bwa mem ${reference} ${paired_reads} | \
        samtools view -b - > ${meta.replicate_id}.bam
        """
}

process BWA_REMAP {
    label 'process_high'
    errorStrategy 'retry'
    publishDir "${params.output_dir}/realigned_bams/", mode: "copy"
    maxRetries 3
    //https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
    cpus {1 * task.attempt}
    // need to account for potentially increasing CPU allocation
    memory { 2G * task.attempt * task.cpus}
    time { 4.h * task.attempt }


    input:
        tuple val(meta),  path(bam), path(reference)

    output:
        tuple val(meta), path("*.remap.bam")

    script:
        """
        bwa index ${reference}
        samtools collate -Oun128 ${bam} -T ./collate.TEMP \
        | samtools fastq - \
        | bwa mem -p ${reference} - \
        | samtools view -b - > ${meta.replicate_id}.remap.bam
        """
}