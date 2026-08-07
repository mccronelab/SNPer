process DEDUPLICATE_READS {
    errorStrategy 'retry'
    maxRetries 3
    tag "${meta.replicate_id}"

    cpus { 1 * task.attempt }
    // need to account for potentially increasing CPU allocation
    memory { 2G * task.attempt * task.cpus }
    time { 4.h * task.attempt }

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("${bam.baseName}.dedup.bam"), path("*.bai")

    script:
    // samtools -r markdup will remove duplicate reads entirely
    if(meta.read_deduplication.toLowerCase() == "umi")
        """
        samtools collate -@ ${task.cpus} -O -u ${bam} -T ./${meta.replicate_id}.TEMP \
        | samtools fixmate -@ ${task.cpus}  -m -u - - \
        | samtools sort -@ ${task.cpus}  -u - \
        | samtools markdup -@ ${task.cpus} --barcode-name -r -f ${meta.replicate_id}.stats.txt - ${bam.baseName}.dedup.bam
        samtools index ${bam.baseName}.dedup.bam
        """

    else if(meta.read_deduplication.toLowerCase() == "positional")
        """
        samtools collate -@ ${task.cpus} -O ${bam} \
        | samtools fixmate -m -@ ${task.cpus} - - \
        | samtools sort -@ ${task.cpus} - \
        | samtools markdup -r -@ ${task.cpus} -s -f ${meta.replicate_id}.stats.txt - ${bam.baseName}.dedup.bam

        samtools index ${bam.baseName}.dedup.bam
        """

    else
        error "read_deduplication value ${meta.read_deduplication} must be umi or positional for read deduplication."
}