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
    if(meta.sequencing_tech.toLowerCase() == "mips") 
        """
        samtools collate -@ ${task.cpus} -O -u ${bam} -T ./${meta.replicate_id}.TEMP \
        | samtools fixmate -@ ${task.cpus}  -m -u - - \
        | samtools sort -@ ${task.cpus}  -u - \
        | samtools markdup -@ ${task.cpus} --barcode-name -r -f ${meta.replicate_id}.stats.txt - ${bam.baseName}.dedup.bam
        samtools index ${bam.baseName}.dedup.bam
        """

    // requires testing
    else if(meta.sequencing_tech.toLowerCase() == "hybrid-capture")
        """
        samtools collate -@ ${task.cpus} -O ${bam} \
        | samtools fixmate -m -@ ${task.cpus} - - \
        | samtools sort -@ ${task.cpus} - \
        | samtools markdup -r -@ ${task.cpus} -s -f ${meta.replicate_id}.stats.txt - ${bam.baseName}.dedup.bam

        samtools index ${bam.baseName}.dedup.bam
        """

    else
        error "Sequencing tech ${meta.sequencing_tech} must be MIPS or hybrid-capture for read deduplication."
}