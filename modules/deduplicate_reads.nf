process DEDUPLICATE_READS {
    publishDir "${params.output_dir}/consensus_bams/", mode: 'copy', pattern: "*.bam"
    errorStrategy 'retry'
    maxRetries 3

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
        samtools sort -n ${bam} | samtools fixmate -mr - ${bam.baseName}.fixed.bam
        samtools sort ${bam.baseName}.fixed.bam | samtools markdup --duplicate-count --barcode-name - ${bam.baseName}.dedup.bam
        samtools index ${bam.baseName}.dedup.bam
        """

    // requires testing
    else if(meta.sequencing_tech.toLowerCase() == "hybrid-capture")
        error "Hybrid Capture deduplication not yet implemented."
    //    """
    //    samtools markdup -r ${bam} ${bam.baseName}.dedup.bam
    //    """

    else
        error "Sequencing tech ${meta.sequencing_tech} must be MIPS or hybrid-capture for read deduplication."
}