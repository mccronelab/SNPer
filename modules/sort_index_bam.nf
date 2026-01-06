process SORT_INDEX_BAM {
    label 'process_low'
    cpus 1
    memory 2G
    time 12.h
    tag "${meta.replicate_id}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.sorted.bam"), path("*.bai")

    script:
    """
    samtools sort -o ${bam.baseName}.sorted.bam ${bam} \\
    && samtools index ${bam.baseName}.sorted.bam
    """
}