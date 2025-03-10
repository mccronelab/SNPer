process SAMTOOLS_SORT_AND_INDEX {
    cpus 1
    memory 20G
    time 12.h

    input:
        tuple val(key), path(trimmed_bam)

    output:
        tuple val(key), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")

    script:
    """
    samtools sort -o ${trimmed_bam.simpleName}.sorted.bam ${trimmed_bam}
    samtools index ${trimmed_bam.simpleName}.sorted.bam
    """
}