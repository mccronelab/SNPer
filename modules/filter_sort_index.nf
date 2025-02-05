process FILTER_SORT_INDEX {
    cpus 1
    memory 2G
    time 12.h

    input:
        tuple val(key), path(sam_file)

    output:
        tuple val(key), path("*.sorted.bam"), path("*.bai")

    script:
    // samtools view -F 4: Exclude unmapped reads from output
    // samtools view -b: output in BAM format
    """
    samtools view -F 4 -b ${sam_file} \\
    | samtools sort -o ${sam_file.baseName}.sorted.bam \\
    && samtools index ${sam_file.baseName}.sorted.bam
    """
}