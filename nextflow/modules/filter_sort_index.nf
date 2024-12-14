process FILTER_SORT_INDEX {
    container = "staphb/samtools:latest"

    input:
        tuple val(key), path(sam_file)

    output:
        tuple val(key), path("*.sorted.bam"), path("*.bai")

    script:
    """
    samtools view -F 4 -Sb ${sam_file} \\
    | samtools sort -o ${sam_file.baseName}.sorted.bam \\
    && samtools index ${sam_file.baseName}.sorted.bam
    """
}