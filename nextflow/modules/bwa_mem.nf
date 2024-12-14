process BWA_MEM {
    container = "staphb/bwa:latest"

    input:
        path reference
        val suffix
        tuple val(key), path(reads_1), path(reads_2)

    output:
        tuple val(key), path("*.sam")

    script:
    """
    bwa index ${reference}
    bwa mem -o ${key}${suffix}.sam ${reference} ${reads_1} ${reads_2}
    """
}