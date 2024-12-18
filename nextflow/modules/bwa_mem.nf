process BWA_MEM {
    container "staphb/bwa:latest"

    input:
        path reference
        val suffix
        tuple val(key), path(paired_reads)

    output:
        tuple val(key), path("*.sam")

    script:
        """
        bwa index ${reference}
        bwa mem -o ${key}${suffix}.sam ${reference} ${paired_reads}
        """
}