process BWA_MEM {
    cpus 1
    memory 2G
    time 12.h


    input:
        tuple val(key), val(replicate_id), path(paired_reads), path(reference)

    output:
        tuple val(key), path("*.sam")

    script:
        """
        bwa index ${reference}
        bwa mem -o ${replicate_id}.sam ${reference} ${paired_reads}
        """
}