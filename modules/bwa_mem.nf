process BWA_MEM {
    errorStrategy 'retry'
    maxRetries 3

    cpus {1 * task.attempt}
    // need to account for potentially increasing CPU allocation
    memory { 2G * task.attempt * task.cpus}
    time { 4.h * task.attempt }


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