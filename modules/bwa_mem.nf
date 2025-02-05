process BWA_MEM {
    cpus 1
    memory 2G
    time 12.h


    input:
        tuple val(key), path(paired_reads), path(reference)

    output:
        tuple val(key), path("*.sam")

    script:
        """
        bwa index ${reference}
        bwa mem -o ${key}_${paired_reads[1].simpleName.split("_")[1]}.sam ${reference} ${paired_reads} # assumes nothing after _ is important for meta data
        """
}