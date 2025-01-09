process BWA_MEM {

    input:
        tuple val(key), path(paired_reads), path(reference)

    output:
        tuple val(key), path("*.sam")

    script:
        """
        bwa index ${reference}
        bwa mem -o ${paired_reads[1].simpleName.split("_")[0]}.sam ${reference} ${paired_reads} # assumes nothing after _ is important for meta data
        """
}