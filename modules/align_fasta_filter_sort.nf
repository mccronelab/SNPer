process ALIGN_FASTA_FILTER_SORT {
    cpus 1
    memory 2G
    time 12.h 

    input:
        tuple val(key), path(consensus), path(fasta)

    output:
        tuple val(key), path("${consensus.simpleName}.fasta.bam"), path(consensus)

    script:
    """
    bwa mem -k 5 -T 16 ${consensus} ${fasta} \\
    | samtools view -bS -F 4 \\
    | samtools sort -o ${consensus.simpleName}.fasta.bam
    """
}