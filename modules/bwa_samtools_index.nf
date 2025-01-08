process BWA_SAMTOOLS_INDEX {
    input:
        tuple val(key), path(consensus_seq)

    output:
        tuple val(consensus_seq.simpleName), path("${consensus_seq.simpleName}.fa.*")

    script:
        """
        bwa index ${consensus_seq}
        samtools faidx ${consensus_seq}
        """
}