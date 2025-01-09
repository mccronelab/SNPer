process BWA_MEM_FILTER_SORT {
    input:
        tuple val(key), path(consensus_seq), path(consensus_bwt)
        path primer_fasta

    output:
        tuple val(consensus_seq.simpleName), path("${consensus_seq.simpleName}_*.bam")

    script:
    // bwa mem -k: minimum seed length
    // bwa mem -T: minimum score to output
    // samtools view -F 4: Exclude unmapped reads from output
    // samtools view -b: output in BAM format
        """
        bwa mem -k 5 -T 16 ${consensus_seq} ${primer_fasta} | samtools view -F 4 -b | samtools sort -o ${consensus_seq.simpleName}_1.bam
        """
}