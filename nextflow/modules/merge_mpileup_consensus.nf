process MERGE_MPILEUP_CONSENSUS {
    publishDir "${params.output_dir}/consensus_seqs/", mode: 'copy'

    input:
        path reference_genome
        tuple val(key), path(bam_A), path(bam_A_index), path(bam_B), path(bam_B_index)

    output:
        tuple val(key), path("${key}.fa")

    script:
    """
    samtools merge ${key}_merged.bam ${bam_A} ${bam_B}
    samtools mpileup -a -A -d 100000 -Q 0 --reference ${reference_genome} ${key}_merged.bam | ivar consensus -p ${key}.fasta -n N -q ${params.min_qual_score} -t ${params.consensus_threshold} -m ${params.min_depth}
    """
}