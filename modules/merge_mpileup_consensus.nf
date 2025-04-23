
process MERGE_MPILEUP_CONSENSUS {
    publishDir "${params.output_dir}/consensus_seqs/", mode: 'copy'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(key), path(bams), path(bais)

    output:
        tuple val(key), path("${key}.fa")

    script:
    """
    samtools merge - ${bams} \
    | samtools sort - \
    | samtools mpileup -a -A -d 100000 -Q 0 - \
    | ivar consensus -p ${key}.fasta -n N -q ${params.consensus_min_qual_score} -t ${params.consensus_threshold} -m ${params.consensus_min_depth}
    """
}