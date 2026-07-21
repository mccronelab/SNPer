
process MERGE_MPILEUP_CONSENSUS {
    label 'process_high'
    // errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    publishDir "${params.output_dir}/consensus_seqs/", mode: 'copy', pattern: "*.fa"
    tag "${sample}_${segment}"

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(sample), val(segment), val(meta), path(bams), path(bais), path(reference), val(suffix)

    output:
        tuple val(sample), val(segment), path("${sample}_${segment}${suffix}.fa")

    script:
    """
    samtools faidx ${reference}
    samtools merge - ${bams} \
    | samtools sort - \
    | samtools mpileup -d 0 -Q 30 -q ${params.variant_min_mapQ} -f ${reference} - \
    | ivar consensus -p ${sample}_${segment}${suffix} -i ${sample}_${segment}${suffix} -n N -q ${params.consensus_min_qual_score} -t ${params.consensus_threshold} -m ${params.consensus_min_depth}
    """
}