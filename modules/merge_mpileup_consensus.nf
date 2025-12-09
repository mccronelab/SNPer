
process MERGE_MPILEUP_CONSENSUS {
    label 'process_high'
    // errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    publishDir "${params.output_dir}/consensus_seqs/rough", mode: 'copy', pattern: "*_rough.fa"
    publishDir "${params.output_dir}/consensus_seqs/polished", mode: 'copy', pattern: "*_polished.fa"

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(meta), path(bams), path(bais), path(reference), val(suffix)

    output:
        tuple val(meta), path("${sample}${suffix}.fa")

    script:
    sample = meta.sample

    """
    samtools faidx ${reference}
    samtools merge - ${bams} \
    | samtools sort - \
    | samtools mpileup -d 0 -Q 30 -q ${params.variant_min_mapQ} -f ${reference} - \
    | ivar consensus -p ${sample}${suffix} -i ${sample}${suffix} -n N -q ${params.consensus_min_qual_score} -t ${params.consensus_threshold} -m ${params.consensus_min_depth}
    """
}