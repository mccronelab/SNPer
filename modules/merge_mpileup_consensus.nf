
process MERGE_MPILEUP_CONSENSUS {
    label 'process_high'
    // errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    // Aliased as both ROUGH_CONSENSUS and POLISH_CONSENSUS. Only the polished round
    // publishes; the rough sequence stays internal as the remap target.
    publishDir "${params.output_dir}/consensus_seqs/", mode: 'copy', pattern: "*.fa",
        saveAs: { fn -> publish ? fn : null }
    tag "${sample}${segment_label}${suffix}"

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(sample), val(segment), val(segment_label), val(meta), path(bams), path(bais), path(reference), val(suffix), val(publish)

    output:
        tuple val(sample), val(segment), path("${sample}${segment_label}${suffix}.fa")

    script:
    """
    samtools faidx ${reference}
    samtools merge - ${bams} \
    | samtools sort - \
    | samtools mpileup -a -d 0 -Q ${params.consensus_min_baseQ} -q ${params.consensus_min_mapQ} -f ${reference} - \
    | ivar consensus -p ${sample}${segment_label}${suffix} -i ${sample}${segment_label}${suffix} -n N -q ${params.consensus_min_qual_score} -t ${params.consensus_threshold} -m ${params.consensus_min_depth}
    """
}