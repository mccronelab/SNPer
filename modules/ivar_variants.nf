process IVAR_VARIANTS {
    label 'process_high'
    // retry if error message indicates a failure due to resource limits
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(meta), path(bam), path(rough_consensus, stageAs: 'rough_fa/*'), path(polished_consensus), path(gff)

    output:
        tuple val(meta), path("*tsv")
        
    script:
    """
    samtools sort ${bam} \
    | samtools mpileup -aa -d 0 -Q 30 -q ${params.variant_min_mapQ}  -f ${rough_consensus} - \
    | ivar variants -p ${bam.simpleName}.variants -q ${params.variant_minQ} -m ${params.variant_min_depth} -t ${params.variant_freq_threshold} -r ${polished_consensus} -g ${gff}
    """
}