process IVAR_VARIANTS {
    label 'process_high'
    // retry if error message indicates a failure due to resource limits
    //errorStrategy { task.exitStatus in 1;137..140 ? 'retry' : 'terminate' }
    //stageInMode 'copy'

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }
    tag "${meta.replicate_id}"

    input:
        tuple val(meta), path(bam), path(consensus), path(gff)

    output:
        tuple val(meta), path("*tsv")
        
    script:
    """
    {
        samtools faidx ${consensus}
        samtools sort ${bam} \
        | samtools mpileup -d 0 -Q 30 -q ${params.variant_min_mapQ}  -f ${consensus} - \
        | ivar variants -p ${bam.simpleName}.variants -q ${params.variant_minQ} -m ${params.variant_min_depth} -t ${params.variant_freq_threshold} -r ${consensus} -g ${gff}
    } 2> var_errors.txt

    if
        grep -qF "[E::faidx_adjust_position]" var_errors.txt; then
        echo "Catch '[E::faidx_adjust_position]' error in variant calling"
        exit 1
    fi
    """
}