process IVAR_VARIANTS {
    publishDir "${params.output_dir}/variants/", mode: 'copy'
    // retry if error message indicates a failure due to resource limits
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(key), path(bam), path(consensus), path(gff)

    output:
        tuple val(key), path("*tsv")
        
    script:
    """
    samtools sort ${bam} \
    | samtools mpileup -aa -A -d 100000 -B -Q 0 -q ${params.variant_min_mapQ}  --reference ${consensus} - \
    | ivar variants -p ${bam.simpleName}.variants -q ${params.variant_minQ} -m ${params.variant_min_depth} -t ${params.variant_freq_threshold} -r ${consensus} -g ${gff}
    """
}