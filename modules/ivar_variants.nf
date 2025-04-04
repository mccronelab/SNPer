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
    | samtools mpileup -aa -A -d 100000 -Q 0  --reference ${consensus} - \
    | ivar variants -p ${bam.simpleName}.variants -q ${params.variant_min_mapQ} -t ${params.variant_freq_threshold} -r ${consensus} -g ${gff}
    """
}