process IVAR_VARIANTS {
    label 'process_high'
    // retry if error message indicates a failure due to resource limits
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    publishDir "${params.output_dir}/raw_variants/", mode: 'copy', pattern: "*.tsv"

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(meta), path(bam), path(consensus), path(gff)

    output:
        tuple val(meta), path("*tsv")
    
    script:
    def sample = meta.sample

    """
    samtools sort ${bam} \
    | samtools mpileup -d 0 -Q 30 -q ${params.variant_min_mapQ}  -f ${consensus} - > ${sample}.mpileup 2> ${sample}.mpileup.err
    
    cat ${sample}.mpileup \
    | ivar variants -p ${bam.simpleName}.variants -q ${params.variant_minQ} -m ${params.variant_min_depth} -t ${params.variant_freq_threshold} -r ${consensus} -g ${gff} 2> ivar.err
    """
}