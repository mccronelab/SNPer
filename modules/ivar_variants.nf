process IVAR_VARIANTS {
    publishDir "${params.output_dir}/variants/", mode: 'copy'

    cpus 1
    memory 2G
    time 12.h

    input:
        tuple val(key), path(bam) , path(consensus), path(reference_gff)
    output:
        tuple val(key), path("*tsv")
    script:
    """
    samtools sort ${bam} \
    | samtools mpileup -aa -A -d 100000 -B -Q 0 - \
    | ivar variants -p ${bam.simpleName}.variants -q ${params.variant_min_mapQ} -t ${params.variant_freq_threshold} -r ${consensus} -g ${reference_gff}
    """
}