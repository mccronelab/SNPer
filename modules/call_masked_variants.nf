process CALL_MASKED_VARIANTS {
    publishDir "${params.output_dir}/called_variant_tsvs/", mode: 'copy'

    input:
        tuple val(key), path(masked_bam), path(masked_index), path (reference_fasta)
        path(reference_gff)

    output:
        tuple val(key), path("*.tsv")

    script:
    """
    bwa index ${reference_fasta}
    samtools mpileup -aa -A -d 100000 -B -Q 0 -q ${params.variant_min_mapQ} --reference ${reference_fasta} ${masked_bam} | ivar variants -p ${masked_bam.simpleName}.variants -q ${params.variant_min_mapQ} -t ${params.variant_freq_threshold} -r ${reference_fasta} -g ${reference_gff}
    """
}