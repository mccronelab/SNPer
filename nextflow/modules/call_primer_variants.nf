process CALL_PRIMER_VARIANTS {
    input:
        tuple val(key_1), path(primer_bam), path(primer_bed)
        tuple val(key_2), path(consensus_file)

    output:
        tuple val(key_1),  path("${primer_bam.simpleName}.tsv"), path(primer_bed)

    script:
        """
        samtools mpileup -aa -A -d 100000 --reference ${consensus_file} -Q ${params.variant_min_qual_score} -q ${params.variant_min_mapQ} -F 0 ${primer_bam} | ivar variants -p ${primer_bam.simpleName}.tsv -t ${params.consensus_threshold}
        """
}