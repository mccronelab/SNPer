process IVAR_PRIMER_VARIANTS {
    input:
        tuple val(key), path(bam), path(consensus)

    output:
        tuple val(key), path("${consensus.simpleName}_primers.bed"), path("${consensus.simpleName}_primer_variants.tsv")

    script:
        """
        samtools mpileup -aa -A -d 10000 --reference ${consensus} -Q ${params.variant_minQ}  -q ${params.variant_min_mapQ} -F 0  ${bam} \\
        | ivar variants -p ${consensus.simpleName}_primer_variants -t ${params.variant_freq_threshold}

        bedtools bamtobed -i ${bam} > ${bam.simpleName}_primers.bed
        """
}