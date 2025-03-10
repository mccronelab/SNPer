process IVAR_PRIMER_VARIANTS {
    input:
        tuple val(key), path(consensus), path(bam)

    output:
        tuple val(key), path("${consensus.simpleName}_primers.bed"), path("${consensus.simpleName}_primer_variants.tsv")

    script:
        """
        samtools mpileup -A -d 0 --reference ${consensus} -Q 0 ${bam} \\
        | ivar variants -p ${consensus}_primer_variants -t 0.02

        bedtools bamtobed -i ${bam} > ${bam.simpleName}_primers.bed
        """
}