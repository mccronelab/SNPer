process BAM_TO_BED {
    input:
        tuple val(key), path(bam_file)

    output:
        tuple val(key), path("${bam_file.simpleName}.bed")

    script:
        """
        bedtools bamtobed -i ${bam_file} > ${bam_file.simpleName}.bed
        """
}