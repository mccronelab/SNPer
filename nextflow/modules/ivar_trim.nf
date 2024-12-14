process IVAR_TRIM {
    container = "staphb/ivar:latest"

    input:
        tuple val(key), path(sorted_bam), path(bam_index)
        path(primer_bedfile)

    output:
        tuple val(key), path("*.primertrim.bam")

    script:
    """
    ivar trim -i ${sorted_bam} -b ${primer_bedfile} -p ${sorted_bam.baseName}.primertrim.bam
    """
}