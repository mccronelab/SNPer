process MASK_PRIMERS {
    input:
        tuple val(key), path(mismatch_tsv), path(primer_bed)
        path (primer_info_tsv)

    output:
        tuple val(key), path("${mismatch_tsv.simpleName}_masked_primers.txt")

    script:
    """
    ivar getmasked -i ${mismatch_tsv} -b ${primer_bed}  -f ${primer_info_tsv} -p ${mismatch_tsv.simpleName}_masked_primers
    """
}