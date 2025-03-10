process MASK_PRIMERS {
    input:
        tuple val(key), path(primer_bed), path(primer_variants_tsv)
        path (primer_pair_tsv)

    output:
        tuple val(key), path(primer_bed), path("${primer_variants_tsv.simpleName}_masked_primers.txt")

    script:
    """
    ivar getmasked -i ${primer_variants_tsv} -b ${primer_bed}  -f ${primer_pair_tsv} -p ${primer_variants_tsv.simpleName}_masked_primers
    """
}