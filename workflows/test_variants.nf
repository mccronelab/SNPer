include { IVAR_VARIANTS } from '../modules/ivar_variants'
include { CONVERT_TSV_COORDS } from '../modules/convert_tsv_coords'

workflow TEST_VARIANTS {
    take:
        bams_with_r_consensus
        polished_consensus

    main:
        reference_fa = file(params.reference_fasta)
        reference_gff = file(params.reference_gff)

        raw_variants = bams_with_r_consensus
            | map { meta, bam, r_consensus -> [meta, bam, r_consensus, reference_gff] }
            | IVAR_VARIANTS

        reference_coordinate_variants = bams_with_r_consensus
            | map {meta, _bam, r_consensus -> [meta, r_consensus, reference_fa] }
            | combine(raw_variants, by:0)
            | CONVERT_TSV_COORDS
}