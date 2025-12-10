include { CONVERT_TSV_COORDS } from '../modules/convert_tsv_coords'
include { IVAR_VARIANTS } from '../modules/ivar_variants'
include { LIFTOFF } from '../modules/liftoff'

workflow TEST_VARIANTS {
    take:
        bams_with_r_consensus

    main:
        reference_fa = file(params.reference_fasta)
        reference_gff = file(params.reference_gff)

        per_sample_gff = bams_with_r_consensus
            | map { meta, _bam, r_consensus -> [meta, reference_fa, r_consensus, reference_gff] }
            | LIFTOFF

        raw_variants = bams_with_r_consensus
            | map { meta, bam, r_consensus -> [meta, bam, r_consensus, reference_gff] }
            | IVAR_VARIANTS
            | map { meta, tsv, _mpileup -> [meta, tsv] }

        reference_coordinate_variants = bams_with_r_consensus
            | map {meta, _bam, r_consensus -> [meta, r_consensus, reference_fa] }
            | combine(raw_variants, by:0)
            | CONVERT_TSV_COORDS
}