include { CONVERT_TSV_COORDS } from '../modules/convert_tsv_coords'
include { IVAR_VARIANTS } from '../modules/ivar_variants'
include { LIFTOFF } from '../modules/liftoff'

workflow CALL_VARIANTS_IVAR {
    take:
        variant_bams //[meta, bam, bai] for each replicate
        consensus // [sample, consensus] for each sample

    main:
        reference_gff = file(params.reference_gff)
        reference_fasta = file(params.reference_fasta)

        // map reference GFF annotations to consensus genome
        per_consensus_gff = consensus.map { sample, consensus_fa -> [sample, reference_fasta, consensus_fa, reference_gff] }
          | LIFTOFF

        // drop index files and call variants
        variants = variant_bams
          | map { meta, bam, bai -> [meta.sample, meta, bam, bai] }
          | combine(consensus, by:0)
          | combine(per_consensus_gff, by:0)
          | map { _sample, meta, bam, _bam_index, consensus_fa, gff -> [meta, bam, consensus_fa, gff] }
          | IVAR_VARIANTS // [meta, variants_tsv]
          | map { meta, variants_tsv -> [meta.sample, meta, variants_tsv] }

        reference_coordinate_variants = consensus.map { sample, consensus_fa -> [sample, consensus_fa, reference_fasta] }
          | combine(variants, by: 0) // [meta, consensus, reference, variant_tsv]
          | map { _sample, consensus_fa, reference_fa, meta, variant_tsv -> [meta, consensus_fa, reference_fa, variant_tsv] }
          | CONVERT_TSV_COORDS

    // emit:
        reference_coordinate_variants
}