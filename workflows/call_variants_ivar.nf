include { CONVERT_TSV_COORDS } from '../modules/convert_tsv_coords'
include { IVAR_VARIANTS } from '../modules/ivar_variants'
include { LIFTOFF } from '../modules/liftoff'

workflow CALL_VARIANTS_IVAR {
    take:
        // bams_with_consensus // [meta, bam, index, consensus] for each replicate
        variant_bams
        rough_consensus
        polished_consensus

    main:
        reference_gff = file(params.reference_gff)
        reference_fasta = file(params.reference_fasta)
        variants = channel.empty()

        polished_consensus.view()

        // check file is at least 1Kb in size
        // filtered_bams = bams_with_consensus.filter{ _meta, _bam, _index, consensus -> consensus.size() >= 1000 }
        // consensus_fastas = filtered_bams.map{ meta, _bam, _index, consensus -> tuple(meta, consensus)}.unique{d -> d[0] }

        // map reference GFF annotations to consensus genome
        per_consensus_gff = polished_consensus.map {meta, consensus -> tuple(meta, reference_fasta, consensus, reference_gff) }
          | LIFTOFF

        // drop index files and call variants
        variant_bams
          | combine(rough_consensus, by:0)
          | combine(polished_consensus, by:0)
          | combine(per_consensus_gff, by:0)
          | map { meta, bam, _bam_index, r_consensus, p_consensus, gff -> tuple(meta, bam, r_consensus, p_consensus, gff) }
          | IVAR_VARIANTS
          | set { variants } // [meta, variant_tsv]

        reference_coordinate_variants = polished_consensus.map { meta, consensus -> tuple(meta, consensus, reference_fasta) }
          | combine(variants, by: 0) // [meta, consensus, reference, variant_tsv]
          | CONVERT_TSV_COORDS

    // emit:
        reference_coordinate_variants
}