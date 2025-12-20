include { CONVERT_TSV_COORDS } from '../modules/convert_tsv_coords'
include { IVAR_VARIANTS } from '../modules/ivar_variants'
include { LIFTOFF } from '../modules/liftoff'

workflow CALL_VARIANTS_IVAR {
    take:
        bams_with_consensus // [meta, bam, index, consensus] for each replicate

    main:
        reference_gff = file(params.reference_gff)
        reference_fasta = file(params.reference_fasta)
        variants = channel.empty()

        // check file is at least 1Kb in size
        consensus_fastas = bams_with_consensus.map{ meta, _bam, _index, consensus -> tuple(meta, consensus)}.unique{d -> d[0] }

        // map reference GFF annotations to consensus genome
        per_consensus_gff = consensus_fastas.map { meta, consensus -> tuple(meta, reference_fasta, consensus, reference_gff) }
          | LIFTOFF

        // drop index files and call variants
        bams_with_consensus.combine(per_consensus_gff, by:0)
          | map { meta, bam, _bam_index, consensus, gff -> tuple(meta, bam, consensus, gff) }
          | IVAR_VARIANTS
          | set { variants } // [meta, variant_tsv]

        reference_coordinate_variants = consensus_fastas.map { meta, consensus -> tuple(meta, consensus, reference_fasta) }
          | combine(variants, by: 0) // [meta, consensus, reference, variant_tsv]
          | CONVERT_TSV_COORDS

    // emit:
        reference_coordinate_variants
}