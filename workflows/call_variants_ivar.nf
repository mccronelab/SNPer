include { CONVERT_TSV_COORDS } from '../modules/convert_tsv_coords'
include { IVAR_VARIANTS } from '../modules/ivar_variants'
include { LIFTOFF } from '../modules/liftoff'

workflow CALL_VARIANTS_IVAR {
    take:
        bams_with_consensus // [key, bam, index, consensus] for each replicate

    main:
        reference_gff = file(params.reference_gff)
        reference_fasta = file(params.reference_fasta)
        variants = channel.empty()

        // check file is at least 1Kb in size
        filtered_bams = bams_with_consensus.filter { _key, _bam, _index, consensus -> consensus.size() >= 1000 }
        consensus_fastas = filtered_bams.map{ key, _bam, _index, consensus -> tuple(key, consensus)}.unique{d -> d[0] }

        // map reference GFF annotations to consensus genome
        per_consensus_gff = consensus_fastas.map {key, consensus -> tuple(key, reference_fasta, consensus, reference_gff) }
          | LIFTOFF

        // drop index files and call variants
        filtered_bams.combine(per_consensus_gff, by:0)
          | map { key, bam, _bam_index, consensus, gff -> tuple(key, bam, consensus, gff) }
          | IVAR_VARIANTS
          | set { variants } // [key, variant_tsv]

        reference_coordinate_variants = consensus_fastas.map { key, consensus -> tuple(key, consensus, reference_fasta) }
          | combine(variants, by: 0) // [key, consensus, reference, variant_tsv]
          | CONVERT_TSV_COORDS

    emit:
        reference_coordinate_variants
}