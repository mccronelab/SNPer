include { IVAR_VARIANTS } from '../modules/ivar_variants.nf'
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
        consensus_fastas = filtered_bams.map{key, _bam, _index, consensus -> tuple(key, consensus)}.unique{d -> d[0]}

        // map reference GFF annotations to consensus genome
        per_consensus_gff = consensus_fastas.map {key, consensus -> tuple(key, reference_fasta, consensus, reference_gff) }
          | LIFTOFF

        filtered_bams.join{ per_consensus_gff }
          | IVAR_VARIANTS
          | set{ variants }

    emit:
        variants 
}