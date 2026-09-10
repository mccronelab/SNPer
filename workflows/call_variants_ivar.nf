include { CONVERT_TSV_COORDS } from '../modules/convert_tsv_coords'
include { IVAR_VARIANTS } from '../modules/ivar_variants'
include { NEXTCLADE } from '../modules/nextclade'

workflow CALL_VARIANTS_IVAR {
    take:
        variant_bams //[meta, bam, bai] for each replicate
        consensus // [sample, segment, consensus] for each (sample, segment)

    main:
        reference_gff = file(params.reference_gff)
        reference_fasta = file(params.reference_fasta)

        // Map reference GFF annotations to consensus genome. `--input-ref` takes one
        // record, so this batches a segment's consensus sequences across samples into a
        // single run -- 8 for influenza, 1 for SARS-CoV-2. Sorted by name so the
        // published MSA's row order doesn't depend on the order upstream tasks finish in.
        nextclade = consensus
          | map { _sample, segment, consensus_fa -> [segment, consensus_fa] }
          | groupTuple(by:0, sort: { it.name })
          | map { segment, consensus_fas -> [segment, consensus_fas, reference_fasta, reference_gff] }
          | NEXTCLADE

        // Re-key the per-segment GFF batch back onto (sample, segment). A consensus and
        // its GFF are both named for the consensus record ID, so they join on basename;
        // failOnMismatch turns any drift into an error rather than a dropped sample.
        gff_by_record = nextclade.segment_gffs
          | flatMap { _segment, gffs ->
                (gffs instanceof Collection ? gffs : [gffs]).collect { gff -> [gff.baseName, gff] }
            }

        per_consensus_gff = consensus
          | map { sample, segment, consensus_fa -> [consensus_fa.baseName, sample, segment] }
          | join(gff_by_record, by:0, failOnMismatch: true, failOnDuplicate: true)
          | map { _record_id, sample, segment, gff -> [sample, segment, gff] }

        // drop index files and call variants
        variants = variant_bams
          | map { meta, bam, bai -> [meta.sample, meta.segment, meta, bam, bai] }
          | combine(consensus, by:[0,1])
          | combine(per_consensus_gff, by:[0,1])
          | map { _sample, _segment, meta, bam, _bam_index, consensus_fa, gff -> [meta, bam, consensus_fa, gff] }
          | IVAR_VARIANTS // [meta, variants_tsv]
          | map { meta, variants_tsv -> [meta.sample, meta.segment, meta, variants_tsv] }

        // Lift variant coordinates back onto the reference off the alignment NEXTCLADE
        // already produced, rather than aligning a second time. Both of its outputs are one
        // item per segment, so they join on segment and fan back out across that segment's
        // samples.
        segment_alignment = nextclade.msa
          | join(nextclade.nextclade_tsv, by:0, failOnMismatch: true, failOnDuplicate: true)

        reference_coordinate_variants = consensus
          | combine(variants, by:[0,1]) // [sample, segment, consensus, meta, variant_tsv]
          | map { _sample, segment, consensus_fa, meta, variant_tsv -> [segment, meta, consensus_fa, variant_tsv] }
          | combine(segment_alignment, by:0) // [segment, meta, consensus, variant_tsv, msa, nextclade_tsv]
          | map { _segment, meta, consensus_fa, variant_tsv, msa, nextclade_tsv -> [meta, consensus_fa, variant_tsv, msa, nextclade_tsv] }
          | CONVERT_TSV_COORDS

    // emit:
        reference_coordinate_variants
}