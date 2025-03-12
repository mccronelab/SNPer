include { ALIGN_FASTA_FILTER_SORT } from '../modules/align_fasta_filter_sort'
include { IVAR_PRIMER_VARIANTS } from '../modules/ivar_primer_variants'
include { MASK_PRIMERS } from '../modules/mask_primers'
include { REMOVE_MASKED_SORT_INDEX } from '../modules/remove_masked_sort_index'

workflow TRIM_AND_MASK {
    take:
        aligned_reads_consensus // key, sorted bam, index, consensus

    main:
        // filter out empty consensus sequences and split into two channels
        // note: the consensus_seqs channel contains duplicates of every consensus,
        // since aligned_reads_consensus has a row for every replicate
        split = aligned_reads_consensus.filter {
            _key, _bam, _index, consensus -> consensus.size() >= 1000
        }.multiMap { 
            key, bam, index, consensus -> 
            consensus_seqs: tuple(key, consensus)
            aligned_reads: tuple(key, bam, index)
        }

        reference_primers_fasta = file(params.primer_fasta)
        primer_pairs_tsv = file(params.primer_pairs)

        consensus_sequence_ref_primers = split.consensus_seqs.map{
                key, consensus -> tuple(key, consensus, reference_primers_fasta)
            }.unique()

        trimmed_with_consensus = ALIGN_FASTA_FILTER_SORT(consensus_sequence_ref_primers) // [key, consensus, sorted_bam]
          | IVAR_PRIMER_VARIANTS // [key, aligned_primer_bed, mismatch_tsv]
          | map { key, bed, tsv -> tuple(key, bed, tsv, primer_pairs_tsv) }
          | MASK_PRIMERS // [key, aligned_primer_bed, mismatch_list_txt]
          | combine ( split.aligned_reads, by: 0 ) // [key, bed, txt, sorted_bam, index]
          | REMOVE_MASKED_SORT_INDEX 
          // call unique() to remove duplicated consensus sequences
          | combine ( split.consensus_seqs.unique(), by: 0 ) // [key, bam, index, consensus]

    emit:
        trimmed_with_consensus
}