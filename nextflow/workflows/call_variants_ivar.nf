include { BWA_SAMTOOLS_INDEX } from "../modules/bwa_samtools_index"
include { BWA_ALIGN_FILTER_SORT } from "../modules/bwa_align_filter_sort"

workflow CALL_VARIANTS_IVAR {
    take:
        reference_fasta
        reference_gff
        primer_bed
        primer_pairs
        primer_fasta
        consensus_seq

    main:
        // tuple (consensus_name, BWT_path)
        consensus_bwts = BWA_SAMTOOLS_INDEX(consensus_seq)

        indexed_seqs = consensus_seq.join(consensus_bwts)

        primer_bams = BWA_ALIGN_FILTER_SORT(indexed_seqs, primer_fasta)

        called_variants = channel.empty()
    emit:
        variants = called_variants
}