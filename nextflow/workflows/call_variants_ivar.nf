include { BWA_SAMTOOLS_INDEX } from "../modules/bwa_samtools_index"
include { BWA_MEM_FILTER_SORT } from "../modules/bwa_mem_filter_sort"
include { BAM_TO_BED } from '../modules/bam_to_bed'
include { CALL_PRIMER_VARIANTS } from '../modules/call_primer_variants'
include { MASK_PRIMERS } from '../modules/mask_primers.nf'
include { REMOVE_MASKED_SORT_INDEX } from '../modules/remove_masked_sort_index'
include { CALL_MASKED_VARIANTS } from '../modules/call_masked_variants'

workflow CALL_VARIANTS_IVAR {
    take:
        reference_fasta
        reference_gff
        primer_bed
        primer_pairs
        primer_fasta
        primer_info
        consensus_seq

    main:
        // TODO: delete me when workflow completed
        called_variants = channel.empty()

        // tuple (consensus_name, [a bunch of index files])
        consensus_bwts = BWA_SAMTOOLS_INDEX(consensus_seq)

        // tuple (key: consensus_name, consensus.fasta, [a bunch of index files])
        indexed_seqs = consensus_seq.join(consensus_bwts)

        // create primer bams and beds to identify primer variants
        // tuple (key: consensus_name, primer.bam)
        primer_bams = BWA_MEM_FILTER_SORT(indexed_seqs, primer_fasta)
        // tuple (key: consensus_name, primer.bed)
        primer_beds = BAM_TO_BED(primer_bams)
        // tuple (key: consensus_name, primer.bam, primer.bed)
        bams_and_beds = primer_bams.join(primer_beds)

        // call and mask primer variants
        // tuple (key: consensus_name, mismatches.tsv, primer_bed)
        mismatches_and_beds = CALL_PRIMER_VARIANTS(bams_and_beds, consensus_seq)

        // tuple (key: consensus_name, masked_primers.txt)
        masks = MASK_PRIMERS(mismatches_and_beds, primer_info)
        // tuple (key: consensus_name, masked_primers.txt, primers.bam)
        masks_and_bams = masks.join(primer_bams)
        // tuple (key: consensus_name, masked.bam, masked.bam.bai)
        masked_reads_removed = REMOVE_MASKED_SORT_INDEX(masks_and_bams, primer_bed)
        masked_reads_removed_consensus = masked_reads_removed.join(consensus_seq)

        // tuple (key: consensus_name, )
        called_variants = CALL_MASKED_VARIANTS(masked_reads_removed_consensus, reference_gff)

    emit:
        variants = called_variants
}