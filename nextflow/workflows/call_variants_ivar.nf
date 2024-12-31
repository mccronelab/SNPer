include { BWA_SAMTOOLS_INDEX } from "../modules/bwa_samtools_index"
include { BWA_MEM_FILTER_SORT } from "../modules/bwa_mem_filter_sort"
include { BAM_TO_BED } from '../modules/bam_to_bed'
include { CALL_PRIMER_VARIANTS } from '../modules/call_primer_variants'
include { MASK_PRIMERS } from '../modules/mask_primers.nf'

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

        masks = MASK_PRIMERS(mismatches_and_beds, primer_info)

        




    emit:
        variants = called_variants
}