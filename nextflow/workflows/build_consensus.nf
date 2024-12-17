// This workflow uses BWA mem, samtools, picard, and iVar to build a consensus genome

include { FASTQC as FASTQC_A } from "../modules/fastqc.nf"
include { FASTQC as FASTQC_B } from "../modules/fastqc.nf"
include { BWA_MEM as BWA_MEM_A } from "../modules/bwa_mem.nf"
include { BWA_MEM as BWA_MEM_B } from "../modules/bwa_mem.nf"
include { FILTER_SORT_INDEX as FILTER_SORT_INDEX_A } from "../modules/filter_sort_index.nf"
include { FILTER_SORT_INDEX as FILTER_SORT_INDEX_B } from "../modules/filter_sort_index.nf"
include { IVAR_TRIM as IVAR_TRIM_A } from "../modules/ivar_trim.nf"
include { IVAR_TRIM as IVAR_TRIM_B } from "../modules/ivar_trim.nf"
include { PICARD_SORT as PICARD_SORT_A } from "../modules/picard_sort.nf"
include { PICARD_SORT as PICARD_SORT_B } from "../modules/picard_sort.nf"
include { GET_COVERAGE as GET_COVERAGE_A } from "../modules/get_coverage.nf"
include { GET_COVERAGE as GET_COVERAGE_B } from "../modules/get_coverage.nf"
include { MERGE_MPILEUP_CONSENSUS } from "../modules/merge_mpileup_consensus.nf"

workflow CONSENSUS_GEN {
    take:
        reference_genome // should be a variable or singleton channel with path to reference workflow
        technical_rep_A
        technical_rep_B
        primer_bed

    main:
        // each fastqc channel tuple contains a key, [path(read1), path(read2)]
        FASTQC_A(technical_rep_A)
        FASTQC_B(technical_rep_B)

        aligned_sam_A = BWA_MEM_A(reference_genome, "_A", technical_rep_A)
        aligned_sam_B = BWA_MEM_B(reference_genome, "_B", technical_rep_B)

        bam_A = FILTER_SORT_INDEX_A(aligned_sam_A)
        bam_B = FILTER_SORT_INDEX_B(aligned_sam_B)

        prim_trim_bam_A = IVAR_TRIM_A(bam_A, primer_bed)
        prim_trim_bam_B = IVAR_TRIM_B(bam_B, primer_bed)

        sorted_prim_trim_bam_A = PICARD_SORT_A(prim_trim_bam_A)
        sorted_prim_trim_bam_B = PICARD_SORT_B(prim_trim_bam_B)

        GET_COVERAGE_A(sorted_prim_trim_bam_A)
        GET_COVERAGE_B(sorted_prim_trim_bam_B)

        paired_replicates = sorted_prim_trim_bam_A.join(sorted_prim_trim_bam_B)
        paired_replicates.view()

        consensus_sequence = MERGE_MPILEUP_CONSENSUS(reference_genome, paired_replicates)

    emit:
        consensus = consensus_sequence
}