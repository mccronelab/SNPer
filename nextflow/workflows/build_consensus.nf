// This workflow uses BWA mem, samtools, and iVar to build a consensus genome

include { FASTQC as FASTQC_A } from "../modules/fastqc.nf"
include { FASTQC as FASTQC_B } from "../modules/fastqc.nf"
include { BWA_MEM as BWA_MEM_A } from "${projectDir}/modules/bwa_mem.nf"
include { BWA_MEM as BWA_MEM_B } from "${projectDir}/modules/bwa_mem.nf"
include { FILTER_SORT_INDEX as FILTER_SORT_INDEX_A } from "${projectDir}/modules/filter_sort_index.nf"
include { FILTER_SORT_INDEX as FILTER_SORT_INDEX_B } from "${projectDir}/modules/filter_sort_index.nf"
include { IVAR_TRIM as IVAR_TRIM_A } from "../modules/ivar_trim.nf"
include { IVAR_TRIM as IVAR_TRIM_B } from "../modules/ivar_trim.nf"

workflow CONSENSUS_GEN {
    take:
        reference_genome // should be a variable or singleton channel with path to reference workflow
        technical_rep_A
        techincal_rep_B
        primer_bed

    main:
    // each fastqc channel tuple contains a key, [path(read1), path(read2)]
    fastqc_A = FASTQC_A(technical_rep_A)
    fastqc_B = FASTQC_B(techincal_rep_B)

    aligned_sam_A = BWA_MEM_A(reference_genome, "_A", fastqc_A)
    aligned_sam_B = BWA_MEM_B(reference_genome, "_B", fastqc_B)

    bam_A = FILTER_SORT_INDEX_A(aligned_sam_A)
    bam_B = FILTER_SORT_INDEX_B(aligned_sam_B)

    prim_trim_bam_A = IVAR_TRIM_A(bam_A, primer_bed)
    prim_trim_bam_B = IVAR_TRIM_B(bam_B, primer_bed)


    consensus_sequence = channel.empty()


    emit:
        consensus = consensus_sequence
}