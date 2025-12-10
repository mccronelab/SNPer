include { BWA_MEM as BWA_MEM_REF; BWA_MEM as BWA_MEM_ROUGH; BWA_MEM as BWA_MEM_POL;
    BWA_REMAP as BWA_REMAP_CON; BWA_REMAP as BWA_REMAP_POL } from "../modules/bwa_mem"
include { SORT_INDEX_BAM as SORT_INDEX_ROUGH; SORT_INDEX_BAM as SORT_INDEX_REF
    SORT_INDEX_BAM as SORT_INDEX_POLISH } from "../modules/sort_index_bam"
include { AMPLICON_CLIP } from "../modules/amplicon_clip"
include { DEDUPLICATE_READS as DEDUP_MIPS_REF; DEDUPLICATE_READS as DEDUP_MIPS_ROUGH;
    DEDUPLICATE_READS as DEDUP_POL } from "../modules/deduplicate_reads"
include { PICARD_SORT as PS_CON; PICARD_SORT as PS_VAR  } from "../modules/picard_sort"
include { GET_CONSENSUS_COVERAGE; GET_VARIANT_READ_DEPTH  } from "../modules/get_coverage"
include { MERGE_MPILEUP_CONSENSUS as ROUGH_CONSENSUS; 
    MERGE_MPILEUP_CONSENSUS as POLISHED_CONSENSUS } from "../modules/merge_mpileup_consensus"

workflow TEST_CONSENSUS {
    take:
        samples

    main:
        reference = file(params.reference_fasta)

        ref_ali_bam = samples.map { meta, reads -> [meta, reads, reference] } 
            | BWA_MEM_REF  // (meta, path(*sam))

        dedup_ref_bam = ref_ali_bam
            | DEDUP_MIPS_REF

        rough_consensus = dedup_ref_bam
            | map { meta, bam, bai -> [meta, bam, bai, reference, "_rough"] }
            | ROUGH_CONSENSUS 
            // filter out consensus with no sequence, which sometimes occurs
            | filter { _meta, consensus_fa -> consensus_fa
                if (!consensus_fa.exists() || consensus_fa.size() == 0) {
                    return false
                }
                def seq = consensus_fa.splitFasta(record: [seqString: true])
                if (seq) {
                    return true
                }
            }

        rough_ali_bam = samples
            | combine(rough_consensus, by:0)
            | BWA_MEM_ROUGH
            | DEDUP_MIPS_ROUGH

        polished_consensus = rough_ali_bam
            | combine(rough_consensus, by:0)
            | map { meta, bam, bai, r_consensus -> [meta, bam, bai, r_consensus, "_polished"] }
            | POLISHED_CONSENSUS

        polished_ali_bam = samples
            | combine(polished_consensus, by:0)
            | BWA_MEM_POL
            | DEDUP_POL

        GET_VARIANT_READ_DEPTH(polished_ali_bam)

        polished_consensus
            | GET_CONSENSUS_COVERAGE
            | filter { _sample, _consensus, coverage -> Float.parseFloat(coverage)>= params.consensus_coverage_cutoff }
            | map { sample, consensus, _coverage -> [sample, consensus] }

        polished_ali_with_p_consensus = polished_ali_bam
            | combine(polished_consensus, by:0)
            | map { meta, bam, _bai, r_con -> [meta, bam, r_con] }

    emit:
        bams_with_consensus = polished_ali_with_p_consensus
}

