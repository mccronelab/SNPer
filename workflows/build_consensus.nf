// This workflow uses BWA mem, samtools, picard, and iVar to build a consensus genome

include { FASTQC} from "../modules/fastqc"
include { BWA_MEM as BWA_MEM_CON  } from "../modules/bwa_mem"
include { BWA_MEM as BWA_MEM_VAR } from "../modules/bwa_mem"
include { FILTER_SORT_INDEX as FSI_CON  } from "../modules/filter_sort_index"
include { FILTER_SORT_INDEX as FSI_VAR  } from "../modules/filter_sort_index"
include { IVAR_TRIM  } from "../modules/ivar_trim"
include { PICARD_SORT  } from "../modules/picard_sort"
include { GET_COVERAGE  } from "../modules/get_coverage"
include { MERGE_MPILEUP_CONSENSUS } from "../modules/merge_mpileup_consensus"

workflow CONSENSUS_GEN {
    take:
        samples // tuple (key,[fastq1,fastq2])

    main:
        primer_bed = file(params.primer_bed)
        reference = file(params.reference_fasta)

        // each fastqc channel tuple contains a key, [path(read1), path(read2)]
        FASTQC(samples)
        
        // TODO add fastp step

        bam = samples.map { key, reads -> tuple(key, reads, reference) } 
          | BWA_MEM_CON  // (key, path(*sam))
          | FSI_CON // tuple val(key), path("*.sorted.bam"), path("*.bai")

        
        polished_bam = bam.map { key, sortedBam, bamIndex -> tuple(key, sortedBam, bamIndex, primer_bed)} 
          | IVAR_TRIM  //  tuple val(key), path("*.primertrim.bam")
          | PICARD_SORT // tuple val(key), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")

        consensus_sequence = polished_bam.groupTuple() // [key, [bams], [bais]]
          | MERGE_MPILEUP_CONSENSUS // key, consensus

        // filter out empty consensus sequences
        consensus_sequence = consensus_sequence.filter { _key, consensus -> consensus.size() >= 1000 }

        variant_bam = samples.combine(consensus_sequence, by:0) // key, [reads], consensus
          | BWA_MEM_VAR // key, sams
          | FSI_VAR // key, sorted bam, index

        GET_COVERAGE(variant_bam)

        // BWA_MEM() doesn't output the consensus, so we rejoin it
        variant_bam_consensus = variant_bam.combine(consensus_sequence, by:0) // key, sorted bam, index, consensus

    emit:
        reads_and_consensus = variant_bam_consensus // (key, sorted bam, bam index, consensus.fa)
}