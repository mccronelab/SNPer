// This workflow uses BWA mem, samtools, picard, and iVar to build a consensus genome

include { FASTQC} from "../modules/fastqc"
include { BWA_MEM  } from "../modules/bwa_mem"
include { FILTER_SORT_INDEX  } from "../modules/filter_sort_index"
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
        |   BWA_MEM  // (key, path(*sam))
        |   FILTER_SORT_INDEX // tuple val(key), path("*.sorted.bam"), path("*.bai")

        
        polished_bam = bam.map { key, sortedBam, bamIndex -> tuple(key,sortedBam,bamIndex,primer_bed)} 
        | IVAR_TRIM  //  tuple val(key), path("*.primertrim.bam")
        | PICARD_SORT // tuple val(key), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")

        //TODO remove dups?

        GET_COVERAGE(polished_bam)

        paired_replicates = polished_bam
                            .groupTuple() // [key, [bams], [bais]]


        consensus_sequence = MERGE_MPILEUP_CONSENSUS(paired_replicates)

    emit:
        consensus = consensus_sequence // (key, consensus.fa)
}