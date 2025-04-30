// This workflow uses BWA mem, samtools, picard, and iVar to build a consensus genome

include { FASTQC} from "../modules/fastqc"
include { BWA_MEM; BWA_REMAP  } from "../modules/bwa_mem"
include { BWA_MEM as BWA_MEM_VAR } from "../modules/bwa_mem"
include { FILTER_SORT_INDEX as FSI_CON  } from "../modules/filter_sort_index"
include { FILTER_SORT_INDEX as FSI_VAR  } from "../modules/filter_sort_index"
include {AMPLICON_CLIP} from "../modules/amplicon_clip"
include { PICARD_SORT as PS_CON  } from "../modules/picard_sort"
include { PICARD_SORT as PS_VAR  } from "../modules/picard_sort"
include { GET_CONSENSUS_COVERAGE; GET_VARIANT_READ_DEPTH  } from "../modules/get_coverage"
include { MERGE_MPILEUP_CONSENSUS } from "../modules/merge_mpileup_consensus"

workflow CONSENSUS_GEN {
    take:
        samples // tuple ([sample:s,replicate_id:r],[fastq1,fastq2])

    main:
       
        reference = file(params.reference_fasta)

        // each fastqc channel tuple contains a meta, replicate, [path(read1), path(read2)]
        // FASTQC(samples)

        bam = samples.map { meta, reads -> [meta, reads, reference] } 
          | BWA_MEM  // (meta, path(*sam))
          | FSI_CON // tuple val(meta), path("*.sorted.bam"), path("*.bai")

        if(params.tiled_amplicons){
          primer_bed = file(params.primer_bed)
          bam = bam.map { meta, sortedBam, bamIndex -> [ meta, sortedBam, bamIndex, primer_bed, reference ]} 
            | filter { _meta, sortedBam, _bamIndex, _primer_bed, reference-> sortedBam.size() >= 1000 } //filter out empty BAMs
            | AMPLICON_CLIP  //  tuple val(meta), path("*.primertrim.bam") // TODO replace is samtools amplicon trim
            | PS_CON // tuple val(meta), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")
        }


        consensus_sequence = bam.map{meta,bam,bai -> tuple( meta.sample, bam, bai)}.groupTuple() // [sample, [bams], [bais]]
          | MERGE_MPILEUP_CONSENSUS // sample, consensus

        // filter out empty consensus sequences
        consensus_sequence = consensus_sequence.filter { _sample, consensus -> consensus.size() >= 1000 }
        | GET_CONSENSUS_COVERAGE
        | filter { sample, consensus, coverage -> Float.parseFloat(coverage)>= params.consensus_coverage_cutoff}
        | map {sample ,consensus, _coverage -> [sample, consensus]}

        variant_bam =  bam.map{meta, bam, bai -> tuple(meta.sample, meta, bam)}
          .combine(consensus_sequence, by:0) // sample, meta, bam, consensus
          .map{_sample,meta,bam,consensus -> [meta, bam, consensus]}
          | BWA_REMAP // meta, sams
          | FSI_VAR // meta, sorted bam, index
        

        GET_VARIANT_READ_DEPTH(variant_bam)

        // BWA_MEM() doesn't output the consensus, so we rejoin it
        variant_bam_consensus = variant_bam.map{meta, bam, bai -> tuple(meta.sample, meta, bam, bai)}
          .combine(consensus_sequence, by:0) // sample, meta, bam, bai, consensus
          .map{_sample,meta,bam,bai,consensus -> tuple( meta, bam, bai, consensus)}  // meta, sorted bam, index, consensus


    emit:
        reads_and_consensus =  variant_bam_consensus // (meta, sorted bam, bam index, consensus.fa)
}