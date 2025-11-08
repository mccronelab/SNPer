// This workflow uses BWA mem, samtools, picard, and iVar to build a consensus genome

include { BWA_MEM; BWA_REMAP as BWA_REMAP_CON; BWA_REMAP as BWA_REMAP_POL } from "../modules/bwa_mem"
include { FILTER_SORT_INDEX as FSI_CON; FILTER_SORT_INDEX as FSI_VAR;
  FILTER_SORT_INDEX as FSI_POL  } from "../modules/filter_sort_index"
include { AMPLICON_CLIP } from "../modules/amplicon_clip"
include { DEDUPLICATE_READS as DEDUP_MIPS;  DEDUPLICATE_READS as DEDUP_HC } from "../modules/deduplicate_reads"
include { PICARD_SORT as PS_CON; PICARD_SORT as PS_VAR  } from "../modules/picard_sort"
include { GET_CONSENSUS_COVERAGE; GET_VARIANT_READ_DEPTH  } from "../modules/get_coverage"
include { MERGE_MPILEUP_CONSENSUS as ROUGH_CONSENSUS; 
  MERGE_MPILEUP_CONSENSUS as REFINE_CONSENSUS } from "../modules/merge_mpileup_consensus"

workflow CONSENSUS_GEN {
  take:
      samples // tuple ([sample:s,replicate_id:r],[fastq1,fastq2])

  main:
    reference = file(params.reference_fasta)

    // split into bams that require trimming and pre-trimmed bams
    consensus_bam = samples.map { meta, reads -> [meta, reads, reference] } 
      | BWA_MEM  // (meta, path(*sam))
      | FSI_CON // tuple val(meta), path("*.sorted.bam"), path("*.bai")
      | branch { meta, _sortedBam, _bamIndex ->
        preprocessed: meta.primer_id.toLowerCase() == "none"
        unprocessed:  true
      }

    // split preprocessed data into each sequencing type, so we can add these back to the workflow at the appropriate time
    // discard indexes- we'll get these back later
    preprocessed_bam = consensus_bam.preprocessed.map { meta, sortedBam, _bamIndex -> [meta, sortedBam]}
      | branch { meta, _sortedBam ->
        amplicon:       meta.sequencing_tech.toLowerCase() == "amplicon"
        mips:           meta.sequencing_tech.toLowerCase() == "mips"
        hybrid_capture: meta.sequencing_tech.toLowerCase() == "hybrid-capture"
        misc: true
        }

    // further split into different types of bam for differential processing. Add reference for amplicon_clip.nf
    unprocessed_bam = consensus_bam.unprocessed.map { meta, sortedBam, bamIndex -> [ meta, sortedBam, bamIndex, reference ]} 
      | filter { _meta, sortedBam, _bamIndex, _reference-> sortedBam.size() >= 1000 } //filter out empty BAMs
      | branch { meta, _sortedBam, _bamIndex, _reference ->
        amplicon:       meta.sequencing_tech.toLowerCase() == "amplicon"
        mips:           meta.sequencing_tech.toLowerCase() == "mips"
        hybrid_capture: meta.sequencing_tech.toLowerCase() == "hybrid-capture"
        misc: true
      }

    amplicon_bam_unp = unprocessed_bam.amplicon
    mips_bam_unp = unprocessed_bam.mips
    hc_bam_unp = unprocessed_bam.hybrid_capture
    misc_bam = unprocessed_bam.misc.concat(preprocessed_bam.misc) // these will throw an error later, so no need to keep them separate

    // if anything is in misc_bam, it's not being handled and sits around. let's throw an error if this is the case:
    misc_bam.take(1)
      | map { meta, _sortedBam, _bamIndex, _reference -> error "Unexpected sequencing tech type: $meta.sequencing_tech" }

    amplicon_bam = amplicon_bam_unp
      | AMPLICON_CLIP  //  tuple val(meta), path("*.primertrim.bam")
      | concat(preprocessed_bam.amplicon) // add back in BAMs derived from pre-trimmed FASTQs
      | PS_CON // tuple val(meta), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")

    mips_bam = mips_bam_unp // we'll add MIPs trimming process later, so preserve this structure for now
      | concat(preprocessed_bam.mips)
      | DEDUP_MIPS

    hc_bam = hc_bam_unp // we'll add hybrid-capture trimming later, so preserve this structure for now
      | concat(preprocessed_bam.hybrid_capture)
      | DEDUP_HC

    processed_bam = amplicon_bam.concat(mips_bam, hc_bam)

    consensus_sequence = processed_bam.map{ meta, bam, bai -> [meta, bam, bai, reference] }
      | ROUGH_CONSENSUS // sample, consensus
      // filter out consensus with no sequence, which sometimes occurs
      | filter { _sample, consensus_fa -> consensus_fa
        if (!consensus_fa.exists() || consensus_fa.size() == 0) {
          return false
          }
        def seq = consensus_fa.splitFasta(record: [seqString: true])
        if (seq) {
          return true
          }
        }

    variant_bam = processed_bam.map{ meta, bam, _bai -> tuple(meta.sample, meta, bam) }
      .combine(consensus_sequence, by:0) // sample, meta, bam, consensus
      .map{ _sample, meta, bam, consensus -> [meta, bam, consensus] }
      | BWA_REMAP_CON // meta, sams
      | FSI_VAR // meta, sorted bam, index

    polished_consensus = variant_bam
      | map {meta, bam, bai -> [meta.sample, meta, bam, bai] }
      | combine(consensus_sequence, by:0)
      | map {_sample, meta, bam, bai, consensus -> [meta, bam, bai, consensus] }
      | REFINE_CONSENSUS
      | GET_CONSENSUS_COVERAGE
      | filter { _sample, _consensus, coverage -> Float.parseFloat(coverage)>= params.consensus_coverage_cutoff }
      | map { sample, consensus, _coverage -> [sample, consensus] }

    variant_bam_polished = variant_bam.map{ meta, bam, _bai -> tuple(meta.sample, meta, bam) }
      | combine(polished_consensus, by:0)
      | map { _sample, meta, bam, polished_con -> [meta, bam, polished_con] }
      | BWA_REMAP_POL
      | FSI_POL

    GET_VARIANT_READ_DEPTH(variant_bam_polished)

    // BWA_MEM() doesn't output the consensus, so we rejoin it
    variant_bam_polished_consensus = variant_bam_polished.map{ meta, bam, bai -> tuple(meta.sample, meta, bam, bai) }
      .combine(polished_consensus, by:0) // sample, meta, bam, bai, consensus
      .map{ _sample, meta, bam, bai, consensus -> tuple( meta, bam, bai, consensus) }  // meta, sorted bam, index, consensus

  emit:
    reads_and_consensus =  variant_bam_polished_consensus // (meta, sorted bam, bam index, consensus.fa)
}