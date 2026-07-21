// This workflow uses BWA mem, samtools, picard, and iVar to build a consensus genome

include { AMPLICON_CLIP } from "../modules/amplicon_clip"
include { BWA_MEM as BWA_MEM_REF; BWA_MEM as BWA_MEM_ROUGH; BWA_MEM as BWA_MEM_POL;
  BWA_REMAP as BWA_REMAP_ROUGH; BWA_REMAP as BWA_REMAP_POL } from "../modules/bwa_mem"
include { DEDUPLICATE_READS as DEDUP_MIPS;  DEDUPLICATE_READS as DEDUP_HC } from "../modules/deduplicate_reads"
include { FILTER_SORT_INDEX as FSI_CON; FILTER_SORT_INDEX as FSI_VAR;
  FILTER_SORT_INDEX as FSI_POL  } from "../modules/filter_sort_index"
include { GET_CONSENSUS_COVERAGE; GET_VARIANT_READ_DEPTH  } from "../modules/get_coverage"
include { MERGE_MPILEUP_CONSENSUS as ROUGH_CONSENSUS; 
  MERGE_MPILEUP_CONSENSUS as POLISH_CONSENSUS } from "../modules/merge_mpileup_consensus"
include { PICARD_SORT as PS_CON; PICARD_SORT as PS_VAR  } from "../modules/picard_sort"
include { SORT_INDEX_BAM as SI_REF; SORT_INDEX_BAM as SI_ROUGH; SORT_INDEX_BAM as SI_POL } \
  from "../modules/sort_index_bam"
include { SPLIT_BAM_BY_SEGMENT } from "../modules/split_bam_by_segment"

workflow CONSENSUS_GEN {
  take:
      samples // tuple ([sample:s,replicate_id:r],[fastq1,fastq2])

  main:
    reference = file(params.reference_fasta)

    // split into bams that require trimming and pre-trimmed bams
    consensus_bam = samples.map { meta, reads -> [meta, reads, reference] } 
      | BWA_MEM_REF  // (meta, path(*sam))
      | SI_REF // tuple val(meta), path("*.sorted.bam"), path("*.bai")
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
    unprocessed_bam = consensus_bam.unprocessed.map { meta, sortedBam, bamIndex -> [meta, sortedBam, bamIndex, reference] } 
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
      | map { meta, bam, index, ref -> [meta, bam, index, ref, meta.primer_bedfile] }
      | AMPLICON_CLIP  //  tuple val(meta), path("*.primertrim.bam")
      | concat(preprocessed_bam.amplicon) // add back in BAMs derived from pre-trimmed FASTQs
      | PS_CON // tuple val(meta), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")

    mips_bam = mips_bam_unp // we'll add MIPs trimming process later, so preserve this structure for now
      | concat(preprocessed_bam.mips)
      | DEDUP_MIPS

    hc_bam = hc_bam_unp // we'll add hybrid-capture trimming later, so preserve this structure for now
      | concat(preprocessed_bam.hybrid_capture)
      | DEDUP_HC

    // put all samples back into same channel, then split each per-replicate BAM
    // into per-segment BAMs. Splitting here (before the consensus grouping and
    // the remap fork) makes single-segment a special case of the multi-segment
    // path: N=1 for a single-record reference. Stamp meta.segment from the
    // parent dir name of each split BAM, building a FRESH map per segment so
    // map aliasing can't overwrite it across the fan-out.
    processed_bam = amplicon_bam.concat(mips_bam, hc_bam)
      | SPLIT_BAM_BY_SEGMENT // [meta, [segment bams], [segment bais]]
      | map { meta, bams, bais -> [meta, [bams].flatten(), [bais].flatten()] } // force lists when N=1
      | transpose() // [meta, segment bam, segment bai]
      | map { meta, bam, bai -> [meta + [segment: bam.parent.name], bam, bai] }
      | map { meta, bam, bai -> [meta.sample, meta, bam, bai] }
    
    // group different replicates of same sample together for consensus calling
    bams_grouped_by_sample = processed_bam
      | groupTuple(by:0)

    consensus_sequence = bams_grouped_by_sample.map{ sample, meta, bam, bai -> [sample, meta, bam, bai, reference, "_rough"] }
      | ROUGH_CONSENSUS // sample, consensus
      // filter out consensus with no sequence, which sometimes occurs
      | filter { sample, consensus_fa ->
        consensus_fa.exists() &&
        // .trim() evaluates to false on a string that contains only whitespace,
        // so if there are any lines in the sequence part of the FASTA that are not
        // only whitespace, the file passes the filter
        consensus_fa.readLines().any { line -> !line.startsWith(">") && line.trim() }
      }
    
    variant_bam = processed_bam
      | combine(consensus_sequence, by:0) // sample, meta, bam, bai, consensus
      | map { _sample, meta, bam, _bai, consensus -> [meta, bam, consensus] }
      | BWA_REMAP_ROUGH // meta, sams
      | SI_ROUGH // meta, sorted bam, index
      | map { meta, bam, bai -> [meta.sample, meta, bam, bai] }

    variant_bams_grouped_by_sample = variant_bam
      | groupTuple(by:0)
      | combine(consensus_sequence, by:0)

    polished_consensus = variant_bams_grouped_by_sample
      | map { sample, meta, bam, bai, consensus -> [sample, meta, bam, bai, consensus, "_polished"] }
      | POLISH_CONSENSUS
      | GET_CONSENSUS_COVERAGE
      | filter { _sample, _consensus, coverage -> Float.parseFloat(coverage)>= params.consensus_coverage_cutoff }
      | map { sample, consensus, _coverage -> [sample, consensus] }

    variant_bam_polished = variant_bam
      | combine(polished_consensus, by:0) // sample, meta, bam, bai, consensus
      | map { _sample, meta, bam, _bai, consensus -> [meta, bam, consensus] }
      | BWA_REMAP_POL
      | SI_POL

    GET_VARIANT_READ_DEPTH(variant_bam_polished)

  emit:
    variant_bams = variant_bam_polished // [meta, bam, bai]
    consensus = polished_consensus // [sample, consensus]
}