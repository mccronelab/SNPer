// This workflow uses BWA mem, samtools, picard, and iVar to build a consensus genome

include { AMPLICON_CLIP } from "../modules/amplicon_clip"
include { BWA_MEM as BWA_MEM_REF; BWA_MEM as BWA_MEM_ROUGH; BWA_MEM as BWA_MEM_POL;
  BWA_REMAP as BWA_REMAP_ROUGH; BWA_REMAP as BWA_REMAP_POL } from "../modules/bwa_mem"
include { COUNT_MAPPED_READS } from "../modules/count_mapped_reads"
include { DEDUPLICATE_READS as DEDUP_UMI;  DEDUPLICATE_READS as DEDUP_POS } from "../modules/deduplicate_reads"
include { FILTER_SORT_INDEX as FSI_CON; FILTER_SORT_INDEX as FSI_VAR;
  FILTER_SORT_INDEX as FSI_POL  } from "../modules/filter_sort_index"
include { GET_CONSENSUS_COVERAGE; GET_VARIANT_READ_DEPTH  } from "../modules/get_coverage"
include { MERGE_MPILEUP_CONSENSUS as ROUGH_CONSENSUS; 
  MERGE_MPILEUP_CONSENSUS as POLISH_CONSENSUS } from "../modules/merge_mpileup_consensus"
include { PICARD_SORT as PS_CON; PICARD_SORT as PS_VAR  } from "../modules/picard_sort"
include { SORT_INDEX_BAM as SI_REF; SORT_INDEX_BAM as SI_ROUGH; SORT_INDEX_BAM as SI_POL } \
  from "../modules/sort_index_bam"
include { SPLIT_BAM_BY_SEGMENT } from "../modules/split_bam_by_segment"

// Name outputs with a segment token only for a multi-record reference; N=1 drops it
// (123.fa rather than 123_MN908947.3.fa).
def segmentLabel(segment, n_segments) {
    n_segments > 1 ? "_${segment}" : ""
}

workflow CONSENSUS_GEN {
  take:
      samples // tuple ([sample:s,replicate_id:r],[fastq1,fastq2])

  main:
    reference = file(params.reference_fasta)
    n_segments = reference.countFasta()

    // split into bams that require trimming and pre-trimmed bams
    consensus_bam = samples.map { meta, reads -> [meta, reads, reference] } 
      | BWA_MEM_REF  // (meta, path(*sam))
      | SI_REF // tuple val(meta), path("*.sorted.bam"), path("*.bai")
      | branch { meta, _sortedBam, _bamIndex ->
        preprocessed: meta.primer_id.toLowerCase() == "none"
        unprocessed:  true
      }

    // split preprocessed data by deduplication method, so we can add these back to the workflow at the appropriate time
    // discard indexes- we'll get these back later
    preprocessed_bam = consensus_bam.preprocessed.map { meta, sortedBam, _bamIndex -> [meta, sortedBam]}
      | branch { meta, _sortedBam ->
        amplicon:       meta.read_deduplication.toLowerCase() == "amplicon"
        umi:            meta.read_deduplication.toLowerCase() == "umi"
        positional:     meta.read_deduplication.toLowerCase() == "positional"
        misc: true
        }

    // further split into different types of bam for differential processing. Add reference for amplicon_clip.nf
    unprocessed_bam = consensus_bam.unprocessed.map { meta, sortedBam, bamIndex -> [meta, sortedBam, bamIndex, reference] }
      | branch { meta, _sortedBam, _bamIndex, _reference ->
        amplicon:       meta.read_deduplication.toLowerCase() == "amplicon"
        umi:            meta.read_deduplication.toLowerCase() == "umi"
        positional:     meta.read_deduplication.toLowerCase() == "positional"
        misc: true
      }

    amplicon_bam_unp = unprocessed_bam.amplicon
    umi_bam_unp = unprocessed_bam.umi
    positional_bam_unp = unprocessed_bam.positional
    misc_bam = unprocessed_bam.misc.concat(preprocessed_bam.misc) // these will throw an error later, so no need to keep them separate

    // if anything is in misc_bam, it's not being handled and sits around. let's throw an error if this is the case:
    misc_bam.take(1)
      | map { meta, _sortedBam, _bamIndex, _reference -> error "Unexpected read_deduplication value: $meta.read_deduplication" }

    amplicon_bam = amplicon_bam_unp
      | map { meta, bam, index, ref -> [meta, bam, index, ref, meta.primer_bedfile] }
      | AMPLICON_CLIP  //  tuple val(meta), path("*.primertrim.bam")
      | mix(preprocessed_bam.amplicon) // add back in BAMs derived from pre-trimmed FASTQs
      | PS_CON // tuple val(meta), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")

    umi_bam = umi_bam_unp // we'll add UMI trimming process later, so preserve this structure for now
      | map { meta, bam, _index, _ref -> [meta, bam] } // DEDUP takes (meta, bam); the branch added an index and reference for AMPLICON_CLIP
      | mix(preprocessed_bam.umi)
      | DEDUP_UMI

    positional_bam = positional_bam_unp
      | map { meta, bam, _index, _ref -> [meta, bam] }
      | mix(preprocessed_bam.positional)
      | DEDUP_POS

    // put all samples back into same channel, then split each per-replicate BAM
    // into per-segment BAMs. Splitting here (before the consensus grouping and
    // the remap fork) makes single-segment a special case of the multi-segment
    // path: N=1 for a single-record reference. Stamp meta.segment from the
    // parent dir name of each split BAM, building a FRESH map per segment so
    // map aliasing can't overwrite it across the fan-out.
    trimmed_bam = amplicon_bam.mix(umi_bam, positional_bam)

    // Removes samples with less mapped reads than --min_mapped_reads. If your primer
    // scheme does not match your samples, ampliconclip will drop almost
    // all reads, which can cause this filter to fail samples with deep coverage
    gated_bam = trimmed_bam
      | COUNT_MAPPED_READS // [meta, mapped read count]
      | map { meta, mapped -> [meta, mapped.trim() as Integer] }
      | join(trimmed_bam) // [meta, mapped, bam, bai]
      | branch { _meta, mapped, _bam, _bai ->
        keep:    mapped > (params.min_mapped_reads as Integer)
        dropped: true
      }

    gated_bam.dropped.subscribe { meta, mapped, _bam, _bai ->
      log.warn "Dropping ${meta.replicate_id}: ${mapped} mapped reads, at or below min_mapped_reads (${params.min_mapped_reads})."
    }

    gated_bam.dropped
      .map { meta, mapped, _bam, _bai -> "${meta.sample}\t${meta.replicate_id}\t${mapped}\n" }
      .collectFile(name: 'dropped_samples.tsv', storeDir: params.output_dir, sort: true,
                   seed: "sample\treplicate_id\tmapped_reads\n")

    processed_bam = gated_bam.keep
      | map { meta, _mapped, bam, bai -> [meta, bam, bai] }
      | SPLIT_BAM_BY_SEGMENT // [meta, [segment bams], [segment bais]]
      | map { meta, bams, bais -> [meta, [bams].flatten(), [bais].flatten()] } // force lists when N=1
      | transpose() // [meta, segment bam, segment bai]
      | map { meta, bam, bai ->
        def seg = bam.parent.name
        [meta + [segment: seg, segment_label: segmentLabel(seg, n_segments)], bam, bai]
      }
      | map { meta, bam, bai -> [meta.sample, meta.segment, meta, bam, bai] }

    // Group replicates of the same sample together for consensus calling, but
    // key on (sample, segment) so replicates of DIFFERENT segments aren't merged
    // back into one consensus. Every join below keys on [0,1] for the same
    // reason. Single-segment is the N=1 special case: one segment per sample, so
    // (sample, segment) collapses to the old per-sample grouping.
    bams_grouped_by_sample = processed_bam
      | groupTuple(by:[0,1]) // sample, segment

    // The rough consensus is an internal remap target: suffixed to keep it distinct
    // from the polished sequence, and not published.
    consensus_sequence = bams_grouped_by_sample.map{ sample, segment, meta, bam, bai -> [sample, segment, segmentLabel(segment, n_segments), meta, bam, bai, reference, "_rough", false] }
      | ROUGH_CONSENSUS // sample, segment, consensus
      // filter out consensus with no sequence, which sometimes occurs
      | filter { _sample, _segment, consensus_fa ->
        consensus_fa.exists() &&
        // .trim() evaluates to false on a string that contains only whitespace,
        // so if there are any lines in the sequence part of the FASTA that are not
        // only whitespace, the file passes the filter
        consensus_fa.readLines().any { line -> !line.startsWith(">") && line.trim() }
      }

    variant_bam = processed_bam
      | combine(consensus_sequence, by:[0,1]) // sample, segment, meta, bam, bai, consensus
      | map { _sample, _segment, meta, bam, _bai, consensus -> [meta, bam, consensus] }
      | BWA_REMAP_ROUGH // meta, sams
      | SI_ROUGH // meta, sorted bam, index
      | map { meta, bam, bai -> [meta.sample, meta.segment, meta, bam, bai] }

    variant_bams_grouped_by_sample = variant_bam
      | groupTuple(by:[0,1])
      | combine(consensus_sequence, by:[0,1])

    polished_consensus = variant_bams_grouped_by_sample
      // No suffix: this is the published sequence, so it carries the plain
      // sample/segment name in both the filename and the FASTA header. That name
      // flows on to the nextclade query rows and the per-sample GFF3s.
      | map { sample, segment, meta, bam, bai, consensus -> [sample, segment, segmentLabel(segment, n_segments), meta, bam, bai, consensus, "", true] }
      | POLISH_CONSENSUS
      | GET_CONSENSUS_COVERAGE
      // Coerce the cutoff: a --consensus_coverage_cutoff given on the command line
      // arrives as a String, which Groovy refuses to compare against a Float.
      | filter { _sample, _segment, _consensus, coverage -> Float.parseFloat(coverage) >= (params.consensus_coverage_cutoff as Float) }
      | map { sample, segment, consensus, _coverage -> [sample, segment, consensus] }

    variant_bam_polished = variant_bam
      | combine(polished_consensus, by:[0,1]) // sample, segment, meta, bam, bai, consensus
      | map { _sample, _segment, meta, bam, _bai, consensus -> [meta, bam, consensus] }
      | BWA_REMAP_POL
      | SI_POL

    // Re-join the polished consensus, dropped by the remap, for depth reporting.
    variant_bam_polished
      | map { meta, bam, bai -> [meta.sample, meta.segment, meta, bam, bai] }
      | combine(polished_consensus, by:[0,1])
      | map { _sample, _segment, meta, bam, bai, consensus -> [meta, bam, bai, consensus] }
      | GET_VARIANT_READ_DEPTH

  emit:
    variant_bams = variant_bam_polished // [meta, bam, bai]
    consensus = polished_consensus // [sample, segment, consensus]
}