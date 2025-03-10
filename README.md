# SNPer v0.3.1-Alpha
Standard Nucleotide Pipeline (emerging resource)

## Description

This is the initial version of SNPer. Heavily based on the Lauring Lab's VOC transmission pipeline as a jumping-off point: https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

## Testing / dev

```
nextflow run ./ -profile test
```

_Requires docker_

## Parameters (last update: v0.3.0-Alpha)

- reference_fasta: A path to the reference genome for the replicon of interest.
- reference_gff: Path to GFF file describing ORFs on reference genome.
- primer_bedfile: Path to .bed file containig ARCTIC primers.
- output_dir: Path where output will be stored.
- consensus_min_qual_score: Minimum score for base to be counted in consensus sequence generation. Default to 0, which somehow relates to indels.
- consensus_threshold: Minimum frequency threshold to call consensus (0-1, default 0).
- consensus_min_depth: Minimum depth to call consensus. `iVar consensus` recommends a default value of 10.
- variant_min_qual_score: Minimum score for base to be counted in variant calling. Default to 30.
- variant_min_mapQ: Minimum quality score to be used in `samtools mpileup` during variant calling. Defaults to 20.
- variant_freq_threshold: Minimum variant frequency to pass `ivar variants`. Defaults to 0.02.
- tiled_amplicons: Boolean variable that indicates whether sequencing data comes from tiled amplicons,
    which requires additional filtering for primers.

## Changelog

### v0.3.0-Alpha

- Add workflow for trimming primers, removing reads with primer mismatches (trim_and_mask).
    - Add process for aligning FASTA with `bwa mem`, which requires specific seed and threshold settings to work well (align_fasta_filter_sort).
    - Add process for calling variants (mismatches) on primer sequences (ivar_primer_variants).
- Update get_coverage process to use variant calling MapQ threshold.
- Add `--reference` to `samtools mpileup` in ivar_variants and enabled BAQ.
- Update parameter names in mask_primers to be more informative.
- Update parameter names in remove_masked_sort_index to be more informative.
- Refactor build_consensus workflow to emit reads aligned to the consensus, indexes, and the consensus itself. This allows us to potentially go straight to variant calling, for datasets that don't require primer trimming and read removal (tiled amplicon sequencing samples).
- Rework call_variants_ivar workflow to remove read trimming steps.
- Add ORF GFF remapping process (liftoff.nf), and add to call_variants_ivar workflow.
- Create samtools_sort process to replace picard_sort, allowing us to remove a dependency.


### v0.2.2-Alpha

- Add default resource allocations to each process based on usage information
    - bwa_mem.nf
    - fastqc.nf
    - filter_sort_index.nf
    - get_coverage.nf
    - ivar_variants.nf
    - merge_mpileup_consensus.nf
    - picard_sort.nf

- Update bwa_mem.nf to work with MIDGE dataset naming scheme
- Add SLURM executor profile to nextflow.config
- Address variable naming conflicts and syntax issues in call_variants_ivar.nf
- Add file size filtering to call_variants_ivar.nf to remove empty consensus genomes or BAM files from workflow
- Update README parameters list


### v0.2.1-Alpha

- refactor workflow to take sample sheet as input
- one consensus for each sample - one variant tsv for each sequencing library
- move all files out of `nextflow` directory
- adds default parameters to the config file

<<<<<<< Updated upstream
### v0.2.0-Alpha
=======

### v0.2.0-Alpha

>>>>>>> Stashed changes
- Add workflow that calls variants with iVar (call_variants_ivar).
    - Add process that maps to a reference, filters and sorts mapped reads (bwa_mem_filter_sort)
    - Add process that generates bwa index, samtools faidx (bwa_samtools_index)
    - Add process that converts BAM files to BED files (bam_to_bed)
    - Add process that calls variants in primer sequences (call_primer_variants)
    - Add process that detects mismatches between primer sequences and consensus sequences, and flags associated reads for masking (mask_primers)
    - Add process that removes masked reads, then sorts and indexes masked BAM files (remove_masked_sort_index)
    - Add process that calls variants against consensus genomes, using Wuhan01 ORFs (call_masked_variants)
- Update workflow parameters
    - Add parameter variant_min_qual_score
    - Add parameter variant_min_mapQ
    - Add parameter variant_min_depth
    - Add parameter variant_freq_threshold
    - Add parameter primer_info_tsv
    - Rename parameter min_qual_score to consensus_min_qual_score, to avoid confusion with a similar parameter used in variant calling.
    - Rename parameter min_depth to consensus_min_depth, to avoid confusion as above.
- Remove deprecated `-S` flag from `bwa mem` processes.
- Move README up one level out of `SNPer/nextflow/`.
- Modify `merge_mpileup_consensus` to emit a tuple with a key and FASTA path, instead of just the path.
- Add new workflow parameters: primer_fasta, primer_pair_tsv, reference_gff. Update reference parameter to reference_fasta
- Create Dockerfile to manage workflow dependencies. 
- Update `picard_sort.nf` to call `PicardCommandLine`.
- Update `merge_mpileup_consensus` to emit tuple with key for later joining, instead of just consensus sequences.
- Set workflow to strict mode. This has a number of effects, but of particular interest to us, it sets the workflow to fail if a `join()` operation is called on a channel with duplicate keys, and will fail if a key in one channel doesn't have a partner in the other. See [Nextflow Documentation](https://www.nextflow.io/docs/latest/reference/feature-flags.html) for the full list of effects.


### v0.1.1-Alpha

- Remove unncessary equal signs from all process directives.
- Add publishDir directive to process that generates consensus.
- Correct error in changelog item placement.
- Change `primer_bedfile` in `main.nf` to be an instance of a file, rather than a Channel of Paths. As a result, it can be provided to more than a single process.

### v0.1.0-Alpha

- Add README.md
- Add process that generates FASTQC quality reports.
- Add process that aligns reads to reference genome with BWA mem.
- Add process that filters, sorts, and indexes BAM files with samtools.
- Add process that uses iVar trim to remove primers from reads.
- Add process that sorts BAM files with Picard.
- Add process that gets coverage information with samtools.
- Adds process that runs samtools merge and mpileup, then generates consensus with iVar consensus.
- Add workflow that manages consensus sequence generation.
- Add main.nf