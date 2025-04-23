# Changelog

## v0.1.0.2-Beta
- Set test profile to faster test with files on GitHub.
    - Add MN908947.3 GFF3 annotation file.
- Update README parameter list.
- Correct double calling of parameter in `ivar_primer_variants.nf`.

## v0.1.0.1-Beta
- Fix typos in README.

## 0.1.0-Beta
- Open up repository to general public.
- Address major primer mismatch detection bug:
    - `trim_and_mask.nf`: Realign primers to consensus genome before trimming primers.
- Tweak README to make some text clearer.

## v0.4.1-Alpha
- Begin reworking documentation in preparation for making the repository public:
    - Move changelog to its own Markdown file.
    - Update README parameter list.
    - Add Quick Start Guide to README.
    - Add Description to README.
    - Add step-by-step descriptions of each subworkflow.
- Replace all double quotes with single quotes in `nextflow.config` to resolve error.
-  Add dynamic resource allocation to resource-intensive processes, triggered when a crash occurs due to exceeding resource limitations:
    - `convert_tsv_coords.nf`
    - `fastqc.nf`
    - `ivar_variants.nf`
    - `merge_mpileup_consensus.nf`
    - `picard_sort.nf`
    - Add maxRetries 3 for all processes in `nextflow.config`.
- Fix spacing issues in `ivar_variants.nf`, `liftoff.nf` process definitions.

## v0.4.0.2-Alpha
- Remove `stageInMode` directive from `convert_tsv_coords.nf`, as it was causing an input bug.
- Rewrite `convert_tsv_coords.py` to combine separate reference, consensus files into one FASTA file.
- Increase default `queueSize`.

## v0.4.0.1-Alpha
- Add MAFFT to image.
- Fix bug related to convert_tsv_coords.py being called as an executable.

## v0.4.0-Alpha
- Update sample sheet to explicitly include a replicate ID field, so we don't have to try to 
extract it from file names during process execution.
    - Update `bwa_mem.nf` and `fastqc.nf` processes to handle new sample sheet configuration.
    - Update `build_consensus.nf` and `process_sample_sheet.nf` workflows to handle new sample
    sheet configuration.
- Add BED file size filter to `trim_and_mask.nf` workflow to filter out completely empty files.
- Start tracking supplemental Python script `generate_sample_sheet.py`. It's not a particularly
robust solution, but it's been sufficient so far.
- Add convert_tsv_coords.py in /bin/.
- Add convert_tsv_coords.nf.
- Add call to convert_tsv_coords.nf to call_variants_ivar.nf
- Ensure consistent use of spaces in files in workflows/.
- Add Python's Bio package to SNPer image.
    

## v0.3.1-Alpha
- Add missing params to test profile.
- Add indexing to align_fasta_filter_sort.nf.
- Change coverage output directory name to make it clear we have coverage for replicates now.
- Fix input order in mask_primers.nf.
- Fix bugs in build_consensus.nf:
    - Rename processes invoked more than once.
    - Filter out empty consensus genomes.
    - Update join() to combine() where we expect 1 consensus genome to match to multiple replicates.
- Fix bugs in call_variants_ivar.nf:
    - join() uses parentheses, not {}.
    - Drop unnecessary index files before calling ivar_variants.nf.
- Fix bugs in trim_and_mask.nf:
    - Fix slightly incorrect input parameter references, rename variables for better clarity.
    - Remove duplicated consensus sequences in split.consensus_seq during combine().
- Fix order of arguments in script block of liftoff.nf.
- Add Liftoff to Docker image.
- Fix extra space throwing off a comma in ivar_variants.nf.
- Rename GFF input parameter in ivar_variants (reference_gff to gff), since we now use GFFs realigned to each consensus.
- Update rhino_test params.
- Fix output name bug in ivar_primer_variants.nf.
- Fix input order mixup in remove_masked_sort_index.nf.

## v0.3.0-Alpha

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


## v0.2.2-Alpha

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


## v0.2.1-Alpha

- refactor workflow to take sample sheet as input
- one consensus for each sample - one variant tsv for each sequencing library
- move all files out of `nextflow` directory
- adds default parameters to the config file


## v0.2.0-Alpha
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


## v0.1.1-Alpha

- Remove unncessary equal signs from all process directives.
- Add publishDir directive to process that generates consensus.
- Correct error in changelog item placement.
- Change `primer_bedfile` in `main.nf` to be an instance of a file, rather than a Channel of Paths. As a result, it can be provided to more than a single process.

## v0.1.0-Alpha

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