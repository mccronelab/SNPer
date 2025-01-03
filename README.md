# SNPer v0.1.1-Alpha
Standard Nucleotide Pipeline (emerging resource)

## Description

This is the initial version of SNPer. Heavily based on the Lauring Lab's VOC transmission pipeline as a jumping-off point: https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

## Parameters (last update: v0.2.0-Alpha)

- reference_fasta: A path to the reference genome for the replicon of interest.
- reference_gff: Path to GFF file describing ORFs on reference genome.
- fastq_dir: Path to directory containing ALL input FASTQ files. Accepts .gz zipped files.
- primer_bedfile: Path to .bed file containig ARCTIC primers.
- primer_pair_tsv: Path to .tsv file where each line has two columns: the left primer and right primer.
- primer_info_tsv: Path to a different .tsv file where each line has two columns: left primer and right primer.
- primer_fasta: Path to FASTA file with sequence of each primer.
- output_dir: Path where output will be stored.
- consensus_min_qual_score: Minimum score for base to be counted in consensus sequence generation. Default to 0, which somehow relates to indels.
- consensus_threshold: Minimum frequency threshold to call consensus (0-1, default 0).
- consensus_min_depth: Minimum depth to call consensus. `iVar consensus` recommends a default value of 10.
- variant_min_qual_score: Minimum score for base to be counted in variant calling. Default to 30.
- variant_min_mapQ: Minimum quality score to be used in `samtools mpileup` during variant calling. Defaults to 20.
- variant_freq_threshold: Minimum variant frequency to pass `ivar variants`. Defaults to 0.02.

## Development TODOs
- Investigate whether using consensus genomes in variant calling process is how this should work (trying the reference Wuhan01 genome, as in the Snakemake workflow, results in errors at `samtools mpileup`)
- Investigate whether the Wuhan01 reference genome GFF file used in the last process of the `call_variants_ivar` workflow (the specific process being `call_masked_variants`) does what we hope it does (accruately describes the ORF boundaries in consensus genomes)

## Changelog

### v0.2.0-Alpha (WIP)
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