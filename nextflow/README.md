# SNPer v0.1.1-Alpha

## Description

This is the initial version of SNPer. Heavily based on the Lauring Lab's VOC transmission pipeline as a jumping-off point: https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

## Parameters

- reference: A path to the reference genome for the replicon of interest.
- fastq_dir: Path to directory containing ALL input FASTQ files. Accepts .gz zipped files.
- primer_bedfile: Path to .bed file containig ARCTIC primers.
- output_dir: Path where output will be stored.
- min_qual_score: Minimum score for base to be counted. Default to 0, which somehow relates to indels.
- consensus_threshold: Minimum frequency threshold to call consensus (0-1, default 0).
- min_depth: Minimum depth to call consensus. `iVar consensus` recommends a default value of 10.

## Changelog

### v0.1.1-Alpha

- Remove unncessary equal signs from all process directives
- Add publishDir directive to process that generates consensus
- Correct error in changelog item placement
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