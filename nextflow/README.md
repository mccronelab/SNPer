# SNPer v0.1.0-Alpha

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
- Add process for FASTQC quality control
- Add process for BWA mem reference alignment
- Add workflow that manages consensus sequence generation
- Add main.nf


### v0.1.0

- Add README.md
