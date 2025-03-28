# SNPer v0.4.0.2-Alpha
Standard Nucleotide Pipeline (emerging resource)

## Description

An in-development, Nextflow-managed viral variant calling workflow. SNPer uses BWA-mem and iVar to align FASTQ files to a reference genome, construct a consensus sequence, and then call variants relative to the consensus. SNPer can either single or replicated samples and includes primer mismatch detection to improve variant calling in tiled amplicon data. SNPer also comes with a Docker image that manages most of its dependencies, which is [available here](https://quay.io/repository/mccronelab/snper). The intitial version of SNPer is heavily based on the Lauring Lab's VOC transmission pipeline as a jumping-off point: https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

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