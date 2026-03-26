# SNPer v2.2.0-Beta

## Description

An in-development, Nextflow-managed viral in-host variant calling workflow. SNPer uses BWA-mem and iVar to align FASTQ files to a reference genome, construct a consensus sequence, and then call variants relative to the consensus. SNPer can handle either single or replicated samples and includes a primer mismatch detection step to improve variant calling in tiled amplicon data. SNPer also comes with a Docker image that manages workflow dependencies, which is [available here](https://quay.io/repository/mccronelab/snper). The intitial version of SNPer is heavily based on the Lauring Lab's VOC transmission pipeline as a jumping-off point: https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

## Quick-Start Guide

### Docker or Apptainer (Recommended)
Required Software:
- Java, any version from 17-21 (to run Nextflow, [available here](https://www.oracle.com/java/technologies/downloads/?er=221886))
- Nextflow (to run the workflow, [available here](https://www.nextflow.io/docs/latest/install.html#install-nextflow))
- Docker or Apptainer (to run the container. Docker is [available here](https://docs.docker.com/get-started/get-docker/), ask your local HPC staff about Apptainer)

To download the latest version of SNPer and run it locally, run:
```
git clone https://github.com/mccronelab/SNPer.git
cd SNPer/
```

To test SNPer, ensure that Docker is running:
```
nextflow run https://github.com/mccronelab/SNPer.git -profile test
```

Alternatively, if you have Apptainer, ensure that it is loaded/running:
```
nextflow run https://github.com/mccronelab/SNPer.git -profile test_apptainer
```

Alternative, you can run the local copy of the main workflow with:
```
nextflow run main.nf -profile test
```

## Workflow

### Sample Sheet Processing
- Input: Sample Sheet (CSV format)
- Output: A tuple containing a hash table with metadata values (meta), FASTQ_1, FASTQ_2. Run `python3 python/generate_sample_sheet.py -h` for more information on sample sheet fields.

### FASTQ Processing.
1. Submit input FASTQ files to `fastqc` for quality report generation.
2. Submit FASTQ files to `trimmomatic`, taking trimmed and paired output reads for further analysis.

### Build Consensus Sequence
- Input: A tuple containing Sample ID, Replicate ID, and FASTQ reads. Reference Sequence, Primer BEDfile
- Output: BAM files where FASTQ reads are aligned to a consensus genome.

1. Align reads to the reference sequence with `bwa mem` and sort/index the output BAM file.
2a. (Tiled Amplicon reads only) Trim primers on aligned reads based on the contents of the primer BEDfile with `samtools ampliconclip`, then sort.
2b. (MIPs reads only) Deduplicate aligned reads with `samtools markdup`.
3. Group reads based on sample ID. Merge reads, including replicates of the same sample, and call a consensus sequence with `iVar consensus`. Samples that failed to generate a consensus sequence are filtered out of the workflow.
4. Remap aligned reads to the consensus genome and sort/index output BAM file.
5. Call a 'polished' consensus sequence from remapped aligned reads.
6. Remap aligned reads to the polished consensus sequence.
5. Get coverage information based on reads aligned to the consensus genome, filtering out any samples that fail to pass a user-defined coverage threshold (we typically use a 75% coverage threshold).

### Calling Variants with iVar
- Input: BAM files paired with BAM indices and their consensus sequence. Reference sequence GFF file. Reference sequence in FASTA format.
- Output: TSV files containing variants. Two types of TSVs are produced: one with variant positions relative to the consensus genome, and one with variant positions aligned to the reference genome.

1. Filter consensus genome FASTA files based on size.
2. Submit the reference FASTA, reference GFF, and consensus FASTA to LiftOff. LiftOff aligns the reference and consensus sequences, then transfers GFF annotations where appropriate. This gives us a consensus GFF file.
3. Call variants relative to the consensus genome with `iVar variants`.
4. Call `convert_tsv_coords.py`, a relatively simple Python script that aligns the reference and consensus genome with MAFFT, then creates a mapping between positions on each genome based on the alignment. Used to convert consensus variant positions to their equivalent positions on the reference genome.

## Testing / dev

Requires NF-Test, [installation instructions here](https://www.nf-test.com/installation/). For local execution in the Docker environment, use `--profile docker`.

```
nf-test test tests/main.nf.test --profile [profile] --ci
```

_Requires docker_

## Parameters 

-   `sample_sheet`: Path to CSV format sample sheet. The sample sheet has 4 fields: sample ID, replicate ID, and 2 paired-end read FASTQ files. Sample ID is used to relate data from separate replicates of the same sample.
- `reference_fasta`: A path to the reference genome for the replicon of interest.
- `reference_gff`: Path to GFF file describing ORFs on reference genome.
- `primer_bed`: Path to .bed file containing short read primers.
- `output_dir`: Path where output will be stored.
- `consensus_min_qual_score`: Minimum score for base to be counted in consensus sequence generation. Default to 0, which somehow relates to indels.
- `consensus_threshold`: Minimum frequency threshold to call consensus (0-1, default 0).
- `consensus_min_depth`: Minimum depth to call consensus. `iVar consensus` recommends a default value of 10.
- `consensus_coverage_cutoff`: Minimum percentage of the consensus genome which must be non-N characters (0-1, default 0.75).
- `variant_minQ`: Minimum score for base to be counted in variant calling. Default to 30.
- `variant_min_mapQ`: Minimum quality score to be used in `samtools mpileup` during variant calling. Defaults to 20.
- `variant_freq_threshold`: Minimum variant frequency to pass `ivar variants`. Defaults to 0.02.
- `variant_min_depth`: Minimum depth for a position to report variants (default 10).
- `remove_unclipped_reads`: Boolean flag that controls whether `samtools ampliconclip` discards reads that are not trimmed (default true).
- `skip_qc`: Boolean flag that skips QC processes. Useful for test runs on large datasets, since FastQC generates large output files.
- `trimmomatic_jarfile`: Path to Trimmomatic's jarfile. You should only change this parameter if you're NOT running SNPer via its Docker container.
- `interleaved`: Boolean flag that sets the default value for sample sheet processing. Use this flag if you've omitted this field from the input sample sheeet (default false).
- `primer_id_default`: Default setting for sample sheet processing, useful if your sample sheet doesn't have this field (and all samples use the same primers). Takes either a string or None.
- `sequencing_technique`: Default value for sample sheet processing, useful if the input sample sheet doesn't include this field. SNPer supports 'amplicon' and 'mips' for this parameter.