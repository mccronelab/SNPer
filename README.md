# SNPer v0.1.0.2-Beta

## Description

An in-development, Nextflow-managed viral in-host variant calling workflow. SNPer uses BWA-mem and iVar to align FASTQ files to a reference genome, construct a consensus sequence, and then call variants relative to the consensus. SNPer can handle either single or replicated samples and includes a primer mismatch detection step to improve variant calling in tiled amplicon data. SNPer also comes with a Docker image that manages workflow dependencies, which is [available here](https://quay.io/repository/mccronelab/snper). The intitial version of SNPer is heavily based on the Lauring Lab's VOC transmission pipeline as a jumping-off point: https://github.com/lauringlab/SARS-CoV-2_VOC_transmission_bottleneck

## Quick-Start Guide

### Docker or Apptainer (Recommended)
Required Software:
- Java, any version from 17-23 (to run Nextflow, [available here](https://www.oracle.com/java/technologies/downloads/?er=221886))
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
- Output: A tuple containing Sample ID, Replicate ID (if no replicates, this is the sample ID), FASTQ_1, FASTQ_2

### Build Consensus Sequence
- Input: A tuple containing Sample ID, Replicate ID, and FASTQ reads. Reference Sequence, Primer BEDfile
- Output: BAM files where FASTQ reads are aligned to a consensus genome.

1. Run FASTQC on reads, saving output.
2. With `BWA mem`, align reads to the reference sequence. Then, filter out unmapped reads and sort the output BAM file.
3. (Tiled Amplicon Only) Using `iVar trim`, mask primers on reads based on the contents of the primer BEDfile, then sort.
4. Group reads based on sample ID. Merge reads, including replicates of the same sample, and call a consensus sequence with `iVar consensus`.
5. Map FASTQ reads to the consensus genome, which will enable variant calling later on. Filter out unmapped reads and sort output BAM file.
6. Get coverage information for reads mapped to the consensus genome.

### Trimming Reads and Masking Primers (Tiled Amplicon Only)
- Input: BAM files where FASTQ reads are aligned to a consensus genome. Consensus Sequence.
- Output: Trimmed FASTQ reads aligned to a consensus genome. Reads with primers that do not match the consensus genome are removed.

1. Filter consensus genome FASTA files based on size. FASTA files smaller than 1kb are presumed to be files where a consensus was not successfully generated and are removed. This prevents further processing of associated reads, avoiding crashes that will occur later.
2. Filter out empty consensus-aligned BAM files, align primers to consensus genome, trim primers with `ivar trim`, and sort trimmed BAMs.
3. Reference primers are aligned to the consensus genome.
4. Using `iVar variants`, primer variants are called. The resulting BAM file is coverted to a BEDfile with `bedtools bamtobed`. Empty BEDfiles (where no primer variants were identified) are removed.
5. Using `iVar getmasked`, generate a list of primers with mismatches (variants) to the consensus genome.
6. Using `ivar removereads`, throw away reads with primer mismatches relative to the consensus sequence. Then, sort and index the filtered BAM file with `samtools`.

### Calling Variants with iVar
- Input: BAM files paired with BAM indices and their consensus sequence. Reference sequence GFF file. Reference sequence in FASTA format.
- Output: TSV files containing variants. Two types of TSVs are produced: one with variant positions relative to the consensus genome, and one with variant positions aligned to the reference genome.

1. Filter consensus genome FASTA files based on size. FASTA files smaller than 1kb are presumed to be files where a consensus was not successfully generated and are removed. This prevents further processing of associated reads, avoiding crashes that will occur later.
2. Submit the reference FASTA, reference GFF, and consensus FASTA to LiftOff. LiftOff aligns the reference and consensus sequences, then transfers GFF annotations where appropriate. This gives us a consensus GFF file.
3. Call variants relative to the consensus genome with `iVar variants`.
4. Call `convert_tsv_coords.py`, a relatively simple Python script that aligns the reference and consensus genome with MAFFT, then creates a mapping between positions on each genome based on the alignment. Used to convert consensus variant positions to their equivalent positions on the reference genome.

## Testing / dev

```
nextflow run ./ -profile test
```

_Requires docker_

## Parameters (last update: v0.3.0-Alpha)

-   `sample_sheet`: Path to CSV format sample sheet. The sample sheet has 4 fields: sample ID, replicate ID, and 2 paired-end read FASTQ files. Sample ID is used to relate data from separate replicates of the same sample.
- `reference_fasta`: A path to the reference genome for the replicon of interest.
- `reference_gff`: Path to GFF file describing ORFs on reference genome.
- `primer_bedfile`: Path to .bed file containig ARCTIC primers.
- `primer_fasta`: Path to file with FASTA sequences for each primer.
- `primer_pairs`: Path to TSV file that lists each pair of left and right primers.
- `output_dir`: Path where output will be stored.
- `consensus_min_qual_score`: Minimum score for base to be counted in consensus sequence generation. Default to 0, which somehow relates to indels.
- `consensus_threshold`: Minimum frequency threshold to call consensus (0-1, default 0).
- `consensus_min_depth`: Minimum depth to call consensus. `iVar consensus` recommends a default value of 10.
- `variant_minQ`: Minimum score for base to be counted in variant calling. Default to 30.
- `variant_min_mapQ`: Minimum quality score to be used in `samtools mpileup` during variant calling. Defaults to 20.
- `variant_freq_threshold`: Minimum variant frequency to pass `ivar variants`. Defaults to 0.02.
- `tiled_amplicons`: Boolean variable that indicates whether sequencing data comes from tiled amplicons,
    which requires additional filtering for primers.
