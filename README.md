# SNPer v2.4.0-Beta

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
nextflow run main.nf -profile test
```

Alternatively, if you have Apptainer, ensure that it is loaded/running:
```
nextflow run main.nf -profile test_apptainer
```

Alternatively, you can let Nextflow manage the repository with:
```
nextflow run https://github.com/mccronelab/SNPer.git -profile test
```

A list of parameters and their defaults is available with:
```
nextflow run main.nf --help
// for more information on a specific parameter, run
nextflow run main.nf --help_param <param>
```

## Workflow

### Sample Sheet Processing
- Input: Sample Sheet (CSV format)
- Output: A tuple containing a hash table with metadata values (meta), FASTQ_1, FASTQ_2. Run `python3 python/generate_sample_sheet.py -h` for more information on sample sheet fields.

### FASTQ Processing.
1. Submit input FASTQ files to `fastqc` for quality report generation (skippable via `--skip_fastqc`).
2. Submit FASTQ files to `fastp` for adapter removal and quality trimming, taking the paired output reads for further analysis. `fastp` also writes its own per-replicate JSON/HTML QC report to `fastp/`.

### Build Consensus Sequence
- Input: A tuple containing Sample ID, Replicate ID, and FASTQ reads. Reference Sequence, Primer BEDfile
- Output: BAM files where FASTQ reads are aligned to a consensus genome.

1. Align reads to the reference sequence with `bwa mem` and sort/index the output BAM file.
2a. (Tiled Amplicon reads only) Trim primers on aligned reads based on the contents of the primer BEDfile with `samtools ampliconclip`, then sort.
2b. (`umi` and `positional` reads only) Deduplicate aligned reads with `samtools markdup`.
3. Group reads based on sample ID. Merge reads, including replicates of the same sample, and call a consensus sequence with `iVar consensus`. Samples that failed to generate a consensus sequence are filtered out of the workflow.
4. Remap aligned reads to the consensus genome and sort/index output BAM file.
5. Call a 'polished' consensus sequence from remapped aligned reads. This is published to `consensus_seqs/`, named for its sample and segment. The first-pass 'rough' consensus is not published.
6. Remap aligned reads to the polished consensus sequence.
5. Get coverage information based on reads aligned to the consensus genome, filtering out any samples that fail to pass a user-defined coverage threshold (we typically use a 75% coverage threshold).

### Calling Variants with iVar
- Input: BAM files paired with BAM indices and their consensus sequence. Reference sequence GFF file. Reference sequence in FASTA format.
- Output: TSV files containing variants. Two types of TSVs are produced: one with variant positions relative to the consensus genome, and one with variant positions aligned to the reference genome.

The reference-coordinate TSV carries iVar's columns plus two added by SNPer:

| Column | Meaning |
| --- | --- |
| `REF_POS` | The variant's position relative to the reference genome, determined by a Nextclade alignment against the reference sequence. |
| `REF_POS_INSERTED` | `TRUE` when the variant sits on a base inserted relative to the reference. Such a base has no reference position of its own, so its `REF_POS` is inherited from the base preceding the insertion and is shared with it; this column is what tells the two apart. |

1. Filter consensus genome FASTA files based on size.
2. Submit the reference FASTA, reference GFF, and every segment consensus FASTA to Nextclade. Nextclade aligns each consensus against the reference and transfers the GFF annotations onto it, producing one combined annotation GFF plus the alignment. `split_nextclade_gff.py` splits that annotation into a consensus GFF file per sample.
3. Call variants relative to the consensus genome with `iVar variants`.
4. Call `convert_tsv_coords.py`, a Python script that reads Nextclade's alignment to create a mapping between positions on the consensus and the reference. Used to convert consensus variant positions to their equivalent positions on the reference genome.

## Testing / dev

Requires NF-Test, [installation instructions here](https://www.nf-test.com/installation/). For local execution in the Docker environment, use `--profile docker`.

```
nf-test test tests/main.nf.test --profile [profile] --ci
```

_Requires docker_

## Parameters 

Parameters are validated against `nextflow_schema.json` before any task is submitted. An unrecognised key — a typo in a `-params-file`, or a parameter removed in an earlier release — stops the run rather than being silently ignored, as do missing required paths and out-of-range values. `sample_sheet`, `reference_fasta`, `reference_gff` and `output_dir` are required; everything else has a default.

`nextflow run main.nf --help` prints every parameter with its type, default and description; `--help_param <name>` prints the full entry for a single parameter.

-   `sample_sheet`: Path to CSV format sample sheet. `sample` and `fastq1` are required, as is `fastq2` unless the row is interleaved. `sample` is used to relate data from separate replicates of the same sample. The optional columns `replicate_id`, `primer_id`, `interleaved` and `read_deduplication` override the corresponding params per row; `replicate_id` defaults to the sample ID. Relative FASTQ paths resolve against the sheet's own directory, not the launch directory.
- `reference_fasta`: A path to the reference genome for the replicon of interest. May hold one record per segment. Its record IDs are the segment names, and must match the `reference_gff` seqids and the primer BED chrom values exactly — the aligned BAM is split by reference name, so a mismatch affects annotation, primer trimming and the coordinate lift.
- `reference_gff`: Path to GFF file describing ORFs on reference genome, transferred onto each consensus by `nextclade run`. Seqids must match the `reference_fasta` record IDs.
- `primer_csv`: Path to CSV mapping `primer_id` values to primer BED files, with columns `primer_bedfile` and `primer_id`. Each sample's `primer_id` in the sample sheet selects its BED file from this table, so different samples in one run may use different primer schemes. The `primer_id` must match a row exactly, including case; `None` selects an empty BED. Primer BED chrom values must match the reference FASTA headers.
- `output_dir`: Path where output will be stored.
- `consensus_min_qual_score`: Minimum score for base to be counted in consensus sequence generation. Default to 0, which somehow relates to indels.
- `consensus_threshold`: Minimum frequency threshold to call consensus (0-1, default 0).
- `consensus_min_depth`: Minimum depth to call consensus. `iVar consensus` recommends a default value of 10.
- `consensus_coverage_cutoff`: Minimum fraction of the consensus genome which must be non-N characters (0-1, default 0.75). Applied per segment, so one poorly covered segment drops only itself.
- `consensus_min_baseQ`: Minimum base quality for a base to enter the `samtools mpileup` used for consensus generation. Defaults to 30. Distinct from `consensus_min_qual_score`, which is passed to `ivar consensus` itself.
- `consensus_min_mapQ`: Minimum mapping quality to be used in `samtools mpileup` during consensus generation. Defaults to 20.
- `variant_minQ`: Minimum score for base to be counted in variant calling. Default to 30.
- `variant_min_baseQ`: Minimum base quality for a base to pass the `samtools mpileup` used for variant calling and replicate coverage. Defaults to 30. Distinct from `variant_minQ`, which is passed to `ivar variants` itself.
- `variant_min_mapQ`: Minimum mapping quality for a read to pass the `samtools mpileup` used for variant calling and replicate coverage. Defaults to 20.
- `variant_freq_threshold`: Minimum variant frequency to pass `ivar variants`. Defaults to 0.01.
- `variant_min_depth`: Minimum depth for a position to report variants (default 10).
- `remove_unclipped_reads`: Boolean flag that controls whether `samtools ampliconclip` discards reads that are not trimmed (default true).
- `skip_qc`: Boolean flag that skips the whole read-QC subworkflow — both FastQC and `fastp`. Note this also skips adapter and quality trimming, not just report generation.
- `skip_fastqc`: Boolean flag that skips only FastQC, keeping `fastp` trimming. Useful for large datasets, since FastQC generates large output files, and `fastp` already reports the same read-level statistics (default false).
- `fastp_min_length`: Reads shorter than this after trimming are discarded (`fastp --length_required`, default 36).
- `interleaved`: Boolean flag that sets the default value for sample sheet processing. Use this flag if you've omitted this field from the input sample sheeet (default false).
- `primer_id_default`: Default setting for sample sheet processing, useful if your sample sheet doesn't have this field (and all samples use the same primers). Takes either a string or None.
- `read_deduplication`: Default value for sample sheet processing, useful if the input sample sheet doesn't include this field (default 'positional'). Selects how duplicate reads are handled: 'amplicon' (no deduplication; primers are clipped instead), 'umi' (UMI-aware `samtools markdup --barcode-name`, for MIPs and other UMI-tagged libraries), or 'positional' (coordinate-based `samtools markdup`, for hybrid-capture and other untagged libraries). Renamed from `sequencing_technique` in v2.3.0-Beta, along with its values.
