# Changelog

## v2.4.0-Beta
- Replace Trimmomatic with `fastp`. **Behaviour change: variant and coverage output shifts, because adapter clipping now runs for the first time** — Trimmomatic's `ILLUMINACLIP` adapter file was never staged, so the step was skipped while the task still exited 0. On the influenza test sample, 4,638 reads / 37,245 bases are now adapter-trimmed.
- Quality filtering follows fastp's defaults (`-q 15`, `-u 40`, `-n 5`); `MINLEN:36` becomes the `fastp_min_length` param. Cost on the influenza test sample: 4,322 of 2,404,726 reads dropped (0.18%). The variant and coverage TSVs change; consensus FASTAs and Liftoff GFF3s do not.
- Fix interleaved input, which was broken whenever `skip_qc` was false.
- Add `skip_fastqc`, which skips FastQC while keeping trimming. fastp's JSON/HTML reports are published to `fastp/`.
- Remove the `trimmomatic_jarfile` param. Params files that still set it are unaffected — the key is ignored.
- Add `tests/fastp.nf.test`, and repoint `testSampleSheet_interleaved.csv` at a committed fixture (it referenced three FASTQs that were never in the repo).
- Resolve sample-sheet FASTQ paths and primer-CSV BED paths relative to the CSV that names them; absolute paths and remote URIs are unaffected. Add a `manifest{}` block with `defaultBranch = 'main'`, required for `nextflow run <github-url>`. Test-profile `output_dir` now uses `${launchDir}`, and `apptainer.autoMounts` is set.
- Invoke `bin/calculate_coverage.py` and `bin/convert_tsv_coords.py` through Nextflow's bin-on-PATH mechanism rather than interpolating `${projectDir}`.
- Move the container base to Ubuntu 26.04 LTS (25.04 is EOL) and drop the arch-specific `JAVA_HOME`, which broke arm64 builds. No output change.
- Report `replicate_coverage/*.tsv` from the same `samtools mpileup` used for variant calling, replacing `samtools depth`. **Behaviour change: published depths fall**, because mpileup counts overlapping mates once where `depth` counted both. `consensus_coverage_cutoff` is unaffected — `calculate_coverage.py` reads the consensus FASTA, not this TSV.
- Fix `variant_min_baseQ` and `variant_min_mapQ` being applied to the wrong filters in `get_coverage.nf`: `samtools depth` reads `-q` as base quality and `-Q` as mapping quality, the reverse of `mpileup`. No effect when both params are set to the same value.
- Add `-a` to the `samtools mpileup` in `ivar_variants.nf` so it matches the coverage pileup. No output change.
- Add a `REF_POS_INSERTED` column to `variants/*.ref_coords.tsv`, `TRUE` where the variant sits on a base inserted relative to the reference and its `REF_POS` is therefore inherited from the preceding base rather than its own. **Output format change: every `.ref_coords.tsv` gains a trailing column.** Readers with a fixed column spec need updating; `REF_POS` itself is unchanged.
- Derive `REF_POS` from Nextclade's alignment rather than a separate MAFFT run per (sample, segment). `CONVERT_TSV_COORDS` no longer takes a reference or runs an aligner, so the coordinate lift and the published GFF3 annotations can no longer disagree about where a consensus sits on the reference. Published values are unchanged on both test sets. **Behaviour change: `REF_POS` moves for variants near an indel**, where MAFFT's placement drifted from the reference base it should have landed on. A variant inside an insertion still carries the preceding reference position, now with a warning naming it.
- Replace Liftoff with Nextclade for transferring reference annotations onto each consensus. One `nextclade run` per segment annotates every sample at once, replacing one Liftoff run per (sample, segment), and publishes each segment's alignment and QC table to a new `msa/` directory. **Behaviour change: the published GFF3s change, and one amino-acid position is corrected.** On a low-coverage SARS-CoV-2 consensus Liftoff placed ORF3a 12 nt late and truncated S to a length not divisible by 3; Nextclade places both correctly, moving `POS_AA` from 246 to 250 for one variant in the test set. Codons and amino-acid calls are unaffected. On segments with overlapping CDS features (M) iVar emits its per-feature rows in a different order, so those rows are reordered but not altered. Note this batches a segment across samples: every sample's consensus for a segment must finish before that segment is annotated, and a consensus Nextclade cannot align now fails the segment rather than only its own sample.
- Add `bin/split_nextclade_gff.py`, which splits Nextclade's combined annotation GFF3 into one file per consensus and normalizes the attribute column. Nextclade serializes attributes from a hash map with a per-run seed and rewrites feature IDs with the query's index in the batch, so without this every GFF3 would change md5 between identical runs and adding one sample would change every other sample's iVar `GFF_FEATURE` values. Covered by `tests/test_split_nextclade_gff.py`.
- Validate params against a new `nextflow_schema.json` (via the `nf-schema` plugin) before any task is submitted. **Behaviour change: an unrecognised param now stops the run.** Strict mode doesn't catch these — a misspelled key in a `-params-file` landed in the params map, nothing read it, and the run proceeded with the default. Params files must drop any deprecated parameter the pipeline no longer reads.
- Change the `read_deduplication` default from `amplicon` to `positional`. **Behaviour change for sample sheets with no `read_deduplication` column**, which switch from primer clipping to coordinate-based `samtools markdup`. Sheets and params files that set the value explicitly are unaffected, including all three nf-test fixtures.
- Remove `tiled_amplicons` and `primer_bed`, neither of which has been in use since v2.0.1-Beta. `read_deduplication` selects primer trimming per sample, and `primer_csv` supplies the primer BEDfile. Params files setting either will stop on schema validation.
- Drop MAFFT and Liftoff, along with Liftoff's `minimap2` and `libparasail-dev` dependencies, from the container image. No process invokes them any more. No output change.
- Rename `data/reference/pr8_multisegment_test.gff3` to `pr8_renamed_ha_corrected.gff3` and correct its HA annotation against EF467821.1. Its seqids are the segment names the reference FASTA uses rather than NCBI accessions. See file header for details. **Behaviour change: `ivar variants` no longer annotates HA positions 1731-1733 as a missense D→A — they are 3'UTR.** All eight published GFF3s change, due to output headers; only HA's annotation content differs. Params files pointing at the old path must be updated.

## v2.3.0-Beta
- Add support for multi-segment viruses (e.g. influenza) from a single, fixed multi-segment reference. Reads map to the full multi-record reference with one `bwa index`, then the aligned BAM is split by reference name so every downstream step runs per segment. Single-segment references are the N=1 special case of this path — no config flag, no separate codepath.
- Add `split_bam_by_segment.nf`, which splits each per-replicate BAM into one BAM per mapped reference name and stamps `meta.segment`. Wire it into `build_consensus.nf` after primer trimming / deduplication, before consensus grouping.
- Retain unplaced reads (`RNAME '*'`) in every per-segment BAM emitted by `split_bam_by_segment.nf`, so `BWA_REMAP` — which rebuilds FASTQ from that BAM — can attempt to place them against the sample's own consensus. They are inert in `samtools mpileup` and `samtools depth`, so consensus and coverage output are unchanged; each unplaced read is duplicated into all N segment BAMs (+2.7% records on the influenza test set, none at N=1).
- Re-key replicate grouping and consensus joins in `build_consensus.nf` and `call_variants_ivar.nf` on `(sample, segment)` (`groupTuple`/`combine` `by:[0,1]`), so replicates of different segments are no longer merged back together at consensus, remap, or variant-calling time.
- Thread the segment into per-segment output names to prevent collisions in shared publish directories: consensus FASTAs in `merge_mpileup_consensus.nf` (`${sample}_${segment}${suffix}`) and remapped BAMs in `bwa_mem.nf` (`${replicate}_${segment}.remap.bam`), which the variant and coverage outputs inherit. Multi-record references only — see the N=1 note below. Add the segment to `get_coverage.nf` (`GET_CONSENSUS_COVERAGE`) and `liftoff.nf` signatures so the join key survives across them.
- Name the variant TSVs in `ivar_variants.nf` and `convert_tsv_coords.nf` from `meta` (`${replicate}${segment_label}`) rather than from the input BAM/TSV filename. `simpleName` truncates at the first dot, so segment names containing dots were silently cut short in those outputs (`MN908947.3` published as `..._MN908947.ref_coords.tsv`) and two segments sharing a pre-dot prefix would have collided. Names are unchanged for dot-free segment names.
- Extract the per-segment reference record in `convert_tsv_coords.nf` before the coordinate lift-back, since `convert_tsv_coords.py` assumes a single-record reference and would otherwise align the wrong pair against a multi-segment reference.
- Suppress the segment token in output names when the reference holds a single record, so single-segment runs keep their pre-refactor filenames and FASTA headers (`123_polished.fa`, `123_A.remap.sorted_coverage.tsv`) despite going through the same split-by-segment path. The record count comes from `countFasta()` on `reference_fasta` in `build_consensus.nf`, which stamps `meta.segment_label` (`""` at N=1, `_${segment}` otherwise) at the split and passes it to `MERGE_MPILEUP_CONSENSUS`; `meta.segment` still carries the real reference name as the join key. Names derived downstream (Liftoff GFF3s, coverage TSVs, `.ref_coords.tsv`) inherit the change.
- Add a PR8 (influenza A, 8 segments) fixture set for exercising the multi-segment path: `data/reference/pHW2000_PR8-N_A.fa`, `data/reference/pr8_multisegment_test.gff3`, `testSampleSheet_multi_segment.csv`, and `multi_segment_test.yaml`.
- Note: in the reference GFF, sub-CDS features such as `mature_protein_region_of_CDS` and `signal_peptide_region_of_CDS` must be parented to their `gene`, not to a `CDS`. Liftoff walks the feature hierarchy expecting an intervening `mRNA` layer and exits with an `IndexError` when a domain feature's parent is a CDS. This applies to any user-supplied `reference_gff`, not just the influenza fixture.
- Add `consensus_min_mapQ`, replacing the use of `variant_min_mapQ` in `merge_mpileup_consensus.nf`. Consensus generation no longer changes when variant-calling parameters are retuned. Defaults to 20, matching the previous effective value; param files that pin `variant_min_mapQ` should pin `consensus_min_mapQ` alongside it to preserve behaviour.
- Add `consensus_min_baseQ` and `variant_min_baseQ`, replacing the hardcoded `samtools mpileup -Q 30` in `merge_mpileup_consensus.nf` and `ivar_variants.nf`, and the hardcoded `samtools depth -Q 30` in `get_coverage.nf`. Both default to 30, so behaviour is unchanged. These are the base-quality floors applied before iVar sees the pileup, and are distinct from `consensus_min_qual_score` and `variant_minQ`, which are passed to `ivar consensus` and `ivar variants` respectively.
- Add `-a` to the `samtools mpileup` in `merge_mpileup_consensus.nf` so the consensus spans the full reference. Without it `mpileup` emits only covered positions, `ivar consensus` silently truncates to the covered span, and each polishing round remaps to the already-truncated sequence — so terminal bases were lost cumulatively and could never be recovered (4-10 bp per segment on the influenza test set). Uncovered positions are now written as `N` instead of dropped. Note this is `-a`, not the `-aa` removed in v2.1.0-Beta: per-segment BAMs retain the full multi-record header, and `-aa` would emit every segment and concatenate them into one consensus record.
- Behaviour change: because `calculate_coverage.py` divides by the consensus length, previously-truncated positions were excluded from the coverage denominator entirely. With the consensus now spanning the reference, `consensus_coverage_cutoff` measures the fraction of the reference actually resolved, and samples missing terminal sequence score lower than they did before.
- Implement hybrid-capture read deduplication in `deduplicate_reads.nf`. The branch previously raised `error "Hybrid Capture deduplication not yet implemented."`, so any sample with `read_deduplication: positional` (formerly `hybrid-capture`) aborted the run. It now mirrors the `umi` branch (`samtools collate | fixmate -m | sort | markdup -r`), writes `${bam.baseName}.dedup.bam` to match the declared output, and records per-replicate markdup stats. First exercised by the multi-segment run; amplicon paths never reach this branch.
- Pass `-@ ${task.cpus}` to the `samtools` stages in both deduplication branches, which previously ran single-threaded regardless of the allocated CPUs.
- Rename the sample-sheet column and `meta` key `sequencing_tech`, and the param `sequencing_technique`, to `read_deduplication`, and rename its values `mips` to `umi` and `hybrid-capture` to `positional` (`amplicon` is unchanged). The field selects how duplicate reads are handled — `amplicon` clips primers and skips deduplication, `umi` runs `samtools markdup --barcode-name`, `positional` runs coordinate-based `samtools markdup` — so it now names the behaviour rather than the library prep that typically calls for it. Breaking change: sample sheets and params files must be updated. `process_sample_sheet.nf` raises an explicit error for a legacy `sequencing_tech` column or a legacy `mips` / `hybrid-capture` value, rather than silently falling back to the param default.
- Rename `data/renamed_fastq/` to `data/fastq/` and `data/sars_reference/` to `data/reference/` to accommodate a multi-segment test set; update all sample sheets, `nextflow.config` profiles, and nf-test/test configs that referenced the old paths.
- Add `data/primers/nCoV-2019.v3.primer.bed`, the ARTIC v3 scheme. `primers.csv` has always mapped `primer_id` `ARTIC_v3` to this path, but the file was never committed, so any sample sheet using that primer set resolved to a missing BED — `process_sample_sheet.nf` does not check existence, so it would have surfaced as a failure inside `AMPLICON_CLIP` rather than as a sample-sheet error. Taken verbatim from `artic-network/primer-schemes` `nCoV-2019/V3/nCoV-2019.primer.bed`: 218 records on `MN908947.3`, including all 22 alt primers, with the strand column `samtools ampliconclip --strand` requires.
- Remove the iVar primer-mismatch masking workflow, which was never wired into the pipeline: `bwa_mem_filter_sort.nf` and `align_fasta_filter_sort.nf` (two attempts at aligning the primer FASTA to the consensus), `bam_to_bed.nf`, `call_primer_variants.nf` and `ivar_primer_variants.nf` (two attempts at calling variants on the primer BAM), `mask_primers.nf` (`ivar getmasked`), `remove_masked_sort_index.nf` (`ivar removereads`), and `call_masked_variants.nf`. Primer handling stays with `AMPLICON_CLIP` (`samtools ampliconclip`), which is coordinate-based and needs none of this. Also remove the `primer_bed`, `primer_fasta`, and `primer_pairs` params the chain would have consumed, and the data files only that chain consumed: the primer-pair tables `data/primers/SARS-CoV-2.primer_pairs.tsv` and `data/primers/v5.3.2_400_primer_pairs.tsv`, read by `ivar getmasked`, and the primer FASTAs `data/primers/SARS-CoV-2.primer.fa` and `data/primers/v5.3.2_400_primers.fa`, which were the alignment step's query. No `.nf` file read those params — primer BEDs reach `AMPLICON_CLIP` through `primer_csv` and each sample's `primer_id` — and none of the modules were included by `main.nf` or any workflow. `call_primer_variants.nf` additionally referenced `params.variant_min_qual_score`, which is not defined anywhere, so it could not have run under `nextflow.enable.strict`. Document `primer_csv` in the README in place of the removed `primer_bed`.
- Remove `data/primers/primer.5.3.2.bed`, a degenerate-expanded derivative of the ARTIC v5.3.2 scheme built for that masking workflow. It splits the scheme's single IUPAC base (the `R` in `SARS-CoV-2_400_84_RIGHT_2`) into explicit A and G primers and drops the `_400_` name infix, which left its primer names matching neither `v5.3.2_400_primer_pairs.tsv` nor `v5.3.2_400_primers.fa` — `ivar getmasked` joins BED to pair file by name, so it would have masked nothing. No effect on primer trimming: `samtools ampliconclip` clips by interval and both files cover `MN908947.3:26048-26072`. `data/primers/SARs-CoV-2_v5.3.2_400.primer.bed`, the file actually used via `primers.csv`, is byte-identical to upstream `quick-lab/SARS-CoV-2` `400/v5.3.2_400/` and is unaffected.
- Point the `rhino_test` profile at `MN908947.fasta` + `SARS_CoV_2.gff3`. It previously referenced `protein_annotations.gff3`, which is not in the repository, and paired it with `NC_045512.fa`, whose seqid `NC_045512.2` matches neither the GFF nor the profile's primer BED (`MN908947.3`) — a desync of the reference/GFF/BED naming invariant. `MN908947.fasta` aligns all three.
- Record that `SARS_CoV_2.gff3` (NCBI annotwriter, seqid `MN908947.3`) is the only supported SARS-CoV-2 annotation, and that Ensembl's `Sars_cov_2.ASM985889v3.101` must not be substituted for it. Several untracked candidates had accumulated in a working copy of `data/reference/` — a byte-identical duplicate of `SARS_CoV_2.gff3`, the Ensembl annotation, and an AGAT-repaired variant of it — and were discarded rather than committed. The Ensembl annotation models ORF1ab as one continuous CDS `266..21555`, ignoring the -1 ribosomal frameshift at 13468, so every amino-acid call downstream of that point is one base out of frame — against `ivar variants` 1.4.3, position 15000 yields codon `GTT` where the NCBI annotation's two-block CDS correctly yields `AGT`. It also emits duplicate variant rows across the overlapping ORF1a/ORF1ab features, and all 12 of its mRNA `Parent` attributes reference gene IDs absent from the file; adding the missing gene layer with AGAT fixes neither problem. `SARS_CoV_2.gff3` (NCBI annotwriter) remains the only SARS-CoV-2 annotation and is what every profile already used.

## v2.2.0-Beta
- Fix MIPS sequence read deduplication in `deduplicate_reads.nf`.
- Add updated GFF3 file.
- Update CI test, snapshot to be compatible with new GFF3 file, latest version of `nf-test.`
- Adjust depth of coverage process definition to match consensus, variant calling.
- Fix small workflow test case.

## v2.1.0-Beta
- Add missing FASTA indexing to `ivar_variants.nf`.
- Add forced error exitcode if FASTA index is missing in `ivar_variants.nf`.
- Remove `-aa` flag from all calls to `samtools mpileup`.
- Add suffixes to consensus FASTA file names to distinguish rough and polished consensus sequences in `merge_mpileup_consensus.nf`.
- Manually set FASTA header lines in `merge_mpileup_consensus.nf` to make clear whether some output is descended from an alignment to a rough or polished consensus sequence.
- Update `merge_mpileup_consensus.nf` to output a meta hashlist and consensus path, rather than a sample ID and consensus path, to simplify workflow channel operations in `build_consensus.nf`.
- Improve clarity of process call naming scheme in `build_consensus.nf`.
- Fix detection of consensus FASTAs with empty sequence blocks in `build_consensus.nf`.
- Add new process, `sort_index_bam.nf`, which retains unmapped reads when sorting/indexing BAM files.
- Remove redundant BAM filtering from `call_ivar_variants.nf`- this should be handled by removing empty FASTA files earlier in the workflow.
- Add NF-Test continuous integration test case.
- Address broken replicate merging ahead of consensus calling (both rough and polished).
- Update outdated workflow syntax in `call_variants_ivar.nf`.
- Change how relative paths to `primer.csv` are handled in `process_sample_sheet.nf`.
- Add replicate ID tags to workflow processes to aid in debugging and improve feedback to users.

## v2.0.2-Beta
- Match `samtools mpileup` settings in consensus calling, variant calling processes.
- Add consensus polishing and remapping to `build_consensus` subworkflow.
- Rework empty consensus FASTA detection in `build_consensus` subworkflow.
- Address output naming error in Trimmomatic process.
- Add reformatted SARS-CoV-2 GFF3 file with simple ORF names in iVar output.
- Remove NGMerge from read processing subworkflow.

## v2.0.1-Beta
- Remove redundant variables in `nextflow.config` param block.

## v2.0-Beta
- Use `samtools ampliconclip` to remove primer sequences.
- Samplesheet changes:
    - Add optional columns: replicate_id, sequencing_tech, primer_id. All can be supplied default values via a params file, which will be applied to all rows where the optional columns are empty.
    - Add support for parsing interleaved FASTQ files, where there is only 1 FASTQ per row.
- Add primer_csv, allowing for input samplesheets with multiple primer schemes to be run simultaneously. Primers are linked to each row on the input samplesheet by a primer ID.
- Add support for sequencing technologies that require read deduplication (MIPs, Hybrid Capture).
    - Support for Hybrid Capture will come in a future update. Currently, we only support pre-processed MIPs reads with UMIs (UMIs already extracted to read names and trimmed). We plan to add full support for MIPs, and unprocessed MIPs reads, in a future update.
- Add Trimmomatic, NG-Merge to SNPer Docker image to support more thorough paired read quality control.
- Add support for pre-trimmed input reads, indicated by a `primer_id` value of 'None' on an input samplesheet.

## v0.1.0.4-Beta
- Update Docker image to fix Liftoff dependencies that have broken since its last push.

## v0.1.0.3-Beta
- Add BAQ to variant calling to address false positives.
- Set 'convert_tsv_coords.nf' to retry on task failure.
- Update default SNPer parameter scheme based on benchmarking results.

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