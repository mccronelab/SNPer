include { CONSENSUS_GEN } from "./workflows/build_consensus"
include { CALL_VARIANTS_IVAR } from './workflows/call_variants_ivar'

nextflow.enable.dsl=2
// strict mode: https://www.nextflow.io/docs/latest/reference/feature-flags.html#config-feature-flags
nextflow.enable.strict = true

workflow {
    // load values from params
    reference_fasta = params.reference_fasta
    reference_gff = params.reference_gff
    fastq_dir = params.fastq_dir
    primer_bed = file(params.primer_bedfile)
    primer_fasta = file(params.primer_fasta)
    primer_pairs = file(params.primer_pair_tsv)

    // produces a triple channel with key, path(read1), path(read2)
    // key is shared text between base names of path(read1) and path(read2)
    // i.e. [SRR001, SRR001_1.fastq.gz, SRR001_2.fastq.gz]
    techinal_rep_A = channel.fromFilePairs("${fastq_dir}/*-A_{1,2}.fastq.gz")
    techinal_rep_B = channel.fromFilePairs("${fastq_dir}/*-B_{1,2}.fastq.gz")

    // tuple (consensus_name, consensus.fasta)
    consensus_seqs = CONSENSUS_GEN(reference_fasta, techinal_rep_A, techinal_rep_B, primer_bed)

    CALL_VARIANTS_IVAR(reference_fasta, reference_gff, primer_bed, primer_pairs, primer_fasta, consensus_seqs)

}