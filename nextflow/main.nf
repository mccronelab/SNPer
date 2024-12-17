include { CONSENSUS_GEN } from "./workflows/build_consensus.nf"

nextflow.enable.dsl=2

workflow {
    // load values from params
    reference = params.reference
    fastq_dir = params.fastq_dir
    primer_bed = channel.fromPath(params.primer_bedfile)

    // produces a triple channel with key, path(read1), path(read2)
    // key is shared text between base names of path(read1) and path(read2)
    // i.e. [SRR001, SRR001_1.fastq.gz, SRR001_2.fastq.gz]
    techinal_rep_A = channel.fromFilePairs("${fastq_dir}/*-A_{1,2}.fastq.gz")
    techinal_rep_B = channel.fromFilePairs("${fastq_dir}/*-B_{1,2}.fastq.gz")

    //techinal_rep_A.view {key, reads -> "Key: ${key}, Reads: ${reads}"}
    //techinal_rep_B.view()

    consensus_seqs = CONSENSUS_GEN(reference, techinal_rep_A, techinal_rep_B, primer_bed)

}