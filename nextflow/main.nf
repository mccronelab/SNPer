include { CONSENSUS_GEN } from "./workflows/build_consensus"
include { CALL_VARIANTS_IVAR } from './workflows/call_variants_ivar'
include {PROCESS_SAMPLE_SHEET} from './workflows/process_sample_sheet'

nextflow.enable.dsl=2
// strict mode: https://www.nextflow.io/docs/latest/reference/feature-flags.html#config-feature-flags
nextflow.enable.strict = true

workflow {
    // params.reference_gff = file(params.reference_gff)
    // params.primer_bed = file(params.primer_bedfile)
    // params.primer_fasta = file(params.primer_fasta)
    // params.primer_pairs = file(params.primer_pair_tsv)
    // params.primer_info = file(params.primer_info_tsv)


    input_ch = Channel.fromPath(params.sample_sheet)
    
    samples  = PROCESS_SAMPLE_SHEET(input_ch) //[key, [fastq1,fastq2]]

    // tuple (consensus_name, consensus.fasta)
    CONSENSUS_GEN(samples).view()

    // CALL_VARIANTS_IVAR(reference_fasta, reference_gff, primer_bed, primer_pairs, primer_fasta, primer_info, consensus_seqs)

}
