include { CONSENSUS_GEN } from "./workflows/build_consensus"
include { CALL_VARIANTS_IVAR } from './workflows/call_variants_ivar'
include { TRIM_AND_MASK } from './workflows/trim_and_mask'
include { PROCESS_SAMPLE_SHEET } from './workflows/process_sample_sheet'

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
    aligned_reads_and_consensus = CONSENSUS_GEN(samples)

    if (params.tiled_amplicons == true) {
        TRIM_AND_MASK(aligned_reads_and_consensus)
          | CALL_VARIANTS_IVAR
    }
    else {
        CALL_VARIANTS_IVAR(aligned_reads_and_consensus)
    }
}
