include { CONSENSUS_GEN } from "./workflows/build_consensus"
include { CALL_VARIANTS_IVAR } from './workflows/call_variants_ivar'
include { TRIM_AND_MASK } from './workflows/trim_and_mask'
include { PROCESS_SAMPLE_SHEET } from './workflows/process_sample_sheet'

nextflow.enable.dsl=2
// strict mode: https://www.nextflow.io/docs/latest/reference/feature-flags.html#config-feature-flags
nextflow.enable.strict = true


workflow {
    // import values and files from params
    interleaved_default = params.interleaved
    default_sequencing_tech = params.sequencing_technique
    primer_csv = Channel.fromPath(params.primer_csv)
    primer_default = params.primer_id_default

    input_ch = Channel.fromPath(params.sample_sheet)
    
    samples  = PROCESS_SAMPLE_SHEET(input_ch, primer_csv, interleaved_default, 
        default_sequencing_tech, primer_default) //metadata, [fastq1,fastq2]]
 
    // tuple (consensus_name, consensus.fasta)
    aligned_reads_and_consensus = CONSENSUS_GEN(samples)

    CALL_VARIANTS_IVAR(aligned_reads_and_consensus)
}
