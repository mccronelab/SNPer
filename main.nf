include { CONSENSUS_GEN } from "./workflows/build_consensus"
include { CALL_VARIANTS_IVAR } from './workflows/call_variants_ivar'
include { PROCESS_SAMPLE_SHEET } from './workflows/process_sample_sheet'
include { READS_QC } from './workflows/process_reads'

nextflow.enable.dsl=2
// strict mode: https://www.nextflow.io/docs/latest/reference/feature-flags.html#config-feature-flags
nextflow.enable.strict = true


workflow {
    // import values and files from params
    interleaved_default = params.interleaved
    default_read_deduplication = params.read_deduplication
    primer_csv = Channel.fromPath(params.primer_csv)
    primer_default = params.primer_id_default

    input_ch = Channel.fromPath(params.sample_sheet)
    
    samples  = PROCESS_SAMPLE_SHEET(input_ch, primer_csv, interleaved_default, 
        default_read_deduplication, primer_default) //metadata, [fastq1,fastq2]]

    if(!params.skip_qc) {
        samples = READS_QC(samples)
    }
 
    // tuple (consensus_name, consensus.fasta)
    CONSENSUS_GEN(samples)

    CALL_VARIANTS_IVAR(CONSENSUS_GEN.out.variant_bams, CONSENSUS_GEN.out.consensus)
}
