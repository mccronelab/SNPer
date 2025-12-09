include { CONSENSUS_GEN } from "./workflows/build_consensus"
include { TEST_CONSENSUS } from "./workflows/test_consensus"
include { CALL_VARIANTS_IVAR } from './workflows/call_variants_ivar'
include { TEST_VARIANTS } from './workflows/test_variants'
include { PROCESS_SAMPLE_SHEET } from './workflows/process_sample_sheet'
include { READS_QC } from './workflows/process_reads'

nextflow.enable.dsl=2
// strict mode: https://www.nextflow.io/docs/latest/reference/feature-flags.html#config-feature-flags
nextflow.enable.strict = true


workflow {
    main:
        // import values and files from params
        interleaved_default = params.interleaved
        default_sequencing_tech = params.sequencing_technique
        primer_csv = Channel.fromPath(params.primer_csv)
        primer_default = params.primer_id_default

        input_ch = Channel.fromPath(params.sample_sheet)
        
        samples  = PROCESS_SAMPLE_SHEET(input_ch, primer_csv, interleaved_default, 
            default_sequencing_tech, primer_default) //metadata, [fastq1,fastq2]]

        if(!params.skip_qc) {
            samples = READS_QC(samples)
        }
    
        TEST_CONSENSUS(samples)
        aligned_reads_and_consensus = TEST_CONSENSUS.out.bams_with_consensus
        polished_consensus = TEST_CONSENSUS.out.polished_consensus

        TEST_VARIANTS(aligned_reads_and_consensus, polished_consensus)

        //CALL_VARIANTS_IVAR(aligned_reads_and_consensus)

    publish:
        bams = aligned_reads_and_consensus
}

output {
    bams {
        path 'bams'
    }
}