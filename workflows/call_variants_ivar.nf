include { BWA_MEM  } from "../modules/bwa_mem"
include { FILTER_SORT_INDEX  } from "../modules/filter_sort_index"
include { IVAR_TRIM  } from "../modules/ivar_trim"
include {IVAR_VARIANTS} from '../modules/ivar_variants.nf'

workflow CALL_VARIANTS_IVAR {
    take:
        samples_with_consensus // [key, [fastq1,fastq2], consensus] for each line in the input csv


    main:
        primer_bed = file(params.primer_bed)
        reference_gff = file(params.reference_gff)
        variants = channel.empty()

    // align to consensus
    // trim primers
    // call variants
        // check file is at least 1Kb in size
        filtered_samples = samples_with_consensus.filter { _key, _reads, consensus -> consensus.size() >= 1000 }
        consensus = filtered_samples.map{key,_fastqs,consensus->tuple(key,consensus)}.unique{d -> d[0]}
        BWA_MEM(filtered_samples)
        |  FILTER_SORT_INDEX
        | map { key, sortedBam, bamIndex -> tuple(key,sortedBam,bamIndex,primer_bed)}
        | filter { _key, bam, _index, _bed -> bam.size() >= 1000 } 
        | IVAR_TRIM  //  tuple val(key), path("*.primertrim.bam")
        | combine(consensus,by:0) //  tuple val(key), path("*.primertrim.bam") consensus.fa
        | map{key, bam, consensus_fa -> tuple(key,bam,consensus_fa,reference_gff)}
        | IVAR_VARIANTS
        | set{variants}

    emit:
        variants 
}