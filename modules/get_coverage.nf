process GET_VARIANT_READ_DEPTH {
    label 'process_high'
    publishDir "${params.output_dir}/replicate_coverage", mode: 'copy'
    tag "${meta.replicate_id}"

    cpus 1
    memory 1G
    time 1.h

    input:
        tuple val(meta), path(bam), path(bam_index), path(consensus)

    output:
        path "*.tsv"

    script:
    // Same mpileup as ivar_variants.nf, so the published depth is the depth iVar saw.
    // This is mpileup, so -Q is the base-quality floor.
        """
        samtools faidx ${consensus}
        samtools mpileup -a -d 0 -Q ${params.variant_min_baseQ} -q ${params.variant_min_mapQ} -f ${consensus} ${bam} \
            | cut -f1,2,4 > ${bam.baseName}_coverage.tsv
        """
}

process GET_CONSENSUS_COVERAGE {
    label 'process_low'
    tag "${sample}_${segment}"

    input:
        tuple val(sample), val(segment), path(consensus)
    output:
        tuple val(sample), val(segment), path(consensus), stdout
    script:
    """
    calculate_coverage.py ${consensus}
    """
    
}