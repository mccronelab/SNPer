process GET_VARIANT_READ_DEPTH {
    label 'process_high'
    publishDir "${params.output_dir}/replicate_coverage", mode: 'copy'
    tag "${meta.replicate_id}"

    cpus 1
    memory 1G
    time 1.h

    input:
        tuple val(meta), path(bam), path(bam_index)

    output:
        path "*.tsv"

    script:
        """
        samtools depth -a -d 0 ${bam} -Q 30 -q ${params.variant_min_mapQ} > ${bam.baseName}_coverage.tsv
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
    python3 ${projectDir}/bin/calculate_coverage.py ${consensus}
    """
    
}