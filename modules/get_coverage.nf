process GET_VARIANT_READ_DEPTH {
    publishDir "${params.output_dir}/replicate_coverage", mode: 'copy'

    cpus 1
    memory 1G
    time 1.h

    input:
        tuple val(key), path(bam), path(bam_index)

    output:
        path "*.tsv"

    script:
        """
        samtools depth -a -d 100000 ${bam} -Q ${params.variant_min_mapQ} > ${bam.baseName}_coverage.tsv
        """
}

process GET_CONSENSUS_COVERAGE {
    input:
        tuple val(meta), path(consensus)
    output:
        tuple val(meta), path(consensus), stdout
    script:
    """
    python3 ${projectDir}/bin/calculate_coverage.py ${consensus}
    """
    
}