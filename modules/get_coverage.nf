process GET_COVERAGE {
    publishDir "${params.output_dir}/replicate_coverage", mode: 'copy'

    cpus 1
    memory 1G
    time 1.h

    input:
        tuple val(key), path(bam), path(bam_index)

    output:
        path "*.csv"

    script:
        """
        samtools depth -a -d 100000 ${bam} -Q ${params.variant_min_mapQ} > ${bam.baseName}_coverage.csv
        """
}