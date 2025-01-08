process GET_COVERAGE {
    publishDir "${params.output_dir}/consensus_coverage", mode: 'copy'
    container "staphb/samtools:latest"

    input:
        tuple val(key), path(sorted_trimmed_bam), path(bam_index)

    output:
        path "*.csv"

    script:
        """
        samtools depth -a -d 100000 ${sorted_trimmed_bam} > ${sorted_trimmed_bam.baseName}_coverage.csv
        """
}