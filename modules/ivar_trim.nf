process IVAR_TRIM {
    label 'process_high'
    cpus 1
    memory 2G
    time 12.h

    input:
        tuple val(key), path(sorted_bam), path(bam_index), path(primer_bed)

    output:
        tuple val(key), path("*.primertrim.bam")

    script:
    """
    ivar trim -i ${sorted_bam} -b ${primer_bed} -p ${sorted_bam.simpleName}.primertrim.bam
    """
}