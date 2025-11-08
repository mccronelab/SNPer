process AMPLICONCLIP {
    cpus 1
    memory 2G
    time 12.h

    input:
        tuple val(key), path(sorted_bam), path(bam_index), path(primer_bedfile)

    output:
        tuple val(key), path("*.primertrim.bam")

    script:
        """
        samtools ampliconclip --hard-clip -b ${primer_bedfile} ${sorted_bam} -o ${sorted_bam.simpleName}.primertrim.bam
        """
}