process REMOVE_MASKED_SORT_INDEX {

    input:
        tuple val(key), path(primer_bed_file), path(mask_txt), path(bam_file), path(index)

    output:
        tuple val(key), path("${bam_file.simpleName}_masked.bam"), path("${bam_file.simpleName}_masked.bam.bai")

    script:
    """
        ivar removereads -i ${bam_file} -p ${bam_file.simpleName}.masked -t ${mask_txt} -b ${primer_bed_file}  
        samtools sort -T ${bam_file.simpleName}.tmp -o ${bam_file.simpleName}_masked.bam ${bam_file.simpleName}.masked.bam        
        samtools index ${bam_file.simpleName}_masked.bam
    """
}