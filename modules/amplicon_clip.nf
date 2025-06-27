process AMPLICON_CLIP {
    cpus 1
    memory 2G
    time 12.h
    label 'process_medium'

    def clipped = params.remove_unclipped_reads? '--clipped': ''

    input:
        tuple val(meta), path(sorted_bam), path(bam_index), path(reference)

    output:
        tuple val(meta), path("*.primertrim.bam")
        

    script:
    """
    samtools ampliconclip \
        -b ${meta.primer_bedfile} \
        --hard-clip \
        --strand \
        ${clipped} \
        -f ${sorted_bam.simpleName}.clipped.log \
        --primer-counts ${sorted_bam.simpleName}.primers.txt \
        -u \
        ${sorted_bam}\
    | samtools collate  -O -u - -T ./${meta.sample}.TEMP \
    | samtools fixmate  -m -u - - \
    | samtools sort  -u - \
    | samtools calmd - ${reference} -b > ${sorted_bam.simpleName}.primertrim.bam
    """
}