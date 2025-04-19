process CUTADAPT {
    input:
        tuple val(key), val(replicate), path(paired_reads)

    output:
        tuple val(key), val(replicate), path("*.cut.fastq.gz")
    script:
        """
        cutadapt -a ${params.adapter_fw} -A ${params.adapter_rev} -o ${paired_reads[0].simpleName}.cut.fastq.gz -p ${paired_reads[1].simpleName}.cut.fastq.gz ${paired_reads}
        """
}