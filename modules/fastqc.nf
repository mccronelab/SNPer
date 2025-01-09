process FASTQC {
    publishDir "${params.output_dir}/fastqc/", mode: 'copy'

    input:
        // this input pattern matches the output of channel.fromFilePairs()
        tuple val(key), path(paired_reads)
    
    output:
        tuple val(key), path("*.html"), emit: html
        tuple val(key), path("*.zip") , emit: zip

    script:
        // supplying 2 paths here produces 2 output files
        """
        fastqc --noextract -f fastq ${paired_reads}
        """
}