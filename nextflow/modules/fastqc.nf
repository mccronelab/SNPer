process FASTQC {
    publishDir = "${output_dir}/fastqc/", mode: 'copy'
    container = "staphb/fastqc:latest"

    input:
        // this input pattern matches the output of channel.fromFilePairs()
        tuple val(key), path(paired_reads)
    
    output:
        tuple val(key), path("*_1_fastqc.zip"), path("*_2_fastqc.zip")

    script:
    // supplying 2 paths here produces 2 output files
    """
    fastqc --noextract -f fastq ${paired_reads}
    """
}