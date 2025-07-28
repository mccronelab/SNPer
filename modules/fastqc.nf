process FASTQC {
    label 'process_high'
    publishDir "${params.output_dir}/fastqc/", mode: 'copy'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(meta), path(paired_reads)
    
    output:
        tuple val(meta), path("*.html"), emit: html
        tuple val(meta), path("*.zip") , emit: zip

    script:
        // supplying 2 paths here produces 2 output files
        """
        fastqc --noextract -f fastq ${paired_reads}
        """
}