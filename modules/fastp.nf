process FASTP {
    tag "${meta.replicate_id}"
    publishDir "${params.output_dir}/fastp/", mode: 'copy', pattern: "*.fastp.{json,html}"

    cpus 4
    memory '2G'
    time '4h'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.replicate_id}_paired.r{1,2}.fastq.gz"), emit: reads
        tuple val(meta), path("${meta.replicate_id}.fastp.json"),             emit: json
        tuple val(meta), path("${meta.replicate_id}.fastp.html"),             emit: html

    script:
    def rep_id = meta.replicate_id
    // Interleaved input was silently broken under the Trimmomatic process this
    // replaces: TrimmomaticPE was handed one file where it needs two. fastp reads
    // interleaved natively and always writes split R1/R2, so READS_QC flips
    // meta.interleaved to false downstream.
    //
    // --detect_adapter_for_pe is deliberately absent from the interleaved branch:
    // it re-opens the read2 *file* to sample it, and in interleaved mode there is
    // no in2, so fastp exits 255 with `ERROR: Failed to open file:`. Adapter
    // removal still happens either way -- for paired data fastp finds adapters by
    // overlap analysis of each pair by default, and that is what does the work
    // here. The flag only adds a second, sequence-based detection pass that helps
    // on non-overlapping pairs.
    def input_args = meta.interleaved
        ? "--interleaved_in --in1 ${reads[0]}"
        : "--in1 ${reads[0]} --in2 ${reads[1]} --detect_adapter_for_pe"
    """
    fastp \\
        ${input_args} \\
        --out1 ${rep_id}_paired.r1.fastq.gz \\
        --out2 ${rep_id}_paired.r2.fastq.gz \\
        --length_required ${params.fastp_min_length} \\
        --thread ${task.cpus} \\
        --json ${rep_id}.fastp.json \\
        --html ${rep_id}.fastp.html
    """
}
