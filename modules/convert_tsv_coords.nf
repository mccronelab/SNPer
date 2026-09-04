process CONVERT_TSV_COORDS {
    label 'process_low'
    publishDir "${params.output_dir}/variants/", mode: 'copy'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    tag "${meta.replicate_id}"

    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.h * task.attempt }

    input:
        tuple val(meta), path(consensus), path(variant_tsv), path(aligned_fasta), path(nextclade_tsv)

    output:
        tuple val(meta), path("${meta.replicate_id}${meta.segment_label}.ref_coords.tsv")

    script:
        // NEXTCLADE' alignment spans every sample in the segment; the script picks its row by the
        // consensus record ID and the reference row by --reference-id. meta.segment is the
        // reference record ID by construction (the BAM split keys on reference name).
        """
        convert_tsv_coords.py ${aligned_fasta} ${nextclade_tsv} ${consensus} ${variant_tsv} \\
            ${meta.replicate_id}${meta.segment_label}.ref_coords.tsv \\
            --reference-id "${meta.segment}"
        """
}