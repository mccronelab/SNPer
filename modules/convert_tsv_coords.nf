process CONVERT_TSV_COORDS {
    label 'process_medium'
    publishDir "${params.output_dir}/variants/", mode: 'copy'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    tag "${meta.replicate_id}"

    cpus {1 * task.attempt}
    // need to account for potentially increasing CPU allocation
    memory { 2G * task.attempt * task.cpus}
    time { 4.h * task.attempt }

    input:
        tuple val(meta), path(consensus), path(reference), path(variant_tsv)

    output:
        tuple val(meta), path("${meta.replicate_id}${meta.segment_label}.ref_coords.tsv")

    script:
        // prevent MAFFT from running into permissions issues on clusters by reassigning $TMPDIR
        //
        // convert_tsv_coords.py assumes a single-record reference: it MAFFT-aligns
        // reference record[0] against consensus record[1]. A multi-segment reference
        // would put two *reference* segments in those slots and lift against the wrong
        // one. The consensus is already single-record (post BAM split), so extract just
        // this segment's reference record and feed that. meta.segment == the reference
        // record ID by construction (the split keys on reference name), so faidx by name
        // is exact.
        """
        export TMPDIR="\$(pwd)/tmp/"
        samtools faidx ${reference}
        samtools faidx ${reference} "${meta.segment}" > segment_reference.fa
        python3 ${projectDir}/bin/convert_tsv_coords.py segment_reference.fa ${consensus} ${variant_tsv} ${meta.replicate_id}${meta.segment_label}.ref_coords.tsv
        """
}