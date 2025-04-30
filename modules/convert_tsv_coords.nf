process CONVERT_TSV_COORDS {
    label 'process_medium'
    publishDir "${params.output_dir}/variants/", mode: 'copy'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    cpus {1 * task.attempt}
    // need to account for potentially increasing CPU allocation
    memory { 2G * task.attempt * task.cpus}
    time { 4.h * task.attempt }

    input:
        tuple val(meta), path(consensus), path(reference), path(variant_tsv)

    output:
        tuple val(meta), path("${variant_tsv.simpleName}.ref_coords.tsv")

    script:
        // prevent MAFFT from running into permissions issues on clusters by reassigning $TMPDIR
        """
        export TMPDIR="\$(pwd)/tmp/"
        python3 ${projectDir}/bin/convert_tsv_coords.py ${reference} ${consensus} ${variant_tsv} ${variant_tsv.simpleName}.ref_coords.tsv
        """
}