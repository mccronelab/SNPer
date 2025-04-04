process CONVERT_TSV_COORDS {
    publishDir "${params.output_dir}/reference_coordinate_variants/", mode: 'copy'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    cpus {1 * task.attempt}
    // need to account for potentially increasing CPU allocation
    memory { 2G * task.attempt * task.cpus}
    time { 4.h * task.attempt }

    input:
        tuple val(key), path(consensus), val(reference), val(variant_tsv)

    output:
        tuple val(key), path("${variant_tsv.simpleName}.ref_coords.tsv")

    script:
        """
        cat ${reference} ${consensus} > ref_and_target.fa
        python3 ${projectDir}/bin/convert_tsv_coords.py ref_and_target.fa ${variant_tsv} ${variant_tsv.simpleName}.ref_coords.tsv
        """
}