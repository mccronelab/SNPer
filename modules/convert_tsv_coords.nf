process CONVERT_TSV_COORDS {
    publishDir "${params.output_dir}/reference_coordinate_variants/", mode: 'copy'

    cpus 1
    memory 2G
    time 12.h

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