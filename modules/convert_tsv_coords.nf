process CONVERT_TSV_COORDS {
    publishDir "${params.output_dir}/reference_coordinate_variants/", mode: 'copy'

    cpus 1
    memory 2G
    time 12.h

    input:
        tuple val(key), path(consensus), path(reference), val(variant_tsv)

    output:
        tuple val(key), path("${variant_tsv.simpleName}.ref_coords.tsv")

    script:
        // prevent MAFFT from running into permissions issues on clusters by reassigning $TMPDIR
        """
        export TMPDIR="\$(pwd)/tmp/"
        python3 ${projectDir}/bin/convert_tsv_coords.py ${reference} ${consensus} ${variant_tsv} ${variant_tsv.simpleName}.ref_coords.tsv
        """
}