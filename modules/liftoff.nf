process LIFTOFF {
    publishDir "${params.output_dir}/gff3/", mode: 'copy'

    input:
        tuple val(key), path(reference), path(target), path(gff_file)

    output:
        tuple val(key), path("${target.simpleName}.gff3")
        
    script:
    """
    liftoff -g ${gff_file} -o ${target.simpleName}.gff3 ${target} ${reference}
    """
}