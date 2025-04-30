process LIFTOFF {
    publishDir "${params.output_dir}/gff3/", mode: 'copy'

    input:
        tuple val(meta), path(reference), path(target), path(gff_file)

    output:
        tuple val(meta), path("${target.simpleName}.gff3")
        
    script:
    """
    liftoff -g ${gff_file} -o ${target.simpleName}.gff3 ${target} ${reference}
    """
}