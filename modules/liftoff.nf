process LIFTOFF {
    label 'process_medium'
    publishDir "${params.output_dir}/gff3/", mode: 'copy'
    tag "${sample}_${segment}"

    input:
        tuple val(sample), val(segment), path(reference), path(target), path(gff_file)

    output:
        tuple val(sample), val(segment), path("${target.simpleName}.gff3")
        
    script:
    """
    liftoff -g ${gff_file} -o ${target.simpleName}.gff3 ${target} ${reference}
    """
}