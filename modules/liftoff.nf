process LIFTOFF {
    conda 'bioconda::liftoff'

    input:
        tuple val(key), path(reference), path(target), path(gff_file)

    output:
        tuple val(key), path("${target.simpleName}.gff3")
    script:
    """
    liftoff -g ${gff_file} -o ${target.simpleName}.gff3 ${reference} ${target}
    """
}