process LIFTOFF {
    publishDir "${params.output_dir}/gff3/", mode: 'copy'
    container 'staphb/liftoff:latest'

    input:
        tuple val(key), path(reference), path(target), path(gff_file)

    output:
        tuple val(key), path("${target.simpleName}.gff3")
    script:
    """
    MAMBA_SKIP_ACTIVATE=False
    source /usr/local/bin/_activate_current_env.sh
    liftoff -g ${gff_file} -o ${target.simpleName}.gff3 ${target} ${reference}
    """
}