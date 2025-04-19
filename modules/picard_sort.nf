process PICARD_SORT {
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    cpus 1
    memory { 2G * task.attempt }
    time { 4.h * task.attempt }

    input:
        tuple val(key), path(trimmed_bam)

    output:
        tuple val(key), path("*.removed.primertrim.sorted.bam"), path("*.removed.primertrim.sorted.bai")

    script:
    """
    PicardCommandLine SortSam SO=coordinate INPUT=${trimmed_bam} OUTPUT=${trimmed_bam.simpleName}.removed.primertrim.sorted.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
    """
}