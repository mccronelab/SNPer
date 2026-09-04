// Counts mapped reads so build_consensus.nf can drop replicates below
// min_mapped_reads. Emits the count on stdout rather than
// re-declaring the BAM to reduce disk use.
// Important: if your primer scheme does not match your samples,
// ampliconclip will drop almost all reads, which can cause
// samples with deep sequences to fail to pass this filter.

process COUNT_MAPPED_READS {
    label 'process_low'
    tag "${meta.replicate_id}"

    cpus 1
    memory 1G
    time 1.h

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), stdout

    script:
    // -F 4 counts secondary and supplementary alignments too, matching flagstat's
    // ' mapped (' line.
    """
    samtools view -c -F 4 ${bam} | tr -d '\\n'
    """
}