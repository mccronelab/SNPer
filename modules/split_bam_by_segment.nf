// Split one per-replicate BAM into per-segment BAMs, one per reference name
// that has mapped reads. This unifies single- and multi-segment paths: N
// reference records yield N per-segment BAMs, and single-segment is simply N=1.
//
// Each segment's BAM is written under segments/<refname>/ so the caller can
// recover the segment from the parent directory name (robust to '.' in names
// like MN908947.3, which trip up baseName/simpleName extension stripping).
//
// The region list ends in '*' so each per-segment BAM also carries the unplaced
// reads. BWA_REMAP rebuilds FASTQ from this BAM, so dropping '*' would remove
// their only chance to align against the sample's own consensus. They are inert
// in mpileup and samtools depth, so consensus and coverage are unchanged.

process SPLIT_BAM_BY_SEGMENT {
    label 'process_low'
    cpus 1
    memory 2G
    time 12.h
    tag "${meta.replicate_id}"

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("segments/*/*.split.bam"), path("segments/*/*.split.bam.bai")

    script:
    """
    mkdir -p segments
    for ref in \$(samtools idxstats ${bam} | awk '\$3 > 0 {print \$1}'); do
        mkdir -p "segments/\${ref}"
        samtools view -b ${bam} "\${ref}" '*' > "segments/\${ref}/${meta.replicate_id}.\${ref}.split.bam"
        samtools index "segments/\${ref}/${meta.replicate_id}.\${ref}.split.bam"
    done
    """
}
