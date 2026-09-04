// One instance per segment, spanning every sample. `--input-ref` takes a single-record
// FASTA, so this fans out by segment and collects across samples -- 8 invocations for
// influenza, 1 for SARS-CoV-2, keeping N=1 the special case of the same path.
//
// Batching is a packaging choice, not a correctness one: nextclade aligns each query
// pairwise against the reference and stacks the rows on reference columns, so per-sample
// results are identical whether or not they share an invocation.
process NEXTCLADE {
    label 'process_medium'
    publishDir "${params.output_dir}/msa/", mode: 'copy', pattern: "*.aligned.fasta"
    publishDir "${params.output_dir}/msa/", mode: 'copy', pattern: "*.nextclade.tsv"
    // flatten per_sample_gff/ so these land beside the names LIFTOFF used to publish
    publishDir "${params.output_dir}/gff3/", mode: 'copy',
        pattern: "per_sample_gff/*.gff3",
        saveAs: { filename -> file(filename).name }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    tag "${segment}"

    cpus { 2 * task.attempt }
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    input:
        tuple val(segment), path(consensus_fas), path(reference), path(reference_annotation)

    output:
        tuple val(segment), path("per_sample_gff/*.gff3"),      emit: segment_gffs
        tuple val(segment), path("${segment}.aligned.fasta"),   emit: msa
        tuple val(segment), path("${segment}.nextclade.tsv"),   emit: nextclade_tsv

    script:
        // The annotation must be subset to this segment: nextclade applies every CDS in
        // --input-annotation to the one reference it was given, and iVar ignores the GFF
        // seqid column entirely, so foreign segments' CDSes would be projected onto the
        // wrong coordinates. Dropping the pragmas takes ##sequence-region with them --
        // it declares NCBI lengths that need not match this reference.
        """
        samtools faidx ${reference}
        samtools faidx ${reference} "${segment}" > segment_reference.fa

        awk -F'\\t' '/^##gff-version/ {print; next} /^#/ {next} \$1 == "${segment}"' \\
            ${reference_annotation} > segment_annotation.gff

        if [ "\$(grep -vc '^#' segment_annotation.gff)" -eq 0 ]; then
            echo "no annotation rows for segment '${segment}' -- the reference FASTA headers and GFF seqids have desynced" >&2
            exit 1
        fi

        # --include-reference puts the reference in the MSA as a gap-free row 1
        # --in-order holds the query rows in input order, so the caller sorts the
        # collected FASTAs to keep the published MSA stable.
        nextclade run \\
            --input-ref segment_reference.fa \\
            --input-annotation segment_annotation.gff \\
            --output-fasta ${segment}.aligned.fasta \\
            --output-tsv ${segment}.nextclade.tsv \\
            --output-annotation-gff ${segment}.annotation.gff \\
            --include-reference \\
            --in-order true \\
            --jobs ${task.cpus} \\
            ${consensus_fas}

        mkdir -p per_sample_gff
        split_nextclade_gff.py ${segment}.annotation.gff ${segment}.nextclade.tsv \\
            segment_annotation.gff per_sample_gff ${consensus_fas}
        """
}
