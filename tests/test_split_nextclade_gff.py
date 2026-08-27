#!/usr/bin/env python3
# test_split_nextclade_gff.py
# Description: Unit tests for bin/split_nextclade_gff.py.
# Author: Conner Copeland
# Contact: ccopelan@fredhutch.org
# Created 2026-08-26

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "split_nextclade_gff.py"


def _load_script():
    """Imports the bin/ script by path; it is a CLI, not an installed module."""
    spec = importlib.util.spec_from_file_location("split_nextclade_gff", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


sng = _load_script()


# ---- Fixture builders ----
#
# Attribute strings mirror what nextclade emits: per-run `CDS-<i>-` / `Gene-<i>-` ID
# prefixes and a `seq_index` bookkeeping key, both of which must not survive the split.

REFERENCE_ANNOTATION = "\n".join(
    [
        "##gff-version 3",
        "##sequence-region HA 1 1775",
        "HA\tRefSeq\tgene\t33\t1730\t.\t+\t.\tID=gene-FLUAVs4gp1;Name=HA",
        "HA\tRefSeq\tCDS\t33\t1730\t.\t+\t0\tID=cds-NP_040980.1;Parent=gene-FLUAVs4gp1;Name=HA",
    ]
) + "\n"


def combined_gff(*, queries=("S1_HA_polished", "S2_HA_polished"), attribute_order="natural"):
    """
    Builds a combined `--output-annotation-gff` for the given queries.

    attribute_order flips the emitted key order to stand in for nextclade's
    hash-seeded attribute serialization, which permutes between runs.
    """
    lines = ["##gff-version 3"]

    for index, query in enumerate(queries):
        if attribute_order == "natural":
            gene_attrs = f"ID=Gene-{index}-gene-FLUAVs4gp1;Name=HA;seq_index={index}"
            cds_attrs = (
                f"ID=CDS-{index}-cds-NP_040980.1;Parent=Gene-{index}-gene-FLUAVs4gp1;"
                f"Name=HA;gene=HA;seq_index={index}"
            )
        else:
            gene_attrs = f"seq_index={index};Name=HA;ID=Gene-{index}-gene-FLUAVs4gp1"
            cds_attrs = (
                f"gene=HA;seq_index={index};Name=HA;"
                f"Parent=Gene-{index}-gene-FLUAVs4gp1;ID=CDS-{index}-cds-NP_040980.1"
            )

        lines.append(f"##sequence-region {query} 1 1775")
        lines.append(f"{query}\tnextclade\tgene\t33\t1730\t.\t+\t.\t{gene_attrs}")
        lines.append(f"{query}\tnextclade\tCDS\t33\t1730\t.\t+\t0\t{cds_attrs}")

    return "\n".join(lines) + "\n"


def nextclade_tsv(rows):
    """Builds a nextclade TSV. `rows` maps seqName to a dict of column overrides."""
    columns = ["seqName", "clade", "errors", "warnings", "failedCdses"]
    lines = ["\t".join(columns)]

    for seq_name, overrides in rows.items():
        record = {"seqName": seq_name, "clade": "A", "errors": "", "warnings": "", "failedCdses": ""}
        record.update(overrides)
        lines.append("\t".join(record[column] for column in columns))

    return "\n".join(lines) + "\n"


@pytest.fixture
def workdir(tmp_path):
    """
    Lays out a two-query run: reference annotation, combined GFF, TSV, consensus FASTAs.

    The FASTA filenames deliberately do NOT match their record IDs, so any test that
    passes proves the record ID was read from the file rather than parsed from the name.
    """
    (tmp_path / "reference.gff").write_text(REFERENCE_ANNOTATION, encoding="utf8")
    (tmp_path / "combined.gff").write_text(combined_gff(), encoding="utf8")
    (tmp_path / "nextclade.tsv").write_text(
        nextclade_tsv({"S1_HA_polished": {}, "S2_HA_polished": {}}), encoding="utf8"
    )

    fastas = []
    for index, query in enumerate(("S1_HA_polished", "S2_HA_polished"), start=1):
        path = tmp_path / f"query{index}.fa"
        path.write_text(f">{query}\nACGTACGTAC\n", encoding="utf8")
        fastas.append(path)

    outdir = tmp_path / "per_sample_gff"
    outdir.mkdir()

    return {
        "root": tmp_path,
        "reference": tmp_path / "reference.gff",
        "combined": tmp_path / "combined.gff",
        "tsv": tmp_path / "nextclade.tsv",
        "fastas": fastas,
        "outdir": outdir,
    }


def run_cli(workdir, **overrides):
    """Runs the script as Nextflow does, returning the CompletedProcess."""
    argv = [
        sys.executable,
        str(SCRIPT),
        str(overrides.get("combined", workdir["combined"])),
        str(overrides.get("tsv", workdir["tsv"])),
        str(overrides.get("reference", workdir["reference"])),
        str(overrides.get("outdir", workdir["outdir"])),
        *[str(path) for path in overrides.get("fastas", workdir["fastas"])],
    ]
    return subprocess.run(argv, capture_output=True, text=True)


# ---- normalize_attributes ----


def test_normalize_strips_per_run_id_prefixes():
    """The CDS-<i>-/Gene-<i>- prefixes are what would churn iVar's GFF_FEATURE."""
    result = sng.normalize_attributes(
        "ID=CDS-3-cds-NP_040980.1;Parent=Gene-3-gene-FLUAVs4gp1;Name=HA"
    )
    assert result == "ID=cds-NP_040980.1;Parent=gene-FLUAVs4gp1;Name=HA"


def test_normalize_drops_seq_index():
    assert "seq_index" not in sng.normalize_attributes("ID=CDS-0-x;seq_index=0;Name=HA")


def test_normalize_orders_id_parent_name_first_then_sorted():
    result = sng.normalize_attributes("zeta=1;Name=HA;alpha=2;Parent=Gene-0-g;ID=CDS-0-c")
    assert result == "ID=c;Parent=g;Name=HA;alpha=2;zeta=1"


def test_normalize_is_order_independent():
    """The determinism guarantee: permuted input must produce identical output."""
    forward = sng.normalize_attributes("ID=CDS-0-c;Parent=Gene-0-g;Name=HA;gene=HA;seq_index=0")
    reverse = sng.normalize_attributes("seq_index=0;gene=HA;Name=HA;Parent=Gene-0-g;ID=CDS-0-c")
    assert forward == reverse


def test_normalize_leaves_unrelated_keys_untouched():
    """iVar reads gene and product straight out of column 9; they must survive verbatim."""
    result = sng.normalize_attributes("ID=CDS-0-c;gene=HA;product=haemagglutinin")
    assert "gene=HA" in result and "product=haemagglutinin" in result


def test_normalize_tolerates_trailing_semicolon():
    assert sng.normalize_attributes("ID=CDS-0-c;Name=HA;") == "ID=c;Name=HA"


def test_normalize_keeps_equals_signs_inside_values():
    """partition() splits on the first '=' only, so a value may contain more."""
    assert sng.normalize_attributes("note=a=b") == "note=a=b"


def test_normalize_only_strips_prefix_at_string_start():
    """The regex is anchored; an interior 'CDS-1-' is part of the real ID."""
    assert sng.normalize_attributes("ID=cds-CDS-1-tail") == "ID=cds-CDS-1-tail"


# ---- parse_gff ----


def test_parse_gff_groups_by_seqid(tmp_path):
    path = tmp_path / "combined.gff"
    path.write_text(combined_gff(), encoding="utf8")

    features, regions = sng.parse_gff(str(path))

    assert set(features) == {"S1_HA_polished", "S2_HA_polished"}
    assert len(features["S1_HA_polished"]) == 2
    assert regions["S1_HA_polished"] == "##sequence-region S1_HA_polished 1 1775"


def test_parse_gff_normalizes_attributes_in_place(tmp_path):
    path = tmp_path / "combined.gff"
    path.write_text(combined_gff(queries=("S1_HA_polished",)), encoding="utf8")

    features, _ = sng.parse_gff(str(path))
    cds_line = [line for line in features["S1_HA_polished"] if line.split("\t")[2] == "CDS"][0]

    assert cds_line.split("\t")[8] == "ID=cds-NP_040980.1;Parent=gene-FLUAVs4gp1;Name=HA;gene=HA"


def test_parse_gff_skips_comments_and_blank_lines(tmp_path):
    path = tmp_path / "combined.gff"
    path.write_text(
        "##gff-version 3\n\n# a comment\n"
        "S1\tnextclade\tCDS\t1\t9\t.\t+\t0\tID=CDS-0-c\n",
        encoding="utf8",
    )

    features, _ = sng.parse_gff(str(path))
    assert features == {"S1": ["S1\tnextclade\tCDS\t1\t9\t.\t+\t0\tID=c"]}


def test_parse_gff_rejects_malformed_line(tmp_path):
    """A short row means the file is not GFF3; guessing at it would corrupt coordinates."""
    path = tmp_path / "broken.gff"
    path.write_text("##gff-version 3\nS1\tnextclade\tCDS\t1\n", encoding="utf8")

    with pytest.raises(SystemExit) as excinfo:
        sng.parse_gff(str(path))

    assert "malformed GFF3 line" in str(excinfo.value)


# ---- cds_ids ----


def test_cds_ids_ignores_non_cds_features():
    lines = [
        "S1\tnextclade\tgene\t1\t9\t.\t+\t.\tID=gene-a",
        "S1\tnextclade\tCDS\t1\t9\t.\t+\t0\tID=cds-a",
    ]
    assert sng.cds_ids(lines) == {"cds-a"}


def test_cds_ids_dedups_joined_cds():
    """A spliced CDS (PA-X, M2, NEP) appears as multiple rows sharing one ID."""
    lines = [
        "S1\tnextclade\tCDS\t25\t597\t.\t+\t0\tID=cds-PA-X",
        "S1\tnextclade\tCDS\t599\t784\t.\t+\t0\tID=cds-PA-X",
    ]
    assert sng.cds_ids(lines) == {"cds-PA-X"}


# ---- read_query_ids ----


def test_read_query_ids_reads_record_id_not_filename(tmp_path):
    path = tmp_path / "arbitrary_name.fa"
    path.write_text(">S1_HA_polished\nACGT\n", encoding="utf8")

    assert sng.read_query_ids([str(path)]) == {"S1_HA_polished": str(path)}


def test_read_query_ids_rejects_duplicate_record_ids(tmp_path):
    first = tmp_path / "a.fa"
    second = tmp_path / "b.fa"
    first.write_text(">S1_HA_polished\nACGT\n", encoding="utf8")
    second.write_text(">S1_HA_polished\nACGT\n", encoding="utf8")

    with pytest.raises(SystemExit) as excinfo:
        sng.read_query_ids([str(first), str(second)])

    assert "duplicate consensus record ID" in str(excinfo.value)


def test_read_query_ids_rejects_multi_record_fasta(tmp_path):
    """SeqIO.read enforces the one-consensus-per-file assumption the split relies on."""
    path = tmp_path / "two.fa"
    path.write_text(">A\nACGT\n>B\nACGT\n", encoding="utf8")

    with pytest.raises(ValueError):
        sng.read_query_ids([str(path)])


# ---- read_nextclade_tsv ----


def test_read_nextclade_tsv_keys_on_seqname(tmp_path):
    path = tmp_path / "nc.tsv"
    path.write_text(nextclade_tsv({"S1": {"clade": "A"}, "S2": {"clade": "B"}}), encoding="utf8")

    rows = sng.read_nextclade_tsv(str(path))

    assert set(rows) == {"S1", "S2"}
    assert rows["S2"]["clade"] == "B"


def test_read_nextclade_tsv_preserves_quote_characters(tmp_path):
    """QUOTE_NONE: nextclade's free-text columns are unquoted and may contain quotes."""
    path = tmp_path / "nc.tsv"
    path.write_text(nextclade_tsv({"S1": {"warnings": 'CDS "HA" is short'}}), encoding="utf8")

    assert sng.read_nextclade_tsv(str(path))["S1"]["warnings"] == 'CDS "HA" is short'


# ---- write_split_gffs / end to end ----


def test_writes_one_gff_per_query(workdir):
    result = run_cli(workdir)

    assert result.returncode == 0, result.stderr
    written = sorted(path.name for path in workdir["outdir"].iterdir())
    assert written == ["S1_HA_polished.gff3", "S2_HA_polished.gff3"]


def test_output_carries_header_and_sequence_region(workdir):
    run_cli(workdir)
    lines = (workdir["outdir"] / "S1_HA_polished.gff3").read_text(encoding="utf8").splitlines()

    assert lines[0] == "##gff-version 3"
    assert lines[1] == "##sequence-region S1_HA_polished 1 1775"
    assert lines[2].split("\t")[0] == "S1_HA_polished"


def test_samples_share_identical_feature_attributes(workdir):
    """
    The point of stripping seq_index and the ID prefixes: two samples' GFFs must differ
    only in the seqid column, so adding a sample cannot change another's GFF_FEATURE.
    """
    run_cli(workdir)

    def attributes(name):
        text = (workdir["outdir"] / name).read_text(encoding="utf8")
        return [line.split("\t")[8] for line in text.splitlines() if not line.startswith("#")]

    assert attributes("S1_HA_polished.gff3") == attributes("S2_HA_polished.gff3")


def test_output_is_stable_across_attribute_permutation(workdir):
    """nextclade permutes column 9 between runs; the split output must not."""
    run_cli(workdir)
    natural = (workdir["outdir"] / "S1_HA_polished.gff3").read_bytes()

    workdir["combined"].write_text(combined_gff(attribute_order="permuted"), encoding="utf8")
    run_cli(workdir)
    permuted = (workdir["outdir"] / "S1_HA_polished.gff3").read_bytes()

    assert natural == permuted


def test_fails_when_a_query_lost_all_features(workdir):
    """A query nextclade dropped must stop the run, not vanish from the output."""
    workdir["combined"].write_text(combined_gff(queries=("S1_HA_polished",)), encoding="utf8")

    result = run_cli(workdir)

    assert result.returncode != 0
    assert "S2_HA_polished: nextclade produced no annotation features" in result.stderr


def test_fails_when_a_query_has_no_cds(workdir):
    """Features without a CDS would leave iVar reporting GFF_FEATURE=NA for everything."""
    gene_only = "\n".join(
        [
            "##gff-version 3",
            "##sequence-region S1_HA_polished 1 1775",
            "S1_HA_polished\tnextclade\tgene\t33\t1730\t.\t+\t.\tID=Gene-0-g;seq_index=0",
            "##sequence-region S2_HA_polished 1 1775",
            "S2_HA_polished\tnextclade\tgene\t33\t1730\t.\t+\t.\tID=Gene-1-g;seq_index=1",
        ]
    ) + "\n"
    workdir["combined"].write_text(gene_only, encoding="utf8")

    result = run_cli(workdir)

    assert result.returncode != 0
    assert "produced no CDS features" in result.stderr


def test_failure_quotes_nextclade_diagnostics(workdir):
    """The operator needs nextclade's own reason, not just our assertion."""
    workdir["combined"].write_text(combined_gff(queries=("S1_HA_polished",)), encoding="utf8")
    workdir["tsv"].write_text(
        nextclade_tsv(
            {"S1_HA_polished": {}, "S2_HA_polished": {"warnings": "low coverage", "failedCdses": "HA"}}
        ),
        encoding="utf8",
    )

    result = run_cli(workdir)

    assert "warnings: low coverage" in result.stderr
    assert "failedCdses: HA" in result.stderr


def test_fails_on_nextclade_error_column(workdir):
    """An unalignable query is absent from the GFF but still listed in the TSV."""
    workdir["tsv"].write_text(
        nextclade_tsv({"S1_HA_polished": {}, "S2_HA_polished": {"errors": "no alignment found"}}),
        encoding="utf8",
    )

    result = run_cli(workdir)

    assert result.returncode != 0
    assert "nextclade failed on 'S2_HA_polished': no alignment found" in result.stderr


def test_fails_on_unexpected_seqid(workdir):
    """Guards against the reference leaking into the annotation GFF, or mis-staged files."""
    leaked = combined_gff() + "HA\tnextclade\tCDS\t33\t1730\t.\t+\t0\tID=CDS-9-c;seq_index=9\n"
    workdir["combined"].write_text(leaked, encoding="utf8")

    result = run_cli(workdir)

    assert result.returncode != 0
    assert "unexpected non-consensus seqid(s)" in result.stderr
    assert "'HA'" in result.stderr


def test_warns_but_succeeds_when_a_reference_cds_is_dropped(workdir):
    """
    nextclade intersects each CDS with the aligned range and drops emptied features.
    That is a data-quality signal, not a pipeline error, so it warns and continues.
    """
    reference_two_cds = REFERENCE_ANNOTATION + (
        "HA\tRefSeq\tCDS\t100\t200\t.\t+\t0\tID=cds-PA-X;Parent=gene-FLUAVs4gp1;Name=PA-X\n"
    )
    workdir["reference"].write_text(reference_two_cds, encoding="utf8")

    result = run_cli(workdir)

    assert result.returncode == 0, result.stderr
    assert "missing 1 of 2 reference CDS features: cds-PA-X" in result.stderr
    assert (workdir["outdir"] / "S1_HA_polished.gff3").exists()


def test_single_query_run(workdir):
    """N=1 -- the SARS-CoV-2 path -- goes through the same code as the 8-segment fan-out."""
    workdir["combined"].write_text(combined_gff(queries=("S1_HA_polished",)), encoding="utf8")
    workdir["tsv"].write_text(nextclade_tsv({"S1_HA_polished": {}}), encoding="utf8")

    result = run_cli(workdir, fastas=[workdir["fastas"][0]])

    assert result.returncode == 0, result.stderr
    assert sorted(path.name for path in workdir["outdir"].iterdir()) == ["S1_HA_polished.gff3"]
