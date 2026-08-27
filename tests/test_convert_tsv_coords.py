#!/usr/bin/env python3
# test_convert_tsv_coords.py
# Description: Unit tests for bin/convert_tsv_coords.py.
# Author: Conner Copeland
# Contact: ccopelan@fredhutch.org
# Created 2026-08-27

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "bin" / "convert_tsv_coords.py"


def _load_script():
    """Imports the bin/ script by path; it is a CLI, not an installed module."""
    spec = importlib.util.spec_from_file_location("convert_tsv_coords", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ctc = _load_script()


# ---- Fixture builders ----
#
# Nextclade's aligned FASTA is reference-framed: every row is the reference's length, the
# reference row is gap-free, and insertions are stripped out and reported in the TSV. These
# fixtures reproduce that layout rather than a general MSA.

REFERENCE = "ACGTACGTACGTACGTACGT"  # 20 bases
REFERENCE_ID = "SEG"
QUERY_ID = "sample_SEG_polished"

TSV_HEADER = "seqName\tclade\tinsertions\ttotalInsertions"


def write_alignment(tmp_path, query_row, *, reference_row=REFERENCE, name="aligned.fasta"):
    path = tmp_path / name
    path.write_text(f">{REFERENCE_ID}\n{reference_row}\n>{QUERY_ID}\n{query_row}\n")
    return path


def write_nextclade_tsv(tmp_path, insertions="", *, name="nc.tsv", query=QUERY_ID):
    path = tmp_path / name
    total = sum(len(item.split(":")[1]) for item in insertions.split(",") if item)
    path.write_text(f"{TSV_HEADER}\n{query}\tnone\t{insertions}\t{total}\n")
    return path


def write_consensus(tmp_path, sequence, *, name="consensus.fa", record_id=QUERY_ID):
    path = tmp_path / name
    path.write_text(f">{record_id}\n{sequence}\n")
    return path


def write_variants(tmp_path, positions, *, name="variants.tsv"):
    """Writes a minimal iVar-shaped TSV; only the second field (POS) is read."""
    path = tmp_path / name
    lines = ["REGION\tPOS\tREF\tALT"]
    for position in positions:
        lines.append(f"{QUERY_ID}\t{position}\tA\tG")
    path.write_text("\n".join(lines) + "\n")
    return path


def run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants, reference_id=REFERENCE_ID):
    output = tmp_path / "out.tsv"
    result = subprocess.run(
        [
            sys.executable, str(SCRIPT), str(aligned), str(nextclade_tsv),
            str(consensus), str(variants), str(output), "--reference-id", reference_id,
        ],
        capture_output=True,
        text=True,
    )
    return result, output


def read_ref_pos(output):
    """Returns [(POS, REF_POS)] from a converted TSV."""
    rows = output.read_text().strip().split("\n")
    assert rows[0].endswith("\tREF_POS")
    return [(int(r.split("\t")[1]), int(r.split("\t")[-1])) for r in rows[1:]]


# ---- build_mapping ----


def test_identity_mapping_when_consensus_matches_reference():
    mapping, inserted, rebuilt = ctc.build_mapping(REFERENCE, [])
    assert rebuilt == REFERENCE
    assert inserted == set()
    assert all(mapping[position] == position for position in range(1, len(REFERENCE) + 1))


def test_deletion_shifts_later_positions_and_leaves_the_preceding_base_alone():
    # consensus is missing reference bases 4-6
    query_row = REFERENCE[:3] + "---" + REFERENCE[6:]
    mapping, inserted, rebuilt = ctc.build_mapping(query_row, [])

    assert rebuilt == REFERENCE[:3] + REFERENCE[6:]
    assert inserted == set()
    # regression: the old MAFFT walk assigned on gap columns too, so consensus 3 reported 6
    assert mapping[3] == 3
    assert mapping[4] == 7


def test_insertion_carries_the_anchor_and_later_positions_resume():
    mapping, inserted, rebuilt = ctc.build_mapping(REFERENCE, [(4, "TTT")])

    assert rebuilt == REFERENCE[:4] + "TTT" + REFERENCE[4:]
    assert inserted == {5, 6, 7}
    assert mapping[4] == 4
    assert [mapping[position] for position in (5, 6, 7)] == [4, 4, 4]
    assert mapping[8] == 5


def test_insertion_before_the_first_reference_base_carries_zero():
    mapping, inserted, rebuilt = ctc.build_mapping(REFERENCE, [(0, "GG")])

    assert rebuilt == "GG" + REFERENCE
    assert inserted == {1, 2}
    assert mapping[1] == 0 and mapping[2] == 0
    assert mapping[3] == 1


def test_insertion_and_deletion_together():
    query_row = REFERENCE[:3] + "---" + REFERENCE[6:]
    mapping, inserted, rebuilt = ctc.build_mapping(query_row, [(10, "CC")])

    assert rebuilt == REFERENCE[:3] + REFERENCE[6:10] + "CC" + REFERENCE[10:]
    assert inserted == {8, 9}
    assert mapping[7] == 10
    assert mapping[8] == 10 and mapping[9] == 10
    assert mapping[10] == 11


# ---- read_insertions ----


def test_read_insertions_parses_and_sorts(tmp_path):
    path = write_nextclade_tsv(tmp_path, "22204:GAGCCAGAA,19194:NNN")
    assert ctc.read_insertions(str(path), QUERY_ID) == [(19194, "NNN"), (22204, "GAGCCAGAA")]


def test_read_insertions_handles_an_empty_field(tmp_path):
    path = write_nextclade_tsv(tmp_path, "")
    assert ctc.read_insertions(str(path), QUERY_ID) == []


def test_read_insertions_rejects_a_missing_query(tmp_path):
    path = write_nextclade_tsv(tmp_path, "", query="someone_else")
    with pytest.raises(SystemExit) as excinfo:
        ctc.read_insertions(str(path), QUERY_ID)
    assert QUERY_ID in str(excinfo.value)


# ---- load_alignment_rows ----


def test_load_alignment_rows_returns_both(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    reference_row, query_row = ctc.load_alignment_rows(str(aligned), REFERENCE_ID, QUERY_ID)
    assert reference_row == REFERENCE and query_row == REFERENCE


def test_load_alignment_rows_rejects_a_missing_reference(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    with pytest.raises(SystemExit) as excinfo:
        ctc.load_alignment_rows(str(aligned), "NOT_THERE", QUERY_ID)
    assert "include-reference" in str(excinfo.value)


def test_load_alignment_rows_rejects_a_missing_query(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    with pytest.raises(SystemExit) as excinfo:
        ctc.load_alignment_rows(str(aligned), REFERENCE_ID, "not_a_sample")
    assert "not_a_sample" in str(excinfo.value)


# ---- validate_mapping ----


def test_validate_rejects_a_gapped_reference_row():
    with pytest.raises(SystemExit) as excinfo:
        ctc.validate_mapping("aln.fa", "AC-GT", "ACGGT", "ACGGT", "ACGGT", QUERY_ID)
    assert "Invalid MSA" in str(excinfo.value)


def test_validate_rejects_rows_of_different_lengths():
    with pytest.raises(SystemExit) as excinfo:
        ctc.validate_mapping("aln.fa", "ACGT", "ACGTA", "ACGTA", "ACGTA", QUERY_ID)
    assert "Invalid MSA" in str(excinfo.value)


def test_validate_rejects_a_consensus_the_alignment_does_not_describe():
    with pytest.raises(SystemExit) as excinfo:
        ctc.validate_mapping("aln.fa", "ACGT", "ACGT", "ACGT", "ACGA", QUERY_ID)
    assert "consensus position 4" in str(excinfo.value)


def test_validate_accepts_a_case_mismatch():
    ctc.validate_mapping("aln.fa", "ACGT", "ACGT", "acgt", "ACGT", QUERY_ID)


# ---- End to end ----


def test_end_to_end_identity(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    consensus = write_consensus(tmp_path, REFERENCE)
    variants = write_variants(tmp_path, [1, 10, 20])

    result, output = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode == 0, result.stderr
    assert read_ref_pos(output) == [(1, 1), (10, 10), (20, 20)]
    assert result.stderr == ""


def test_end_to_end_with_a_deletion(tmp_path):
    query_row = REFERENCE[:3] + "---" + REFERENCE[6:]
    aligned = write_alignment(tmp_path, query_row)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    consensus = write_consensus(tmp_path, REFERENCE[:3] + REFERENCE[6:])
    variants = write_variants(tmp_path, [3, 4])

    result, output = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode == 0, result.stderr
    assert read_ref_pos(output) == [(3, 3), (4, 7)]


def test_end_to_end_variant_inside_an_insertion_warns(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path, "4:TTT")
    consensus = write_consensus(tmp_path, REFERENCE[:4] + "TTT" + REFERENCE[4:])
    variants = write_variants(tmp_path, [4, 5, 6, 8])

    result, output = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode == 0, result.stderr
    assert read_ref_pos(output) == [(4, 4), (5, 4), (6, 4), (8, 5)]
    assert "2 of 4 variants" in result.stderr
    assert "5, 6" in result.stderr


def test_end_to_end_no_warning_when_no_variant_is_inside_an_insertion(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path, "4:TTT")
    consensus = write_consensus(tmp_path, REFERENCE[:4] + "TTT" + REFERENCE[4:])
    variants = write_variants(tmp_path, [4, 8])

    result, _ = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode == 0, result.stderr
    assert result.stderr == ""


def test_end_to_end_preserves_every_original_column(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    consensus = write_consensus(tmp_path, REFERENCE)
    variants = write_variants(tmp_path, [7])

    _, output = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    header, row = output.read_text().strip().split("\n")
    assert header.split("\t") == ["REGION", "POS", "REF", "ALT", "REF_POS"]
    assert row.split("\t") == [QUERY_ID, "7", "A", "G", "7"]


def test_end_to_end_empty_variant_tsv_still_writes_a_header(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    consensus = write_consensus(tmp_path, REFERENCE)
    variants = write_variants(tmp_path, [])

    result, output = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode == 0, result.stderr
    assert output.read_text() == "REGION\tPOS\tREF\tALT\tREF_POS\n"


def test_end_to_end_rejects_a_position_past_the_consensus(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    consensus = write_consensus(tmp_path, REFERENCE)
    variants = write_variants(tmp_path, [21])

    result, _ = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode != 0
    assert "past the end" in result.stderr


def test_end_to_end_rejects_a_consensus_that_disagrees_with_the_alignment(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    # the insertion is real but absent from the TSV, so the rebuild comes up short
    consensus = write_consensus(tmp_path, REFERENCE[:4] + "TTT" + REFERENCE[4:])
    variants = write_variants(tmp_path, [1])

    result, _ = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode != 0
    assert "does not match" in result.stderr


def test_end_to_end_rejects_a_multi_record_consensus(tmp_path):
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    consensus = tmp_path / "two.fa"
    consensus.write_text(f">{QUERY_ID}\n{REFERENCE}\n>extra\n{REFERENCE}\n")
    variants = write_variants(tmp_path, [1])

    result, _ = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode != 0


def test_end_to_end_reads_the_record_id_from_the_file_not_the_filename(tmp_path):
    # the filename deliberately disagrees with the record ID inside it
    aligned = write_alignment(tmp_path, REFERENCE)
    nextclade_tsv = write_nextclade_tsv(tmp_path)
    consensus = write_consensus(tmp_path, REFERENCE, name="unrelated_name.fa")
    variants = write_variants(tmp_path, [2])

    result, output = run_cli(tmp_path, aligned, nextclade_tsv, consensus, variants)

    assert result.returncode == 0, result.stderr
    assert read_ref_pos(output) == [(2, 2)]
