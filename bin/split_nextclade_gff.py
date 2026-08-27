#!/usr/bin/env python3
# split_nextclade_gff.py
# Description: Splits the combined GFF3 that `nextclade run --output-annotation-gff`
# writes for a batch of queries into one file per consensus sequence, and normalizes
# the attribute column so the output is byte-stable across runs.
# Author: Conner Copeland
# Contact: ccopelan@fredhutch.org
# Created 2026-08-24

import argparse
import csv
import os
import re
import sys

from Bio import SeqIO
from typing import *


# nextclade rewrites feature IDs to `CDS-<seq_index>-<original>` / `Gene-<seq_index>-<original>`,
# where <seq_index> is the query's position within the run. Left in place, adding one sample
# would change every other sample's iVar GFF_FEATURE values.
ID_PREFIX = re.compile(r"^(?:CDS|Gene)-\d+-")

# nextclade serializes attributes out of a hash map with a per-process seed, so the order
# permutes on every run. Emit a fixed order instead: these first, then the rest sorted.
ATTRIBUTE_ORDER = ["ID", "Parent", "Name"]

# carries the query's index within the run, which is meaningless once the file is split
DROPPED_ATTRIBUTES = {"seq_index"}


def read_query_ids(consensus_paths: List[str]) -> Dict[str, str]:
    """
    Maps each consensus FASTA's record ID to its path.

    The record ID is read from the file rather than parsed from the filename:
    `<sample>_<segment>_polished` cannot be split on "_" unambiguously (segment `N_A`,
    sample `HS10498_A`). merge_mpileup_consensus.nf passes the same string to `ivar
    consensus -p` and `-i`, so the record ID and the file's basename agree by construction.
    """
    query_ids = {}

    for path in consensus_paths:
        # SeqIO.read raises if the file does not hold exactly one record
        record = SeqIO.read(path, "fasta")

        if record.id in query_ids:
            sys.exit(f"duplicate consensus record ID '{record.id}' ({query_ids[record.id]}, {path})")

        query_ids[record.id] = path

    return query_ids


def read_nextclade_tsv(tsv_path: str) -> Dict[str, Dict[str, str]]:
    """
    Reads nextclade's per-query TSV into a dict keyed on seqName.

    QUOTE_NONE because nextclade's free-text columns (errors, warnings) are not quoted
    and may contain quote characters, which would otherwise be read as field quoting.
    """
    with open(tsv_path, "r", encoding="utf8") as tsv_file:
        rows = {}

        for row in csv.DictReader(tsv_file, delimiter="\t", quoting=csv.QUOTE_NONE):
            rows[row.get("seqName", "")] = row

    return rows


def normalize_attributes(attributes: str) -> str:
    """
    Rewrites a GFF3 column 9 into a deterministic form: drops nextclade's bookkeeping
    attributes, strips its per-run ID prefixes, and emits the rest in a fixed order.
    Keys and values are otherwise untouched, so iVar's GFF_FEATURE is unchanged.
    """
    parsed = {}

    for pair in attributes.split(";"):
        if not pair:
            continue

        key, _, value = pair.partition("=")

        if key in DROPPED_ATTRIBUTES:
            continue

        # the prefix appears on ID and on the Parent that references it
        if key in ("ID", "Parent"):
            value = ID_PREFIX.sub("", value)

        parsed[key] = value

    ordered = [key for key in ATTRIBUTE_ORDER if key in parsed]
    ordered += sorted(key for key in parsed if key not in ATTRIBUTE_ORDER)

    return ";".join(f"{key}={parsed[key]}" for key in ordered)


def parse_gff(gff_path: str) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """
    Groups a combined GFF3's data lines by seqid, normalizing each line's attributes.

    Returns:
        Tuple of (seqid -> normalized data lines, seqid -> ##sequence-region pragma).
    """
    features = {}
    regions = {}

    with open(gff_path, "r", encoding="utf8") as gff_file:
        for line in gff_file:
            line = line.rstrip("\n")

            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) > 1:
                    regions[fields[1]] = line
                continue

            if line.startswith("#") or not line.strip():
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                sys.exit(f"malformed GFF3 line in {gff_path}: {line!r}")

            fields[8] = normalize_attributes(fields[8])
            features.setdefault(fields[0], []).append("\t".join(fields))

    return features, regions


def cds_ids(lines: List[str]) -> Set[str]:
    """Collects the distinct CDS ID values among a seqid's feature lines."""
    ids = set()

    for line in lines:
        fields = line.split("\t")

        if fields[2] != "CDS":
            continue

        for pair in fields[8].split(";"):
            key, _, value = pair.partition("=")
            if key == "ID":
                ids.add(value)

    return ids


def fail_for_query(query_id: str, reason: str, tsv_rows: Dict[str, Dict[str, str]]) -> NoReturn:
    """Aborts, quoting nextclade's own diagnostics for the offending query."""
    row = tsv_rows.get(query_id, {})
    details = [f"{query_id}: {reason}"]

    for column in ("errors", "warnings", "failedCdses"):
        value = row.get(column, "")
        if value:
            details.append(f"  {column}: {value}")

    sys.exit("\n".join(details))


def write_split_gffs(
    outdir: str,
    query_ids: Dict[str, str],
    features: Dict[str, List[str]],
    regions: Dict[str, str],
    tsv_rows: Dict[str, Dict[str, str]],
    reference_cds: Set[str],
) -> None:
    """Writes one GFF3 per consensus record, failing if any query lost its annotation."""
    for query_id in sorted(query_ids):
        lines = features.get(query_id, [])

        if not lines:
            fail_for_query(query_id, "nextclade produced no annotation features", tsv_rows)

        found_cds = cds_ids(lines)
        if not found_cds:
            fail_for_query(query_id, "nextclade produced no CDS features", tsv_rows)

        # nextclade intersects each CDS with the alignment range and drops emptied
        # features, so poor terminal coverage silently costs a CDS. iVar would then
        # report GFF_FEATURE=NA for variants inside it rather than failing.
        missing = reference_cds - found_cds
        if missing:
            print(
                f"warning: {query_id} is missing {len(missing)} of {len(reference_cds)} "
                f"reference CDS features: {', '.join(sorted(missing))}",
                file=sys.stderr,
            )

        out_path = os.path.join(outdir, f"{query_id}.gff3")

        with open(out_path, "w", encoding="utf8") as out_file:
            out_file.write("##gff-version 3\n")

            if query_id in regions:
                out_file.write(regions[query_id] + "\n")

            for line in lines:
                out_file.write(line + "\n")


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Splits nextclade's combined query annotation GFF3 into one file per "
        "consensus sequence, normalizing the attribute column so output is stable across runs."
    )

    parser.add_argument(
        "annotation_gff",
        type=str,
        help="Combined GFF3 from `nextclade run --output-annotation-gff`.",
    )

    parser.add_argument(
        "nextclade_tsv",
        type=str,
        help="Per-query TSV from `nextclade run --output-tsv`, used for error reporting.",
    )

    parser.add_argument(
        "reference_annotation",
        type=str,
        help="The single-segment reference annotation passed to `nextclade run "
        "--input-annotation`, used to report CDS features nextclade dropped.",
    )

    parser.add_argument(
        "outdir",
        type=str,
        help="Directory to write per-consensus GFF3 files into.",
    )

    parser.add_argument(
        "consensus_fasta",
        type=str,
        nargs="+",
        help="Consensus FASTA files that were passed to nextclade as queries. Each must "
        "contain exactly one record; its record ID names the output file.",
    )

    return parser.parse_args(sys_args)


def _main():
    args = parse_args(sys.argv[1:])

    query_ids = read_query_ids(args.consensus_fasta)
    tsv_rows = read_nextclade_tsv(args.nextclade_tsv)

    # a query nextclade failed to align is absent from the alignment but still listed
    # here; without this it would silently vanish from the run instead of stopping it
    for seq_name, row in sorted(tsv_rows.items()):
        if row.get("errors"):
            sys.exit(f"nextclade failed on '{seq_name}': {row['errors']}")

    features, regions = parse_gff(args.annotation_gff)

    # --include-reference keeps the reference out of the annotation GFF, so every seqid
    # here should name a consensus; anything else means the wrong files were staged
    unexpected = set(features) - set(query_ids)
    if unexpected:
        sys.exit(
            f"unexpected non-consensus seqid(s) in {args.annotation_gff}: {sorted(unexpected)}"
        )

    # the expected CDS set comes from the annotation we asked nextclade for, not from
    # its own output, so a feature it dropped entirely is still counted as missing
    reference_features, _ = parse_gff(args.reference_annotation)
    reference_cds = set()
    for lines in reference_features.values():
        reference_cds |= cds_ids(lines)

    write_split_gffs(args.outdir, query_ids, features, regions, tsv_rows, reference_cds)


if __name__ == "__main__":
    _main()
