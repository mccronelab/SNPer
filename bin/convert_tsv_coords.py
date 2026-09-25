#!/usr/bin/env python3
# convert_tsv_coords.py
# Description: Accepts nextclade's per-segment alignment and a TSV file describing variants
# in iVar output format, and outputs a TSV where each variant carries the equivalent
# position on the reference genome. The goal here is to make it possible to directly compare the
# locations of variations across target genomes with possible indels.
# Author: Conner Copeland
# Contact: ccopelan@fredhutch.org
# Created 2025-03-12
# Updated 2026-08-27

import argparse
import csv
import sys
from typing import *

from Bio import SeqIO


GAP = "-"


def load_alignment_rows(aligned_path: str, reference_id: str, query_id: str) -> Tuple[str, str]:
    """
    Pulls the reference and one consensus row out of nextclade's per-segment alignment.

    Args:
        aligned_path (str): Path to nextclade's aligned FASTA. Holds every sample in the
            segment, so this stops reading once both rows are in hand.
        reference_id (str): Record ID of the reference, present because nextclade ran with
            --include-reference.
        query_id (str): Record ID of the consensus to convert.

    Returns:
        Tuple[str, str]: The reference row and the query row, both gapped to the alignment.
    """
    reference_row = None
    query_row = None

    for record in SeqIO.parse(aligned_path, "fasta"):
        if record.id == reference_id:
            reference_row = str(record.seq)
        elif record.id == query_id:
            query_row = str(record.seq)

        if reference_row is not None and query_row is not None:
            break

    if reference_row is None:
        sys.exit(
            f"reference record '{reference_id}' is missing from {aligned_path}; "
            "nextclade must run with --include-reference"
        )
    if query_row is None:
        sys.exit(f"consensus record '{query_id}' is missing from {aligned_path}")

    return reference_row, query_row


def read_insertions(nextclade_tsv: str, query_id: str) -> List[Tuple[int, str]]:
    """
    Reads one consensus's insertions out of nextclade's TSV.

    Nextclade strips insertions from the aligned FASTA and reports them here as
    `<pos>:<seq>`, anchored after 1-based reference position `pos`.

    Args:
        nextclade_tsv (str): Path to nextclade's per-segment TSV.
        query_id (str): Record ID of the consensus to look up.

    Returns:
        List[Tuple[int, str]]: (anchor position, inserted sequence), ascending.
    """
    with open(nextclade_tsv, newline="", encoding="utf8") as nextclade_file:
        for row in csv.DictReader(nextclade_file, delimiter="\t", quoting=csv.QUOTE_NONE):
            if row["seqName"] != query_id:
                continue

            insertions = []
            for item in (row.get("insertions") or "").split(","):
                item = item.strip()
                if not item:
                    continue
                position, sequence = item.split(":")
                insertions.append((int(position), sequence))

            return sorted(insertions)

    sys.exit(f"consensus record '{query_id}' is missing from {nextclade_tsv}")


def build_mapping(
    query_row: str, insertions: List[Tuple[int, str]]
) -> Tuple[Dict[int, int], Set[int], str]:
    """
    Maps consensus coordinates onto reference coordinates using the aligned row.

    Nextclade's alignment is reference-framed, so alignment column i is reference position
    i+1 and the reference coordinate falls out of the walk directly. Inserted bases have no
    reference position of their own and carry their anchor, matching the previous MAFFT
    behaviour; they are returned separately so the caller can report them.

    Args:
        query_row (str): The consensus's row of the alignment.
        insertions (List[Tuple[int, str]]): Insertions stripped from that row.

    Returns:
        Tuple[Dict[int,int], Set[int], str]: consensus->reference positions, the consensus
            positions that sit inside an insertion, and the consensus sequence rebuilt from
            the row plus its insertions.
    """
    # keyed by the reference position each insertion sits after, so the walk below can
    # splice it back into the consensus on reaching that position
    insertions_after = {}
    for position, sequence in insertions:
        insertions_after.setdefault(position, []).append(sequence)

    consensus_to_reference = {}
    inserted_positions = set()
    rebuilt = []
    consensus_position = 0

    # position 0 covers an insertion anchored before the first reference base
    for reference_position in range(0, len(query_row) + 1):
        if reference_position > 0:
            base = query_row[reference_position - 1]
            if base != GAP:
                consensus_position += 1
                consensus_to_reference[consensus_position] = reference_position
                rebuilt.append(base)

        for sequence in insertions_after.get(reference_position, []):
            for character in sequence:
                consensus_position += 1
                consensus_to_reference[consensus_position] = reference_position
                inserted_positions.add(consensus_position)
                rebuilt.append(character)

    return consensus_to_reference, inserted_positions, "".join(rebuilt)


def validate_mapping(
    aligned_path: str,
    reference_row: str,
    query_row: str,
    rebuilt: str,
    consensus: str,
    query_id: str,
) -> None:
    """
    Fails loudly unless the alignment supports the coordinates read off it.

    Two independent invariants. A gap-free reference row of the same length as the query is
    what makes column index equal reference position; rebuilding the consensus byte for byte
    is what makes the consensus numbering and the insertion parsing trustworthy. Neither
    covers the other, and a silent break in either would shift every REF_POS downstream.
    """
    if GAP in reference_row:
        sys.exit(
            f"Invalid MSA: the reference row of {aligned_path} contains gaps. Nextclade "
            "writes one alignment column per reference base, which is what lets a column "
            "number be read as a reference position. This alignment does not, so the "
            "variant positions in it cannot be converted."
        )

    if len(reference_row) != len(query_row):
        sys.exit(
            f"Invalid MSA: rows of {aligned_path} disagree in length -- reference "
            f"{len(reference_row)} bases, '{query_id}' {len(query_row)}. Nextclade should pad "
            "every row to the reference's length. These positions do not line up and "
            "cannot be converted."
        )

    if rebuilt.upper() != consensus.upper():
        detail = f"they are {len(rebuilt)} and {len(consensus)} bases long"
        for index, (left, right) in enumerate(zip(rebuilt.upper(), consensus.upper())):
            if left != right:
                detail = f"they first differ at consensus position {index + 1}"
                break
        sys.exit(
            f"'{query_id}' as rebuilt from {aligned_path} plus its insertions does not match "
            f"the consensus FASTA -- {detail}. The alignment does not describe this "
            "consensus, so its variant positions would be converted against the wrong "
            "sequence."
        )


def convert_coordinates(
    tsv_path: str,
    output_path: str,
    consensus_to_reference: Dict[int, int],
    inserted_positions: Set[int],
    query_id: str,
) -> None:
    """
    Reads in original TSV in target genome coordinates, then uses consensus_to_reference to
    look up equivalent coordinates on the reference genome. Writes to a new output TSV with
    positions in reference genome coordinates.

    Adds two columns. REF_POS is the reference coordinate. REF_POS_INSERTED marks the rows
    where that coordinate is inherited rather than the variant's own, because the base sits
    inside an insertion relative to the reference; those rows share a REF_POS with the base
    preceding the insertion, so without the flag the collision is invisible.

    Args:
        tsv_path (str): Tab-separated file describing variants on target genome. Expects the
            second field to be the location of the SNP on a chromosome.
        output_path (str): Path to output TSV file with coordinates on reference genome.
        consensus_to_reference (Dict[int,int]): Target genome coordinate keys and equivalent
            reference genome coordinate values, allowing us to map from target to reference.
        inserted_positions (Set[int]): Target coordinates with no reference position of their
            own, flagged per row and reported once at the end.
        query_id (str): Record ID of the consensus, used in messages.
    """
    variants = 0
    inside_insertions = []

    with open(tsv_path, "r", encoding="utf8") as variants_tsv:
        with open(output_path, "w", encoding="utf8") as reference_coords_tsv:
            # we expect there to be a header line in iVar TSV output
            header = variants_tsv.readline().strip()
            header = header + "\tREF_POS\tREF_POS_INSERTED\n"
            reference_coords_tsv.write(header)

            for line in variants_tsv:
                split_line = line.strip().split("\t")
                variants += 1

                # iVar's second field is POS on the consensus. Every consensus base is a key
                # in the mapping, so a miss means this TSV was called against a different
                # sequence than the one aligned.
                target_position = int(split_line[1])
                if target_position not in consensus_to_reference:
                    sys.exit(
                        f"variant at POS {target_position} lies past the end of consensus "
                        f"'{query_id}' ({len(consensus_to_reference)} bases)"
                    )
                is_inserted = target_position in inserted_positions
                if is_inserted:
                    inside_insertions.append(target_position)

                # get target coordinate and write; TRUE/FALSE matches iVar's PASS column
                split_line.append(str(consensus_to_reference[target_position]))
                split_line.append("TRUE" if is_inserted else "FALSE")

                updated_line = "\t".join(split_line) + "\n"
                reference_coords_tsv.write(updated_line)

    if inside_insertions:
        positions = ", ".join(str(position) for position in sorted(set(inside_insertions)))
        print(
            f"warning: {len(inside_insertions)} of {variants} variants in '{query_id}' fall "
            f"inside insertions relative to the reference; REF_POS carries the preceding "
            f"reference position (consensus {positions})",
            file=sys.stderr,
        )


def parse_args(sys_args: str) -> argparse.Namespace:
    """
    Parses command line arguments.

    Args:
        sys_args (str): Command line input with script name cut off.

    Returns:
        argparse.Namespace: Object containing argument values.
    """

    parser = argparse.ArgumentParser(
        sys_args,
        description = "Reads nextclade's per-segment alignment to create a mapping between " \
        "a consensus genome and the reference that accounts for indels. This mapping is then " \
        "used to convert the genomic coordinates of a TSV file describing variants on the " \
        "consensus genome to the equivalent coordinates in the reference genome. Outputs a " \
        "TSV file with converted coordinates that are relative to the reference genome."
    )

    parser.add_argument(
        "aligned_fasta",
        type = str,
        help = "Path to nextclade's alignment for this segment. Must include the " \
            "reference row, which nextclade emits under --include-reference."
    )

    parser.add_argument(
        "nextclade_tsv",
        type = str,
        help = "Path to nextclade's insertion TSV for this segment."
    )

    parser.add_argument(
        "consensus_fasta",
        type = str,
        help = "Path to consensus genome file. Must be in FASTA format and contain only one " \
            "sequence. Used to confirm the coordinates read off the alignment."
    )

    parser.add_argument(
        "target_variants",
        type = str,
        help = "Path to variant TSV file with coordinates relative to the target."
    )

    parser.add_argument(
        "output",
        type = str,
        help = "Path to output TSV file with an additional row, coordinates relative to the reference."
    )

    parser.add_argument(
        "--reference-id",
        type = str,
        required = True,
        help = "Record ID of the reference sequence in the alignment. This is the segment " \
            "name, which is the reference FASTA record ID by construction."
    )

    return parser.parse_args()


def _main():
    # input arguments
    args = parse_args(sys.argv[1:])

    consensus_record = SeqIO.read(args.consensus_fasta, "fasta")
    query_id = consensus_record.id

    reference_row, query_row = load_alignment_rows(
        args.aligned_fasta, args.reference_id, query_id
    )
    insertions = read_insertions(args.nextclade_tsv, query_id)

    consensus_to_reference, inserted_positions, rebuilt = build_mapping(query_row, insertions)
    validate_mapping(
        args.aligned_fasta, reference_row, query_row, rebuilt, str(consensus_record.seq), query_id
    )

    convert_coordinates(
        args.target_variants, args.output, consensus_to_reference, inserted_positions, query_id
    )


if __name__ == "__main__":
    _main()