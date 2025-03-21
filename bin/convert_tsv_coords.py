# convert_gff_coords.py
# Description: Accepts a reference genome, a target genome, and a TSV file describing variants
# in iVar output format. This script aligns the reference and target genome to create a genomic 
# coordinate mapping between them and outputs a TSV file where target variants occur at the
# equivalent position in the reference genome. The goal here is to make it possible to directly 
# compare the locations of variations across target genomes with possible indels.
# Author: Conner Copeland
# Contact: ccopelan@fredhutch.org
# Created 2025-03-12
# Updated 2025-03-13

import argparse
import subprocess
import sys

from Bio import SeqIO
from typing import *


def convert_coordinates(tsv_path: str, output_path: str, target_to_ref: Dict[int,int]) -> None:
    """
    Reads in original TSV in target genome coordinates, then uses a dict (target_to_ref) to look
    up equivalent coordinates on the reference genome. Writes to a new output TSV with positions
    in reference genome coordinates.

    Args:
        tsv_path (str): Tab-separated file describing variants on target genome. Expects the
            second field to be the location of the SNP on a chromosome.
        output_path (str): Path to output TSV file with coordinates on reference genome.
        target_to_ref (Dict[int,int]): Dictionary with target genome coordinate keys and equivalent
            reference genome coordinate values, allowing us to map from target to reference.
    """
    with open(tsv_path, "r", encoding="utf8") as variants_tsv:
        with open(output_path, "w", encoding="utf8") as reference_coords_tsv:
            # we expect there to be a header line in iVar TSV output
            header = variants_tsv.readline()
            reference_coords_tsv.write(header)

            for line in variants_tsv:
                # tsv are tab-separated
                split_line = line.split("\t")
                # get target coordinate, replace, and write
                split_line[1] = str(target_to_ref[int(split_line[1])])

                updated_line = "\t".join(split_line)
                reference_coords_tsv.write(updated_line)



def mafft_alignment(genomes_path: str) -> SeqIO.SeqRecord:
    """
    Aligns genomes in genomes_path with MAFFT.

    Args:
        genomes_path (str): Path to file containing reference and target genome.

    Returns:
        str: Uses subprocess to run MAFFT, which aligns the reference and target genomes.
    """
    temp_path = "alignment.msa"

    with open("alignment.msa", "w", encoding="utf8") as alignment_file:
        subprocess.run(["mafft", "--auto", f"{genomes_path}" ], stdout=alignment_file, check=True)

    # convert to list so we can access reference, target later   
    return list(SeqIO.parse(temp_path, format="fasta"))


def create_mapping(genomes_path: str) -> Dict[int, int]:
    """
    Takes in a path to a FASTA file containing reference and target genomes. Aligns the genomes
    and returns a mapping from the target genome coordinates to the reference genome
    coordinates.

    Args:
        genomes_path (str): Path to FASTA file with reference sequence, then target sequence.

    Returns:
        Dict: Dictionary where target coordinate keys map to reference coordinate values.
    """
    alignment_list = mafft_alignment(genomes_path)

    reference_ali = alignment_list[0].seq
    target_ali = alignment_list[1].seq

    # here, we want to create a mapping from the target genome to the reference genome in
    # nucleotide space
    target_to_ref_dict = {}

    ref_counter = 0
    target_counter = 0

    for i in range(len(reference_ali)):
        if reference_ali[i] != '-':
            ref_counter += 1

        if target_ali[i] != '-':
            target_counter += 1

        target_to_ref_dict[target_counter] = ref_counter

    return target_to_ref_dict


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
        description = "Aligns a target and reference genome, creating a mapping between the two " \
        "that accounts for indels. This mapping is then used to convert the genomic coordinates " \
        "of a TSV file describing variants on the target genome to the equivalent coordinates in " \
        "the reference genome. Outputs a TSV file with converted coordinates that are relative to " \
        "the reference genome."
    )

    parser.add_argument(
        "genomes_fasta",
        type = str,
        help = "Path to file containing both reference and target genomes. File is in FASTA format " \
        "and should list the reference sequence first."
    )

    parser.add_argument(
        "target_variants",
        type = str,
        help = "Path to variant TSV file with coordinates relative to the target."
    )

    parser.add_argument(
        "output",
        type = str,
        help = "Path to output TSV file with coordinates relative to the reference."
    )

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    genomes_path = args.genomes_fasta
    tsv_path = args.target_variants
    output_path = args.output

    target_to_ref_dict = create_mapping(genomes_path)
    convert_coordinates(tsv_path, output_path, target_to_ref_dict)


if __name__ == "__main__":
    _main()