#!/usr/bin/env python3
# calculate_coverage.py
# Description: We accept a fasta file with one sequence and write the proportion of Ns found in the file to stdout.
# Samples files with more than 1 sample will fail

import argparse
import subprocess
import sys

from Bio import SeqIO


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
        description = "Description: We accept a fasta file with one sequence and write the proportion of Ns found in the file to stdout.  " \
        "Samples files with more than 1 sample will fail"
    )

    parser.add_argument(
        "fasta",
        type = str,
        help = "Path to file containing the sequence in question. File must be in FASTA " \
            "format and contain only one sequence."
    )
    return parser.parse_args()

def get_coverage(record):
    length = len(record.seq)
    n = record.seq.upper().count("N")
    if(length==0):
        return 0
    else:
        return (length-n)/length

def _main():
    # input arguments
    args = parse_args(sys.argv[1:])
    fasta = args.fasta

    records = SeqIO.parse(fasta, "fasta")

    coverage = get_coverage(next(records))
    
    if next(records,None) is not None:
        ValueError(f'More than one entry found in ${fasta}')
    
    sys.stdout.write(str(coverage))


if __name__ == "__main__":
    _main()