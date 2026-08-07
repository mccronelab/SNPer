#!/usr/bin/env python3

# generate_sample_sheet.py
# Description: This script accepts an input text file with one replicate per line and a path 
# to a directory containing FASTQ files. It then generates an output CSV format sample sheet.
# Replicates are assumed to have paired end reads.
# Author: Conner Copeland
# Contact: ccopelan@fredhutch.org
# Created: 2025-01-31
# Updated: 2025-07-29

import argparse
import sys
from typing import Union


def parse_args(sys_args: str) -> argparse.Namespace:
    """
    Parses arguments from command line.
    
    Args:
        -sys_args (str): Sys args with first entry (script name) cut off, containing arguments
            detailed on the command line by the user.
            
    Returns:
        -argparse.Namespace: Object containing variable defined below, associated with values fed
            into on the command line.
    """
    
    parser = argparse.ArgumentParser(
        sys_args,
        description="Generates a sample sheet based on a list of replicates and a path to a directory" \
            "containing FASTQ files. Output sample sheet is in CSV format."
    )

    parser.add_argument(
        "identifier_list",
        type=str,
        help="Path to a text file with one sample or replicate identifier per line. If a sample is \
            replicated, this should include a replicate identifier. For example, if sample 1234 \
            has replicates A and B, then this file should include 1234_A and 1234_B."
    )

    parser.add_argument(
        "fastq_dir",
        type=str,
        help="Path to directory containing FASTQ files"
    )

    parser.add_argument(
        "primer_id",
        type=str,
        default="",
        help="String value that will populate the primer_id for all rows. The value " \
        "should match a key in a primer CSV file supplied to the SNPer workflow."
    )

    parser.add_argument(
        "output_sheet",
        type=str,
        help="Path to output sample sheet in CSV format"
    )

    parser.add_argument(
        "--replicated",
        action="store_true",
        help="Set if input data has replicates. Default is False."
    )

    parser.add_argument(
        "--delimiter",
        type=str,
        default="_",
        help="Character to use to separate the sample ID from the replicate ID in input file names " \
             "(from identifier_list). For example, if a replicate is named 1234_A, this " \
             "should be '_'. If a replicate is named 1234-A, this should be '-'. Default is '_'."
    )

    parser.add_argument(
        "--protocol",
        type=str,
        default="",
        choices=["amplicon", "umi", "positional"],
        help="If set, creates an optional read_deduplication column in the output sample sheet. Any " \
        "value supplied to this option will populate all columns in the output sample sheet. Value " \
        "be one of amplicon, umi, or positional."
    )

    parser.add_argument(
        "--interleaved",
        action="store_true",
        help="Set if input reads are interleaved. Default is false."
    )

    parser.add_argument(
        "--suppress-header",
        action="store_true",
        help="Set to prevent the header line from being printed in the output sample sheet (useful "\
        "for concatenating multiple samplesheets together). Default is false."
    )

    parser.add_argument(
        "--r1-suffix",
        type=str,
        default="1",
        help="What to use as the suffix at the end of a file name to indicate r1. In 123_1.fastq.gz, " \
        "this would be the character 1 (another common form is R1, as in 123_R1.fastq.gz). Defaults " \
        "to 1."
    )

    parser.add_argument(
        "--r2-suffix",
        type=str,
        default="2",
        help="What to use as the suffix at the end of a file name to indicate r2. In 123_2.fastq.gz, " \
        "this would be the character 2 (another common form is R2, as in 123_R2.fastq.gz). Defaults " \
        "to 2."
    )

    return parser.parse_args()

def if_true_concat(condition: Union[str, bool], base_string: str, add_on: str) -> str:
    return_str = base_string

    if condition:
        return_str += add_on

    return return_str


def _main():
    args = parse_args(sys.argv[1:])
    identifier_list_path = args.identifier_list
    fastq_dir_path = args.fastq_dir
    output_path = args.output_sheet
    has_replicates = args.replicated
    delimiter = args.delimiter
    primer_id = args.primer_id
    read_deduplication = args.protocol
    interleaved = args.interleaved
    suppress_header = args.suppress_header
    r1_suffix = args.r1_suffix
    r2_suffix = args.r2_suffix

    file_ID_list = []

    with open(identifier_list_path, "r", encoding="utf8") as identifier_file:
        for file_ID in identifier_file:
            file_ID_list.append(file_ID.strip())

    with open(output_path, "w", encoding="utf8") as sampleSheet:
        if(not suppress_header):
            header = "sample,fastq1,fastq2,primer"

            # add optional columns, if specified
            header = if_true_concat(has_replicates, header, ",replicate-id")
            header = if_true_concat(read_deduplication, header, ",read_deduplication")
            header = if_true_concat(interleaved, header, ",interleaved")

            sampleSheet.write(f"{header}\n")

        for file_ID in file_ID_list:

            row = f"{file_ID}, {fastq_dir_path}/{file_ID}_{r1_suffix}.fastq.gz, {fastq_dir_path}/{file_ID}_{r2_suffix}.fastq.gz, {primer_id}"

            if has_replicates:
                # in this case, there's a replicate ID we want to split off from sample ID
                sample_list = file_ID.split(delimiter)[:-1]
                sample = delimiter.join(sample_list)
                row = f"{sample}, {fastq_dir_path}/{file_ID}_{r1_suffix}.fastq.gz, {fastq_dir_path}/{file_ID}_{r2_suffix}.fastq.gz, {primer_id}, {file_ID}"

            row = if_true_concat(read_deduplication, row, f", {read_deduplication}")
            row = if_true_concat(interleaved, row, ", True")

            sampleSheet.write(f"{row}\n")

if __name__ == "__main__":
    _main()