#!/usr/bin/env python3

# generate_sample_sheet.py
# Description: This script accepts an input text file with one replicate per line and a path 
# to a directory containing FASTQ files. It then generates an output CSV format sample sheet.
# Replicates are assumed to have paired end reads.
# Author: Conner Copeland
# Contact: ccopelan@fredhutch.org
# Created: 2025-01-31
# Updated: 2024-01-31

import argparse
import re
import sys
from typing import TextIO


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
        "output_sheet",
        type=str,
        help="Path to output sample sheet in CSV format"
    )

    parser.add_argument(
        "--replicated",
        action="store_true",
        help="Set if input data has replicates. Default is False."
    )

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    identifier_list_path = args.identifier_list
    fastq_dir_path = args.fastq_dir
    output_path = args.output_sheet
    has_replicates = args.replicated

    file_ID_list = []

    with open(identifier_list_path, "r", encoding="utf8") as identifier_file:
        for file_ID in identifier_file:
            file_ID_list.append(file_ID.strip())

    with open(output_path, "w", encoding="utf8") as sampleSheet:
        sampleSheet.write(f"sample,replicate_id,fastq1,fastq2\n")

        for file_ID in file_ID_list:
            if has_replicates:
                # in this case, there's a replicate ID we want to split off from sample ID
                # currently, this will only work for underscores. If need be, use re.split() to support multiple break characters
                sample_list = file_ID.split("_")[:-1]
                sample = "_".join(sample_list)
                sampleSheet.write(f"{sample}, {file_ID}, {fastq_dir_path}/{file_ID}_1.fastq.gz, {fastq_dir_path}/{file_ID}_2.fastq.gz\n")

            else:
                # in this case, sample ID and file ID are the same
                sampleSheet.write(f"{file_ID}, {file_ID}, {fastq_dir_path}/{file_ID}_1.fastq.gz, {fastq_dir_path}/{file_ID}_2.fastq.gz\n")
    

if __name__ == "__main__":
    _main()