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
        "replicate_list",
        type=str,
        help="Path to a text file with one replicate name per line."
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

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    replicate_list_path = args.replicate_list
    fastq_dir_path = args.fastq_dir
    output_path = args.output_sheet

    replicate_list = []

    with open(replicate_list_path, "r", encoding="utf8") as replicate_file:
        for replicate in replicate_file:
            replicate_list.append(replicate.strip())

    with open(output_path, "w", encoding="utf8") as sampleSheet:
        sampleSheet.write(f"sample,replicate_id,fastq1,fastq2\n")

        for replicate in replicate_list:
            # currently, this will only work for underscores. If need be, use re.split() to support multiple break characters
            sample_id = replicate.split("_")[0]
            sampleSheet.write(f"{sample_id}, {replicate}, {fastq_dir_path}/{replicate}_1.fastq.gz, {fastq_dir_path}/{replicate}_2.fastq.gz\n")
    

if __name__ == "__main__":
    _main()