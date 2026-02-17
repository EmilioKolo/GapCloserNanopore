#!/usr/bin/env python3
"""
extract_fasta_regions.py

Extract one or more genomic subsections from a FASTA file and write them
as separate sequences to a new FASTA file.

Each requested region is written as an independent FASTA record, and the
region string (chr:start-end) is included in the output sequence ID.

USAGE
-----
python extract_fasta_regions.py \
    --input genome.fasta \
    --region chr1:1000-2000 \
    --region chr2:500-900 \
    --output-folder results/ \
    --output-name genome_subsections.fasta

ARGUMENTS
---------
--input
    Input FASTA file.

--region
    Genomic region in the format chrN:start-end.
    Can be provided multiple times. Each region becomes a separate output
    FASTA sequence.

--output-folder
    Folder where the output FASTA file is created.
    If not provided, uses the same folder as the input FASTA.

--output-name
    Name of the output FASTA file.
    If ".fasta" is not present, it is appended automatically.
    If not provided, defaults to:
        <input_basename>_subsections.fasta
"""

import argparse
import os
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract subsections from a FASTA file."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input FASTA file."
    )
    parser.add_argument(
        "--region",
        action="append",
        required=True,
        help="Genomic region chr:start-end (can be used multiple times)."
    )
    parser.add_argument(
        "--output-folder",
        default=None,
        help="Output folder (default: same as input FASTA)."
    )
    parser.add_argument(
        "--output-name",
        default=None,
        help="Output FASTA file name."
    )

    return parser.parse_args()


def ensure_fasta_extension(name):
    """Ensure output file name ends with .fasta."""
    if not name.endswith(".fasta"):
        return name + ".fasta"
    return name


def parse_region(region_str):
    """
    Parse a region string of the form chr:start-end.

    Returns:
        chrom (str), start (int, 0-based), end (int, exclusive)
    """
    try:
        chrom, coords = region_str.split(":")
        start, end = coords.split("-")
        return chrom, int(start) - 1, int(end)
    except ValueError:
        raise ValueError(f"Invalid region format: {region_str}")


def main():
    args = parse_args()

    print("[INFO] Starting FASTA subsection extraction.")
    print(f"[INFO] Input FASTA: {args.input}")

    # Determine output folder
    if args.output_folder is None:
        output_folder = os.path.dirname(os.path.abspath(args.input))
    else:
        output_folder = args.output_folder

    os.makedirs(output_folder, exist_ok=True)

    # Determine output file name
    if args.output_name is None:
        base = os.path.splitext(os.path.basename(args.input))[0]
        output_name = f"{base}_subsections.fasta"
    else:
        output_name = ensure_fasta_extension(args.output_name)

    output_path = os.path.join(output_folder, output_name)
    print(f"[INFO] Output FASTA: {output_path}")

    print("[INFO] Loading FASTA sequences.")
    fasta_index = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

    output_records = []

    print(f"[INFO] Extracting {len(args.region)} region(s).")
    for region_str in args.region:
        chrom, start, end = parse_region(region_str)

        if chrom not in fasta_index:
            print(
                f"[ERROR] Chromosome '{chrom}' not found in FASTA.",
                file=sys.stderr
            )
            sys.exit(1)

        full_seq = fasta_index[chrom].seq
        subseq = full_seq[start:end]

        record_id = f"{chrom}:{start + 1}-{end}"
        record = SeqRecord(
            subseq,
            id=record_id,
            description=f"subsection_from_{chrom}"
        )

        output_records.append(record)

    print("[INFO] Writing output FASTA.")
    SeqIO.write(output_records, output_path, "fasta")

    print(f"[INFO] Finished. Regions written: {len(output_records)}")
    print("[INFO] Done.")


if __name__ == "__main__":
    main()