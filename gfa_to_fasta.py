#!/usr/bin/env python3
"""
gfa_to_fasta.py

Extract sequences from a GFA file and write them to a FASTA file.
This script reads segment (S) lines from the GFA and converts each
segment into a FASTA record.

USAGE
-----
python gfa_to_fasta.py \
    --input assembly.gfa \
    --output-name assembly.fasta \
    [--output-folder results/]

ARGUMENTS
---------
--input
    Input GFA file.

--output-name
    Name of the output FASTA file.
    If ".fasta" is not present, it is appended automatically.

--output-folder
    Folder where the output FASTA file is created.
    If not provided, uses the same folder as the input GFA file.
"""

import argparse
import os
import sys


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert GFA segments to FASTA."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input GFA file."
    )
    parser.add_argument(
        "--output-name",
        required=True,
        help="Output FASTA file name."
    )
    parser.add_argument(
        "--output-folder",
        default=None,
        help="Output folder (default: same as input GFA)."
    )

    return parser.parse_args()


def ensure_fasta_extension(name):
    """Ensure output file name ends with .fasta."""
    if not name.endswith(".fasta"):
        return name + ".fasta"
    return name


def main():
    args = parse_args()

    print("[INFO] Starting GFA to FASTA conversion.")
    print(f"[INFO] Input GFA: {args.input}")

    # Determine output folder
    if args.output_folder is None:
        output_folder = os.path.dirname(os.path.abspath(args.input))
    else:
        output_folder = args.output_folder

    os.makedirs(output_folder, exist_ok=True)

    # Determine output file path
    output_name = ensure_fasta_extension(args.output_name)
    output_path = os.path.join(output_folder, output_name)

    print(f"[INFO] Output FASTA: {output_path}")

    sequences_written = 0

    print("[INFO] Reading GFA and extracting sequences.")
    with open(args.input, "r") as gfa, open(output_path, "w") as fasta:
        for line in gfa:
            if not line.startswith("S\t"):
                continue

            fields = line.rstrip().split("\t")

            # GFA S-line format: S <name> <sequence> [optional tags]
            if len(fields) < 3:
                continue

            seq_id = fields[1]
            sequence = fields[2]

            if sequence == "*" or not sequence:
                continue

            fasta.write(f">{seq_id}\n")
            fasta.write(f"{sequence}\n")
            sequences_written += 1

    if sequences_written == 0:
        print(
            "[WARNING] No sequences were written. "
            "The GFA may not contain segment sequences.",
            file=sys.stderr
        )
    else:
        print(f"[INFO] Sequences written: {sequences_written}")

    print("[INFO] Done.")


if __name__ == "__main__":
    main()
