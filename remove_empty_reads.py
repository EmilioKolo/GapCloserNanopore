#!/usr/bin/env python3
"""
remove_empty_reads.py

Remove empty reads from a FASTQ file.

A read is considered empty if:
  - The sequence length is zero, OR
  - The quality string is empty

Both uncompressed (.fastq) and compressed (.fastq.gz) files are supported.

Usage:
  python remove_empty_reads.py \
      --input input.fastq.gz \
      --output cleaned.fastq.gz
"""

import argparse
import gzip
import sys


def open_maybe_gzip(path, mode):
    """
    Open a file that may be gzip-compressed based on its extension.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def remove_empty_reads(input_fastq, output_fastq):
    """
    Stream through a FASTQ file and write only non-empty reads.

    Parameters
    input_fastq: str
        Path to input FASTQ (optionally gzipped).
    output_fastq: str
        Path to output FASTQ (optionally gzipped).
    """
    kept = 0
    removed = 0

    with open_maybe_gzip(input_fastq, "rt") as fin, \
         open_maybe_gzip(output_fastq, "wt") as fout:

        while True:
            header = fin.readline()
            if not header:
                break

            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            # Basic FASTQ integrity check
            if not (seq and plus and qual):
                raise ValueError("Malformed FASTQ file: incomplete record detected")

            # Strip newline characters for checks
            seq_stripped = seq.rstrip("\n")
            qual_stripped = qual.rstrip("\n")

            if len(seq_stripped) == 0 or len(qual_stripped) == 0:
                removed += 1
                continue

            fout.write(header)
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)
            kept += 1

    print(f"Reads kept   : {kept}", file=sys.stderr)
    print(f"Reads removed: {removed}", file=sys.stderr)


def main():
    """
    Parse command-line arguments and run the filtering.
    """
    parser = argparse.ArgumentParser(
        description="Remove empty reads from a FASTQ file"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input FASTQ file (.fastq or .fastq.gz)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTQ file (.fastq or .fastq.gz)"
    )

    args = parser.parse_args()
    remove_empty_reads(args.input, args.output)


if __name__ == "__main__":
    main()
