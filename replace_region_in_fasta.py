#!/usr/bin/env python3
"""
replace_region_in_fasta.py

Replace a specific genomic interval in a reference FASTA sequence with an
arbitrary insert sequence taken from another FASTA file.

Key properties:
- Coordinates are 1-based and inclusive
- The replacement sequence may be longer or shorter than the original region
- FASTA line wrapping is normalized on output
- Designed for correctness and transparency, not speed

Typical use case:
- Gap filling
- Manual sequence correction
- Assembly polishing
"""

import argparse


def read_fasta(path):
    """
    Read a FASTA file into a dictionary mapping headers to sequences.

    Assumptions:
    - Each record has a unique header
    - Sequence lines may be wrapped arbitrarily
    """
    sequences = {}
    header = None
    chunks = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # Save the previous record
                if header is not None:
                    sequences[header] = "".join(chunks)

                # Start a new record
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

        # Save the final record
        if header is not None:
            sequences[header] = "".join(chunks)

    return sequences


def write_fasta(sequences, path, line_width=60):
    """
    Write a FASTA dictionary to disk using fixed line width.
    """
    with open(path, "w") as f:
        for header, seq in sequences.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


def replace_region(sequence, start, end, insert):
    """
    Replace a 1-based inclusive interval [start, end] in a sequence.

    Internally converts to 0-based indexing.
    """
    # convert to Python indexing
    start0 = start - 1
    return sequence[:start0] + insert + sequence[end:]


def main():
    parser = argparse.ArgumentParser(
        description="Replace a genomic region in a FASTA sequence with an arbitrary insert"
    )

    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA file"
    )
    parser.add_argument(
        "--insert",
        required=True,
        help="FASTA file containing the sequence to insert (must contain exactly one record)"
    )
    parser.add_argument(
        "--chrom",
        required=True,
        help="Chromosome / contig name to modify (e.g. ChrC)"
    )
    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="Start coordinate (1-based, inclusive)"
    )
    parser.add_argument(
        "--end",
        type=int,
        required=True,
        help="End coordinate (1-based, inclusive)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA file"
    )

    args = parser.parse_args()

    # Load FASTA files
    reference_seqs = read_fasta(args.reference)
    insert_seqs = read_fasta(args.insert)

    # Basic validation
    if args.chrom not in reference_seqs:
        raise ValueError(f"Chromosome '{args.chrom}' not found in reference FASTA")

    if len(insert_seqs) != 1:
        raise ValueError("Insert FASTA must contain exactly one sequence")

    if args.start < 1 or args.end < args.start:
        raise ValueError("Invalid coordinate range")

    # Extract sequences
    reference_seq = reference_seqs[args.chrom]
    insert_seq = next(iter(insert_seqs.values()))

    if args.end > len(reference_seq):
        raise ValueError("End coordinate exceeds chromosome length")

    # Report length change
    original_length = len(reference_seq)

    # Perform replacement
    reference_seqs[args.chrom] = replace_region(
        sequence=reference_seq,
        start=args.start,
        end=args.end,
        insert=insert_seq
    )

    new_length = len(reference_seqs[args.chrom])

    # Write output
    write_fasta(reference_seqs, args.output)

    # High-level logging (non-verbose)
    print("Replacement completed")
    print(f"Chromosome: {args.chrom}")
    print(f"Original length: {original_length}")
    print(f"New length: {new_length}")
    print(f"Length delta: {new_length - original_length}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
