#!/usr/bin/env python3
"""
align_fastas.py

Run local BLAST and MAFFT on two FASTA files and generate an alignment.

USAGE
-----
python align_fastas.py \
    --fasta1 seq1.fasta \
    --fasta2 seq2.fasta \
    --out-prefix alignment_results \
    [--blast-type blastn] \
    [--mafft-threads 4]

REQUIRED ARGUMENTS
------------------
--fasta1
    First input FASTA file.

--fasta2
    Second input FASTA file.

--out-prefix
    Prefix for all output files.

OPTIONAL ARGUMENTS
------------------
--blast-type
    BLAST program to use (default: blastn).

--mafft-threads
    Number of threads to pass to MAFFT (default: 1).

OUTPUT FILES
------------
<out-prefix>.blast.tsv
    Tabular BLAST output (outfmt 6).

<out-prefix>.mafft.fasta
    MAFFT alignment of both sequences.
"""

import argparse
import subprocess
from pathlib import Path


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Align two FASTA files using BLAST and MAFFT."
    )

    parser.add_argument(
        "--fasta1",
        required=True,
        help="First FASTA file."
    )
    parser.add_argument(
        "--fasta2",
        required=True,
        help="Second FASTA file."
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output file prefix."
    )
    parser.add_argument(
        "--blast-type",
        default="blastn",
        help="BLAST program to use (default: blastn)."
    )
    parser.add_argument(
        "--mafft-threads",
        type=int,
        default=1,
        help="Number of threads for MAFFT (default: 1)."
    )

    return parser.parse_args()

def run_blast(fasta1, fasta2, blast_type, out_tsv):
    """
    Run local BLAST between two FASTA files.

    The second FASTA is used as a temporary BLAST database.
    """
    print("[INFO] Running local BLAST.")

    db_prefix = Path(out_tsv).with_suffix("")

    subprocess.run(
        [
            "makeblastdb",
            "-in", fasta2,
            "-dbtype", "nucl",
            "-out", str(db_prefix),
        ],
        check=True
    )

    subprocess.run(
        [
            blast_type,
            "-query", fasta1,
            "-db", str(db_prefix),
            "-outfmt", "6",
            "-out", out_tsv,
        ],
        check=True
    )

def run_mafft(fasta1, fasta2, threads, out_fasta):
    """
    Run MAFFT on two FASTA files and write an alignment.
    """
    print("[INFO] Running MAFFT alignment.")

    combined_fasta = out_fasta + ".tmp.fasta"

    with open(combined_fasta, "w") as out:
        with open(fasta1) as f1, open(fasta2) as f2:
            out.write(f1.read())
            out.write(f2.read())

    subprocess.run(
        [
            "mafft",
            "--thread", str(threads),
            combined_fasta,
        ],
        stdout=open(out_fasta, "w"),
        stderr=subprocess.DEVNULL,
        check=True
    )

def main():
    args = parse_args()

    fasta1 = args.fasta1
    fasta2 = args.fasta2
    prefix = args.out_prefix

    blast_out = f"{prefix}.blast.tsv"
    mafft_out = f"{prefix}.mafft.fasta"

    print("[INFO] Starting two-FASTA alignment workflow.")
    print(f"[INFO] FASTA1: {fasta1}")
    print(f"[INFO] FASTA2: {fasta2}")

    run_blast(fasta1, fasta2, args.blast_type, blast_out)
    run_mafft(fasta1, fasta2, args.mafft_threads, mafft_out)

    print("[INFO] All steps completed successfully.")
    print(f"[INFO] BLAST output: {blast_out}")
    print(f"[INFO] MAFFT alignment: {mafft_out}")


if __name__ == "__main__":
    main()
