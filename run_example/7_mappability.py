#!/usr/bin/env python3

"""
Calculate the mappability of a DNA sequence by generating k-mers 
and aligning them to a reference genome using Bowtie2. 

The script counts how many k-mers map uniquely to the reference, 
providing an estimate of the mappability of the sequence.

IMPORTANT: You must have indexed the genome with Bowtie2 before 
running this script.
"""

import argparse
import subprocess

def read_fasta(path):
    """
    Read a FASTA file and return the concatenated sequence.
    Assumes a single sequence.
    """
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def generate_kmers(seq, k, out_fasta):
    """
    Generate all k-mers from a sequence and write them to a FASTA file.
    """
    with open(out_fasta, "w") as out:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            out.write(f">kmer_{i}\n{kmer}\n")


def run_bowtie2(index, kmers_fasta, sam_out):
    """
    Align kmers using Bowtie2.
    """
    cmd = [
        "bowtie2",
        "-x", index,
        "-f", kmers_fasta,
        "-S", sam_out,
        "--very-sensitive",
        "-k", "2"  # allow detection of multi-mappers
    ]

    subprocess.run(cmd, check=True)


def count_unique_mappings(sam_file):
    """
    Count unique alignments from a SAM file.
    A read is considered uniquely mapped if:
    - it is mapped
    - it appears only once
    """
    counts = {}

    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.split("\t")
            qname = fields[0]
            flag = int(fields[1])

            unmapped = flag & 4

            if not unmapped:
                counts[qname] = counts.get(qname, 0) + 1

    unique = sum(1 for v in counts.values() if v == 1)
    total = len(counts)

    return unique, total


def main():
    parser = argparse.ArgumentParser(description="Calculate mappability of a DNA sequence.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA sequence")
    parser.add_argument("-k", "--kmer", required=True, type=int, help="k-mer length")
    parser.add_argument("-x", "--index", required=True, help="Bowtie2 reference index prefix")
    parser.add_argument("-o", "--output", default="mappability.txt", help="Output file")

    args = parser.parse_args()

    print("Reading FASTA...")
    seq = read_fasta(args.input)

    kmers_fasta = "kmers.fa"
    sam_out = "kmers.sam"

    print("Generating k-mers...")
    generate_kmers(seq, args.kmer, kmers_fasta)

    print("Running Bowtie2 alignment...")
    run_bowtie2(args.index, kmers_fasta, sam_out)

    print("Counting unique mappings...")
    unique, total = count_unique_mappings(sam_out)

    mappability = unique / total if total else 0

    with open(args.output, "w") as out:
        out.write(f"Total_kmers\t{total}\n")
        out.write(f"Unique_mapped\t{unique}\n")
        out.write(f"Mappability\t{mappability:.6f}\n")

    print("Done.")
    print(f"Mappability: {mappability:.6f}")


if __name__ == "__main__":
    main()