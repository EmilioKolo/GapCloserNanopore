#!/usr/bin/env python3

import argparse
import gzip
import subprocess


def run_minimap2(index, fastq, sam_out, threads):
    """
    Run minimap2 allowing split alignments.
    """
    cmd = [
        "minimap2",
        "-a",
        "-x", "map-ont",
        "-t", str(threads),
        "--secondary=yes",
        "-N", "10",
        index,
        fastq
    ]

    with open(sam_out, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def get_split_reads(sam_file):
    """
    Identify reads with supplementary alignments (split reads).
    """
    split_reads = set()

    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.split("\t")
            qname = fields[0]
            flag = int(fields[1])

            supplementary = flag & 0x800

            if supplementary:
                split_reads.add(qname)

    return split_reads


def main():
    parser = argparse.ArgumentParser(
        description="Count split reads using minimap2."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-x", "--index", required=True, 
                        help="Reference FASTA or .mmi")
    parser.add_argument("--output-prefix", default="split_reads")
    parser.add_argument("-t", "--threads", type=int, default=4)

    args = parser.parse_args()

    sam_file = f"{args.output_prefix}.sam"
    out_fastq = f"{args.output_prefix}.split_reads.fastq"

    print("Running minimap2...")
    run_minimap2(args.index, args.input, sam_file, args.threads)

    print("Identifying split reads...")
    split_reads = get_split_reads(sam_file)

    print(f"Found {len(split_reads)} split reads.")

    print("Extracting reads to FASTQ...")
    # Detect gzip-compressed input
    if args.input.endswith(".gz"):
        fin_handle = gzip.open(args.input, "rt")
    else:
        fin_handle = open(args.input, "r")

    with fin_handle as fin, open(out_fastq, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            read_id = header.split()[0][1:]

            if read_id in split_reads:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)

    print("Done.")
    print(f"Output FASTQ: {out_fastq}")


if __name__ == "__main__":
    main()