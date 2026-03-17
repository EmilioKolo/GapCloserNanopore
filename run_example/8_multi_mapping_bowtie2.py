#!/usr/bin/env python3

import argparse
import subprocess


def run_bowtie2(index, fastq, sam_out, threads):
    """
    Run Bowtie2 allowing multiple alignments per read.
    """
    cmd = [
        "bowtie2",
        "-x", index,
        "-U", fastq,
        "-S", sam_out,
        "-k", "2",  # detect multi-mappers
        "--very-sensitive",
        "-p", str(threads)
    ]
    subprocess.run(cmd, check=True)


def get_multimappers(sam_file):
    """
    Identify read names that have multiple alignments.
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
            if unmapped:
                continue

            counts[qname] = counts.get(qname, 0) + 1

    multimappers = {r for r, c in counts.items() if c > 1}
    return multimappers


def extract_fastq(input_fastq, read_ids, output_fastq):
    """
    Extract reads from FASTQ whose IDs are in read_ids.
    """
    with open(input_fastq) as fin, open(output_fastq, "w") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            read_id = header.split()[0][1:]

            if read_id in read_ids:
                fout.write(header)
                fout.write(seq)
                fout.write(plus)
                fout.write(qual)


def main():
    parser = argparse.ArgumentParser(description="Extract multi-mapping reads using Bowtie2.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-x", "--index", required=True, help="Bowtie2 index prefix")
    parser.add_argument("--output-prefix", default="multimappers")
    parser.add_argument("-t", "--threads", type=int, default=4)

    args = parser.parse_args()

    sam_file = f"{args.output_prefix}.sam"
    out_fastq = f"{args.output_prefix}.multimappers.fastq"

    print("Running Bowtie2...")
    run_bowtie2(args.index, args.input, sam_file, args.threads)

    print("Identifying multi-mapping reads...")
    multimappers = get_multimappers(sam_file)

    print(f"Found {len(multimappers)} multi-mapping reads.")

    print("Extracting reads...")
    extract_fastq(args.input, multimappers, out_fastq)

    print("Done.")
    print(f"Output FASTQ: {out_fastq}")


if __name__ == "__main__":
    main()