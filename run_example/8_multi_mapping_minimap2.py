#!/usr/bin/env python3

import argparse
import subprocess


def run_minimap2(index, fastq, sam_out, threads):
    """
    Run minimap2 allowing multiple alignments per read.
    """
    cmd = [
        "minimap2",
        "-a",              # SAM output
        "-x", "map-ont",   # preset for ONT long reads (change if needed)
        "-t", str(threads),
        "-N", "5",         # report up to 5 alignments per read
        index,
        fastq
    ]

    with open(sam_out, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


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
    parser = argparse.ArgumentParser(description="Extract multi-mapping reads using minimap2.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-x", "--index", required=True, help="Reference FASTA (minimap2 builds index on the fly)")
    parser.add_argument("--output-prefix", default="multimappers")
    parser.add_argument("-t", "--threads", type=int, default=4)

    args = parser.parse_args()

    sam_file = f"{args.output_prefix}.sam"
    out_fastq = f"{args.output_prefix}.multimappers.fastq"

    print("Running minimap2...")
    run_minimap2(args.index, args.input, sam_file, args.threads)

    print("Identifying multi-mapping reads...")
    multimappers = get_multimappers(sam_file)

    print(f"Found {len(multimappers)} multi-mapping reads.")

    print("Extracting reads...")
    extract_fastq(args.input, multimappers, out_fastq)

    print("Done.")
    print(f"Output FASTQ: {out_fastq}")


if __name__ == "__main__":
    main()