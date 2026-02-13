#!/usr/bin/env python3
"""
extract_reads.py

Extract reads from one or more genomic regions in a BAM file and write them
to a compressed FASTQ file (.fastq.gz). Optionally include unmapped reads.

USAGE
-----
python extract_reads.py \
    --input sample.bam \
    --region chr1:1000-2000 \
    --region chr2:500-900 \
    --include-unmapped \
    --output-folder results/ \
    --output-name sample_regions.fastq.gz

ARGUMENTS
---------
--input
    Path to an input BAM file (must be readable; index recommended).

--region
    Genomic region in the format chr:start-end.
    Can be provided multiple times to specify multiple regions.

--include-unmapped
    Boolean flag. If set, unmapped reads are also written to the output FASTQ.

--output-folder
    Folder where the output FASTQ file will be created.
    If not provided, the output is written to the same folder as the input BAM.

--output-name
    Name of the output FASTQ file.
    If ".fastq.gz" is not present, it will be appended automatically.
    If not provided, the default is:
        <input_basename>_region.fastq.gz
"""

import argparse
from argparse import Namespace
import gzip
import os
import sys

import pysam


def parse_args() -> Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract reads from BAM regions and write FASTQ output."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input BAM file."
    )

    parser.add_argument(
        "--region",
        action="append",
        help=(
            "Genomic region (chr:start-end). "
            "Can be specified multiple times."
        )
    )

    parser.add_argument(
        "--include-unmapped",
        action="store_true",
        help="Include unmapped reads in the output."
    )

    parser.add_argument(
        "--output-folder",
        default=None,
        help="Folder for output fastq.gz (default: same as input BAM)."
    )

    parser.add_argument(
        "--output-name",
        default=None,
        help=(
            "Output FASTQ file name. '.fastq.gz' will be appended if missing."
        )
    )

    return parser.parse_args()


def ensure_fastq_gz(name:str) -> str:
    """Ensure the output file name ends with .fastq.gz."""
    if name.endswith(".fastq"):
        return name + ".gz"
    elif not name.endswith(".fastq.gz"):
        return name + ".fastq.gz"
    return name


def write_read_fastq(handle, read):
    """
    Write a single pysam.AlignedSegment to a FASTQ handle.

    If quality scores are missing, '*' is written as a placeholder.
    """
    seq = read.query_sequence or ""
    qual = read.qual if read.qual is not None else "*"

    handle.write(f"@{read.query_name}\n")
    handle.write(f"{seq}\n")
    handle.write("+\n")
    handle.write(f"{qual}\n")


def pipeline_read_extraction(
    input_bam: str,
    folder_out: str,
    name_out: str,
    selected_region,
    inc_unmapped: bool
):
    """
    Pipeline to perform extraction of the selected reads.
    """
    print("[INFO] Starting read extraction.")
    print(f"[INFO] Input BAM: {input_bam}")

    # Determine output folder
    if folder_out is None:
        output_folder = os.path.dirname(os.path.abspath(input_bam))
    else:
        output_folder = folder_out

    os.makedirs(output_folder, exist_ok=True)

    # Determine output file name
    if name_out is None:
        base = os.path.splitext(os.path.basename(input_bam))[0]
        output_name = f"{base}_region.fastq.gz"
    else:
        output_name = ensure_fastq_gz(name_out)

    output_path = os.path.join(output_folder, output_name)
    print(f"[INFO] Output FASTQ: {output_path}")

    if not selected_region and not inc_unmapped:
        print(
            "[ERROR] No regions specified and --include-unmapped not set. "
            "Nothing to extract.",
            file=sys.stderr
        )
        sys.exit(1)

    print("[INFO] Opening BAM file.")
    bam = pysam.AlignmentFile(input_bam, "rb")

    written_reads = set()
    total_written = 0

    print("[INFO] Writing fastq.gz output.")
    with gzip.open(output_path, "wt") as out_fq:

        # Extract reads from specified regions
        if selected_region:
            print(f"[INFO] Extracting reads from {len(selected_region)} region(s).")
            for region in selected_region:
                print(f"[INFO] Processing region: {region}")
                for read in bam.fetch(region=region):
                    if read.query_name in written_reads:
                        continue
                    write_read_fastq(out_fq, read)
                    written_reads.add(read.query_name)
                    total_written += 1

        # Optionally include unmapped reads
        if inc_unmapped:
            print("[INFO] Including unmapped reads.")
            for read in bam.fetch(until_eof=True):
                if not read.is_unmapped:
                    continue
                if read.query_name in written_reads:
                    continue
                write_read_fastq(out_fq, read)
                written_reads.add(read.query_name)
                total_written += 1

    bam.close()

    print(f"[INFO] Finished. Total reads written: {total_written}")
    print("[INFO] Done.")


def main():
    args = parse_args()

    pipeline_read_extraction(
        args.input,
        args.output_folder,
        args.output_name,
        args.region,
        args.include_unmapped
    )


if __name__ == "__main__":
    main()