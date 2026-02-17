#!/usr/bin/env python3
"""
plot_alignments.py

Generate alignment plots from BLAST and/or MAFFT outputs produced by
align_fastas.py.

If either input file is missing, the corresponding plot is skipped and
the script proceeds with the available data.

USAGE
-----
python plot_alignments.py \
    --blast-tsv alignment.blast.tsv \
    --mafft-fasta alignment.mafft.fasta \
    --out-prefix alignment_plots \
    [--anchor-id-substring Chr] \
    [--fasta-contigs contigs.fasta]

ARGUMENTS
---------
--blast-tsv
    BLAST tabular output (outfmt 6).

--mafft-fasta
    MAFFT alignment FASTA file.

--out-prefix
    Prefix for output plot files.

--anchor-id-substring
    Substring used to identify anchor sequences in the MAFFT alignment.

--fasta-contigs
    FASTA file with ungapped query contigs (used if MAFFT is not provided).

OUTPUT FILES
------------
<out-prefix>.blast.png
    BLAST alignment plot (query vs subject coordinates).

<out-prefix>.mafft_identity.png
    Percent identity across MAFFT alignment.
"""

import argparse
import os

import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO


def get_query_lengths_dict(mafft_fasta=None, fasta_contigs=None):
    """
    Build a {query_id: ungapped_length} dictionary.

    Priority:
      1) MAFFT alignment (gaps removed)
      2) FASTA contigs (raw length)
    """

    lengths = {}

    if mafft_fasta:
        print("[INFO] Extracting query lengths from MAFFT alignment.")
        alignment = AlignIO.read(mafft_fasta, "fasta")
        for rec in alignment:
            ungapped_len = len(str(rec.seq).replace("-", ""))
            lengths[rec.id] = ungapped_len
        return lengths

    if fasta_contigs:
        print("[INFO] Extracting query lengths from FASTA contigs.")
        for rec in SeqIO.parse(fasta_contigs, "fasta"):
            lengths[rec.id] = len(rec.seq)
        return lengths

    print(
        "[WARNING] No MAFFT or FASTA contigs provided. "
        "Cannot determine query lengths."
    )
    return None


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot BLAST and MAFFT alignment results."
    )

    parser.add_argument(
        "--blast-tsv",
        default=None,
        help="BLAST output file (outfmt 6)."
    )

    parser.add_argument(
        "--mafft-fasta",
        default=None,
        help="MAFFT alignment FASTA file."
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Prefix for output plot files."
    )
    parser.add_argument(
        "--anchor-id-substring",
        default='Chr',
        help="Substring used to identify anchor sequences in MAFFT alignment."
    )
    parser.add_argument(
        "--fasta-contigs",
        default=None,
        help="FASTA file with ungapped query contigs (used if MAFFT is not provided)."
    )

    return parser.parse_args()


def plot_blast(blast_tsv, out_png, query_lengths_dict):
    """
    Generate BLAST alignment plots showing the top 5 hits per subject.

    X-axis represents the full subject coordinate space.

    Subject lengths must be provided explicitly as a dict:
        {query_id: length}
    """

    print("[INFO] Generating BLAST top-hit plots.")

    if query_lengths_dict is None:
        print(
            "[WARNING] Subject lengths not provided. "
            "Cannot plot BLAST alignments meaningfully. Skipping."
        )
        return

    by_subject = {}

    with open(blast_tsv) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue

            qid = fields[0]
            sid = fields[1]
            s_start = int(fields[8])
            s_end = int(fields[9])
            bitscore = float(fields[11])

            s_lo = min(s_start, s_end)
            s_hi = max(s_start, s_end)

            by_subject.setdefault(sid, []).append(
                (qid, s_lo, s_hi, bitscore)
            )

    if not by_subject:
        print("[WARNING] BLAST file contains no alignments. Skipping plot.")
        return

    for sid, hits in by_subject.items():
        if sid not in query_lengths_dict:
            print(
                f"[WARNING] No length found for subject {sid}. Skipping."
            )
            continue

        # Group hits by anchor (query ID)
        by_anchor = {}
        for qid, s_lo, s_hi, bitscore in hits:
            by_anchor.setdefault(qid, []).append((s_lo, s_hi, bitscore))

        # Keep top 5 hits PER anchor
        top_hits = []
        for qid, ahits in by_anchor.items():
            ahits.sort(key=lambda x: x[2], reverse=True)
            for s_lo, s_hi, bitscore in ahits[:5]:
                top_hits.append((qid, s_lo, s_hi, bitscore))

        subject_len = query_lengths_dict[sid]

        plt.figure(figsize=(12, 3))

        # Subject baseline
        plt.hlines(
            y=0,
            xmin=1,
            xmax=subject_len,
            linewidth=8,
            color="black"
        )

        # Anchor hits
        y_labels = []
        for i, (qid, s_lo, s_hi, bitscore) in enumerate(top_hits, start=1):
            plt.hlines(
                y=i,
                xmin=s_lo,
                xmax=s_hi,
                linewidth=8
            )
            y_labels.append(qid)

        plt.yticks(
            [0] + list(range(1, len(top_hits) + 1)),
            ["subject"] + y_labels
        )

        padding = int(0.02 * subject_len)

        plt.xlim(1 - padding, subject_len + padding)
        plt.ylim(-1, len(top_hits) + 1)
        plt.xlabel("Subject coordinate")
        plt.ylabel("Alignments")
        plt.title(f"Top 5 BLAST hits on {sid}")
        plt.tight_layout()

        out_file = f"{out_png}.{sid}.png"
        plt.savefig(out_file)
        plt.close()

        print(f"[INFO] Saved BLAST plot: {out_file}")
    return None


def plot_mafft_identity(
    mafft_fasta,
    out_png,
    anchor_id_substring,
    window=50,
    identity_percent=50.0
):
    """
    For each non-anchor sequence, generate a plot showing where each 
    anchor aligns with it at high identity.

    Anchors are sequences whose IDs contain 'anchor_id_substring'.
    One plot is produced per non-anchor sequence.
    """

    print("[INFO] Generating MAFFT anchor vs non-anchor plots.")

    alignment = AlignIO.read(mafft_fasta, "fasta")

    anchors = [
        rec for rec in alignment
        if anchor_id_substring in rec.id
    ]
    others = [
        rec for rec in alignment
        if anchor_id_substring not in rec.id
    ]

    if len(anchors) != 2:
        print(
            f"[WARNING] Expected 2 anchors, found {len(anchors)}. "
            "Proceeding anyway."
        )

    if not anchors or not others:
        print(
            "[WARNING] Insufficient anchors or non-anchor sequences. "
            "Skipping MAFFT plot."
        )
        return

    aln_len = alignment.get_alignment_length()
    identity_threshold = identity_percent * window / 100.0 / 2.0

    for other in others:
        # Map alignment columns to non-anchor coordinates
        col_to_coord = {}
        coord = 0
        for i in range(aln_len):
            if other.seq[i] != "-":
                coord += 1
                col_to_coord[i] = coord

        max_coord = coord
        if max_coord == 0:
            continue

        plt.figure(figsize=(12, 3))

        # Y positions
        y_other = 0
        y_anchors = list(range(1, 1 + len(anchors)))

        # Draw non-anchor baseline
        plt.hlines(
            y=y_other,
            xmin=1,
            xmax=max_coord,
            linewidth=6,
            color="black",
            label=other.id
        )

        for anchor, y in zip(anchors, y_anchors):
            identity = []

            for i in range(aln_len):
                if i not in col_to_coord:
                    identity.append(None)
                    continue

                a = anchor.seq[i]
                b = other.seq[i]

                if a == "-" or b == "-":
                    identity.append(0)
                else:
                    identity.append(100 if a == b else 0)

            # Sliding window averaging
            smoothed = []
            for i in range(len(identity)):
                window_vals = [
                    v for v in identity[max(0, i - window // 2): i + window // 2]
                    if v is not None
                ]
                smoothed.append(
                    sum(window_vals) / len(window_vals)
                    if window_vals else 0
                )

            # Extract high-identity regions
            start = None
            for i, val in enumerate(smoothed):
                if val >= identity_threshold and start is None:
                    start = i
                elif val < identity_threshold and start is not None:
                    x1 = col_to_coord.get(start)
                    x2 = col_to_coord.get(i - 1)
                    if x1 and x2:
                        plt.hlines(y, x1, x2, linewidth=6)
                    start = None

            if start is not None:
                x1 = col_to_coord.get(start)
                x2 = col_to_coord.get(len(smoothed) - 1)
                if x1 and x2:
                    plt.hlines(y, x1, x2, linewidth=6)

        plt.yticks(
            [y_other] + y_anchors,
            [other.id] + [a.id for a in anchors]
        )

        plt.xlabel("Non-anchor sequence coordinate")
        plt.ylabel("Sequences")
        plt.xlim(-max_coord * 0.02, max_coord * 1.02)
        plt.ylim(-1, max(y_anchors) + 1)
        plt.title(f"Anchor alignment regions (≥{identity_percent}% " +\
                  "identity, 50 bp window)")
        plt.tight_layout()

        out_file = f"{out_png}.{other.id}.png"
        plt.savefig(out_file)
        plt.close()

        print(f"[INFO] Saved MAFFT plot: {out_file}")


def main():
    args = parse_args()

    print("[INFO] Starting alignment plotting workflow.")
    print(f"[INFO] Output prefix: {args.out_prefix}")

    query_lengths_dict = get_query_lengths_dict(
        mafft_fasta=args.mafft_fasta if args.mafft_fasta and os.path.exists(args.mafft_fasta) else None,
        fasta_contigs=args.fasta_contigs if args.fasta_contigs and os.path.exists(args.fasta_contigs) else None
    )

    if args.blast_tsv and os.path.exists(args.blast_tsv):
        plot_blast(
            args.blast_tsv,
            f"{args.out_prefix}.blast",
            query_lengths_dict
        )
    else:
        print(
            "[INFO] BLAST input not provided or not found.",
            "Skipping BLAST plot."
        )

    if args.mafft_fasta and os.path.exists(args.mafft_fasta):
        plot_mafft_identity(
            args.mafft_fasta,
            f"{args.out_prefix}.mafft_identity",
            args.anchor_id_substring
        )
    else:
        print(
            "[INFO] MAFFT input not provided or not found.",
            "Skipping MAFFT plot."
        )

    print("[INFO] Done.")
    return None


if __name__ == "__main__":
    main()
