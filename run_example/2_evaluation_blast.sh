#!/bin/bash

# Align the gap region to the reference assembly and visualize the alignments
# (Blast version)

# Make the script stop if any command fails
set -e

# Create fasta files with selected regions for gap and telomeres
python extract_fasta_regions.py \
    --input ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name gap_5k.fasta

# Align the gap region to the reference assembly
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3511_gap_unmap.fasta \
    --out-prefix ./results/alignments/3511_gap_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3511_gap_no_unmap.fasta \
    --out-prefix ./results/alignments/3511_gap_no_unmap \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3512_gap_unmap.fasta \
    --out-prefix ./results/alignments/3512_gap_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3512_gap_no_unmap.fasta \
    --out-prefix ./results/alignments/3512_gap_no_unmap \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3567_gap_unmap.fasta \
    --out-prefix ./results/alignments/3567_gap_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3567_gap_no_unmap.fasta \
    --out-prefix ./results/alignments/3567_gap_no_unmap \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3569_gap_unmap.fasta \
    --out-prefix ./results/alignments/3569_gap_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/3569_gap_no_unmap.fasta \
    --out-prefix ./results/alignments/3569_gap_no_unmap \
    --mafft-threads 4

# Visualize the alignments

python plot_alignments.py \
    --blast-tsv ./results/alignments/3511_gap_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3511_gap_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3511_gap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/3511_gap_no_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3511_gap_no_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3511_gap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/alignments/3512_gap_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3512_gap_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3512_gap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/3512_gap_no_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3512_gap_no_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3512_gap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/alignments/3567_gap_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3567_gap_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3567_gap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/3567_gap_no_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3567_gap_no_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3567_gap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/alignments/3569_gap_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3569_gap_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3569_gap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/3569_gap_no_unmap.blast.tsv \
    --mafft-fasta ./results/alignments/3569_gap_no_unmap.mafft.fasta \
    --out-prefix ./results/alignments/3569_gap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
