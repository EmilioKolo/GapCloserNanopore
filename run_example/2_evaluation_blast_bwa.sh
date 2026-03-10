#!/bin/bash

# To be run with the BWA alignments instead of Minimap2
# Align the gap region to the reference assembly and visualize the alignments
# (Blast version)

# Make the script stop if any command fails
set -e

# Create fasta files with selected regions for gap and telomeres
#python extract_fasta_regions.py \
#    --input ./input/C_glabrata_CBS138_current_chromosomes.fasta \
#    --region ChrC_C_glabrata_CBS138:95600-100700 \
#    --region ChrC_C_glabrata_CBS138:101600-106700 \
#    --output-folder results/ \
#    --output-name gap_5k.fasta

# Align the gap region to the reference assembly
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3511_bwa_gap_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3511_bwa_gap_unmap_gfa \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3511_bwa_gap_no_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3511_bwa_gap_no_unmap_gfa \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3512_bwa_gap_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3512_bwa_gap_unmap_gfa \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3512_bwa_gap_no_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3512_bwa_gap_no_unmap_gfa \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3567_bwa_gap_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3567_bwa_gap_unmap_gfa \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3567_bwa_gap_no_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3567_bwa_gap_no_unmap_gfa \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3569_bwa_gap_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3569_bwa_gap_unmap_gfa \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/runs_3569_bwa_gap_no_unmap_gfa.fasta \
    --out-prefix ./results/alignments/runs_3569_bwa_gap_no_unmap_gfa \
    --mafft-threads 4

# Visualize the alignments

python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3511_bwa_gap_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3511_bwa_gap_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3511_bwa_gap_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3511_bwa_gap_no_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3511_bwa_gap_no_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3511_bwa_gap_no_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3512_bwa_gap_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3512_bwa_gap_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3512_bwa_gap_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3512_bwa_gap_no_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3512_bwa_gap_no_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3512_bwa_gap_no_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3567_bwa_gap_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3567_bwa_gap_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3567_bwa_gap_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3567_bwa_gap_no_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3567_bwa_gap_no_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3567_bwa_gap_no_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3569_bwa_gap_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3569_bwa_gap_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3569_bwa_gap_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/alignments/runs_3569_bwa_gap_no_unmap_gfa.blast.tsv \
    --mafft-fasta ./results/alignments/runs_3569_bwa_gap_no_unmap_gfa.mafft.fasta \
    --out-prefix ./results/alignments/runs_3569_bwa_gap_no_unmap_gfa_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
