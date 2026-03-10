#!/bin/bash

# To be run with the BWA alignments instead of Minimap2
# Transform gfa files into fasta files

# Make the script stop if any command fails
set -e

## OPTIONAL: Turn gfa files into fasta files
python gfa_to_fasta.py \
    --input ./results/runs_3511_bwa_gap_no_unmap.gfa \
    --output-name runs_3511_bwa_gap_no_unmap_gfa.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3511_bwa_gap_unmap.gfa \
    --output-name runs_3511_bwa_gap_unmap_gfa.fasta \
    --output-folder results/

python gfa_to_fasta.py \
    --input ./results/runs_3512_bwa_gap_no_unmap.gfa \
    --output-name runs_3512_bwa_gap_no_unmap_gfa.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3512_bwa_gap_unmap.gfa \
    --output-name runs_3512_bwa_gap_unmap_gfa.fasta \
    --output-folder results/

python gfa_to_fasta.py \
    --input ./results/runs_3567_bwa_gap_no_unmap.gfa \
    --output-name runs_3567_bwa_gap_no_unmap_gfa.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3567_bwa_gap_unmap.gfa \
    --output-name runs_3567_bwa_gap_unmap_gfa.fasta \
    --output-folder results/

python gfa_to_fasta.py \
    --input ./results/runs_3569_bwa_gap_no_unmap.gfa \
    --output-name runs_3569_bwa_gap_no_unmap_gfa.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3569_bwa_gap_unmap.gfa \
    --output-name runs_3569_bwa_gap_unmap_gfa.fasta \
    --output-folder results/
