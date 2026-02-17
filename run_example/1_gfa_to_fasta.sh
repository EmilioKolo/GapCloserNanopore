#!/bin/bash

# Make the script stop if any command fails
set -e

## OPTIONAL: Turn gfa files into fasta files
python gfa_to_fasta.py \
    --input ./results/runs_3511_minimap_gap_no_unmap.gfa \
    --output-name 3511_gap_no_unmap.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3511_minimap_gap_unmap.gfa \
    --output-name 3511_gap_unmap.fasta \
    --output-folder results/

python gfa_to_fasta.py \
    --input ./results/runs_3512_minimap_gap_no_unmap.gfa \
    --output-name 3512_gap_no_unmap.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3512_minimap_gap_unmap.gfa \
    --output-name 3512_gap_unmap.fasta \
    --output-folder results/

python gfa_to_fasta.py \
    --input ./results/runs_3567_minimap_gap_no_unmap.gfa \
    --output-name 3567_gap_no_unmap.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3567_minimap_gap_unmap.gfa \
    --output-name 3567_gap_unmap.fasta \
    --output-folder results/

python gfa_to_fasta.py \
    --input ./results/runs_3569_minimap_gap_no_unmap.gfa \
    --output-name 3569_gap_no_unmap.fasta \
    --output-folder results/
python gfa_to_fasta.py \
    --input ./results/runs_3569_minimap_gap_unmap.gfa \
    --output-name 3569_gap_unmap.fasta \
    --output-folder results/
