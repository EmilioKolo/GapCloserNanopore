#!/bin/bash

# Align the gap region to the reference assembly and visualize the alignments 
# (minimap2 version)

# Make the script stop if any command fails
set -e

# Create bam files
minimap2 -ax map-ont \
    ./results/3511_gap_no_unmap.fasta \
    ./results/runs_3511_minimap_gap_no_unmap_clean.fastq \
    | samtools sort -o ./results/3511_gap_no_unmap_minimap.bam
minimap2 -ax map-ont \
    ./results/3511_gap_unmap.fasta \
    ./results/runs_3511_minimap_gap_unmap_clean.fastq \
    | samtools sort -o ./results/3511_gap_unmap_minimap.bam

minimap2 -ax map-ont \
    ./results/3512_gap_no_unmap.fasta \
    ./results/runs_3512_minimap_gap_no_unmap_clean.fastq \
    | samtools sort -o ./results/3512_gap_no_unmap_minimap.bam
minimap2 -ax map-ont \
    ./results/3512_gap_unmap.fasta \
    ./results/runs_3512_minimap_gap_unmap_clean.fastq \
    | samtools sort -o ./results/3512_gap_unmap_minimap.bam

minimap2 -ax map-ont \
    ./results/3567_gap_no_unmap.fasta \
    ./results/runs_3567_minimap_gap_no_unmap_clean.fastq \
    | samtools sort -o ./results/3567_gap_no_unmap_minimap.bam
minimap2 -ax map-ont \
    ./results/3567_gap_unmap.fasta \
    ./results/runs_3567_minimap_gap_unmap_clean.fastq \
    | samtools sort -o ./results/3567_gap_unmap_minimap.bam

minimap2 -ax map-ont \
    ./results/3569_gap_no_unmap.fasta \
    ./results/runs_3569_minimap_gap_no_unmap_clean.fastq \
    | samtools sort -o ./results/3569_gap_no_unmap_minimap.bam
minimap2 -ax map-ont \
    ./results/3569_gap_unmap.fasta \
    ./results/runs_3569_minimap_gap_unmap_clean.fastq \
    | samtools sort -o ./results/3569_gap_unmap_minimap.bam

# Index bam files
samtools index ./results/3511_gap_no_unmap_minimap.bam
samtools index ./results/3511_gap_unmap_minimap.bam
samtools index ./results/3512_gap_no_unmap_minimap.bam
samtools index ./results/3512_gap_unmap_minimap.bam
samtools index ./results/3567_gap_no_unmap_minimap.bam
samtools index ./results/3567_gap_unmap_minimap.bam
samtools index ./results/3569_gap_no_unmap_minimap.bam
samtools index ./results/3569_gap_unmap_minimap.bam

# Coverage plots (done using CoveragePlotter: https://github.com/EmilioKolo/CoveragePlotter )
PATH_TO_MAIN=../CoveragePlotter/main.py
python $PATH_TO_MAIN \
    --bam ./results/3511_gap_no_unmap_minimap.bam \
    --fasta ./results/3511_gap_no_unmap.fasta \
    --out ./results/3511_gap_no_unmap_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/3511_gap_unmap_minimap.bam \
    --fasta ./results/3511_gap_unmap.fasta \
    --out ./results/3511_gap_unmap_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/3512_gap_no_unmap_minimap.bam \
    --fasta ./results/3512_gap_no_unmap.fasta \
    --out ./results/3512_gap_no_unmap_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/3512_gap_unmap_minimap.bam \
    --fasta ./results/3512_gap_unmap.fasta \
    --out ./results/3512_gap_unmap_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/3567_gap_no_unmap_minimap.bam \
    --fasta ./results/3567_gap_no_unmap.fasta \
    --out ./results/3567_gap_no_unmap_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/3567_gap_unmap_minimap.bam \
    --fasta ./results/3567_gap_unmap.fasta \
    --out ./results/3567_gap_unmap_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/3569_gap_no_unmap_minimap.bam \
    --fasta ./results/3569_gap_no_unmap.fasta \
    --out ./results/3569_gap_no_unmap_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/3569_gap_unmap_minimap.bam \
    --fasta ./results/3569_gap_unmap.fasta \
    --out ./results/3569_gap_unmap_coverage.png --binsize 10
