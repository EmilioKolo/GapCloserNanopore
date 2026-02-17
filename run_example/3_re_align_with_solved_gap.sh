#!/bin/bash

# Make the script stop if any command fails
set -e

# Extract sequence between anchors from contig
python extract_fasta_regions.py \
    --input ./results/3511_gap_unmap.fasta \
    --region Utg34:60744-62854 \
    --output-folder results/ \
    --output-name 3511_gap_unmap_solved1.fasta
python extract_fasta_regions.py \
    --input ./results/3511_gap_unmap.fasta \
    --region Utg34:82156-84272 \
    --output-folder results/ \
    --output-name 3511_gap_unmap_solved2.fasta

# Paste each extracted sequence into the gap region of the original genome
# Gap region: ChrC_C_glabrata_CBS138:100663-101662
python replace_region_in_fasta.py \
    --reference ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --insert ./results/3511_gap_unmap_solved1.fasta \
    --chrom ChrC_C_glabrata_CBS138 \
    --start 100663 --end 101662 \
    --output ./results/reference_3511_gap_unmap_solved1.fasta
python replace_region_in_fasta.py \
    --reference ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --insert ./results/3511_gap_unmap_solved2.fasta \
    --chrom ChrC_C_glabrata_CBS138 \
    --start 100663 --end 101662 \
    --output ./results/reference_3511_gap_unmap_solved2.fasta

# Re-align to create bam files
minimap2 -ax map-ont \
    ./results/reference_3511_gap_unmap_solved1.fasta \
    ./input/runs_3511.fastq.gz \
    | samtools sort -o ./results/remap_3511_gap_unmap_solved1.bam
minimap2 -ax map-ont \
    ./results/reference_3511_gap_unmap_solved2.fasta \
    ./input/runs_3511.fastq.gz \
    | samtools sort -o ./results/remap_3511_gap_unmap_solved2.bam

# Index bam files
samtools index ./results/remap_3511_gap_unmap_solved1.bam
samtools index ./results/remap_3511_gap_unmap_solved2.bam

# Coverage plots
PATH_TO_MAIN=../CoveragePlotter/main.py
python $PATH_TO_MAIN \
    --bam ./results/remap_3511_gap_unmap_solved1.bam \
    --fasta ./results/reference_3511_gap_unmap_solved1.fasta \
    --out ./results/remap_3511_gap_unmap_solved1_coverage.png \
    --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/remap_3511_gap_unmap_solved2.bam \
    --fasta ./results/reference_3511_gap_unmap_solved2.fasta \
    --out ./results/remap_3511_gap_unmap_solved2_coverage.png \
    --binsize 10
