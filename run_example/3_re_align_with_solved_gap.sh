#!/bin/bash

# Re-align to the reference assembly after solving the gap and visualize the alignments

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

python extract_fasta_regions.py \
    --input ./results/3512_gap_unmap.fasta \
    --region Ctg18:49492-51594 \
    --output-folder results/ \
    --output-name 3512_gap_unmap_solved1.fasta

python extract_fasta_regions.py \
    --input ./results/3567_gap_unmap.fasta \
    --region Utg28:78514-80624 \
    --output-folder results/ \
    --output-name 3567_gap_unmap_solved1.fasta
python extract_fasta_regions.py \
    --input ./results/3567_gap_unmap.fasta \
    --region Utg28:99903-102012 \
    --output-folder results/ \
    --output-name 3567_gap_unmap_solved2.fasta

python extract_fasta_regions.py \
    --input ./results/3569_gap_unmap.fasta \
    --region Utg44:51614-53724 \
    --output-folder results/ \
    --output-name 3569_gap_unmap_solved1.fasta
python extract_fasta_regions.py \
    --input ./results/3569_gap_unmap.fasta \
    --region Utg44:73003-75118 \
    --output-folder results/ \
    --output-name 3569_gap_unmap_solved2.fasta

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

python replace_region_in_fasta.py \
    --reference ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --insert ./results/3512_gap_unmap_solved1.fasta \
    --chrom ChrC_C_glabrata_CBS138 \
    --start 100663 --end 101662 \
    --output ./results/reference_3512_gap_unmap_solved1.fasta

python replace_region_in_fasta.py \
    --reference ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --insert ./results/3567_gap_unmap_solved1.fasta \
    --chrom ChrC_C_glabrata_CBS138 \
    --start 100663 --end 101662 \
    --output ./results/reference_3567_gap_unmap_solved1.fasta
python replace_region_in_fasta.py \
    --reference ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --insert ./results/3567_gap_unmap_solved2.fasta \
    --chrom ChrC_C_glabrata_CBS138 \
    --start 100663 --end 101662 \
    --output ./results/reference_3567_gap_unmap_solved2.fasta

python replace_region_in_fasta.py \
    --reference ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --insert ./results/3569_gap_unmap_solved1.fasta \
    --chrom ChrC_C_glabrata_CBS138 \
    --start 100663 --end 101662 \
    --output ./results/reference_3569_gap_unmap_solved1.fasta
python replace_region_in_fasta.py \
    --reference ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --insert ./results/3569_gap_unmap_solved2.fasta \
    --chrom ChrC_C_glabrata_CBS138 \
    --start 100663 --end 101662 \
    --output ./results/reference_3569_gap_unmap_solved2.fasta

# Re-align to create bam files
minimap2 -ax map-ont \
    ./results/reference_3511_gap_unmap_solved1.fasta \
    ./input/runs_3511.fastq.gz \
    | samtools sort -o ./results/remap_3511_gap_unmap_solved1.bam
minimap2 -ax map-ont \
    ./results/reference_3511_gap_unmap_solved2.fasta \
    ./input/runs_3511.fastq.gz \
    | samtools sort -o ./results/remap_3511_gap_unmap_solved2.bam

minimap2 -ax map-ont \
    ./results/reference_3512_gap_unmap_solved1.fasta \
    ./input/runs_3512.fastq.gz \
    | samtools sort -o ./results/remap_3512_gap_unmap_solved1.bam

minimap2 -ax map-ont \
    ./results/reference_3567_gap_unmap_solved1.fasta \
    ./input/runs_3567.fastq.gz \
    | samtools sort -o ./results/remap_3567_gap_unmap_solved1.bam
minimap2 -ax map-ont \
    ./results/reference_3567_gap_unmap_solved2.fasta \
    ./input/runs_3567.fastq.gz \
    | samtools sort -o ./results/remap_3567_gap_unmap_solved2.bam

minimap2 -ax map-ont \
    ./results/reference_3569_gap_unmap_solved1.fasta \
    ./input/runs_3569.fastq.gz \
    | samtools sort -o ./results/remap_3569_gap_unmap_solved1.bam
minimap2 -ax map-ont \
    ./results/reference_3569_gap_unmap_solved2.fasta \
    ./input/runs_3569.fastq.gz \
    | samtools sort -o ./results/remap_3569_gap_unmap_solved2.bam

# Index bam files
samtools index ./results/remap_3511_gap_unmap_solved1.bam
samtools index ./results/remap_3511_gap_unmap_solved2.bam

samtools index ./results/remap_3512_gap_unmap_solved1.bam

samtools index ./results/remap_3567_gap_unmap_solved1.bam
samtools index ./results/remap_3567_gap_unmap_solved2.bam

samtools index ./results/remap_3569_gap_unmap_solved1.bam
samtools index ./results/remap_3569_gap_unmap_solved2.bam

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

python $PATH_TO_MAIN \
    --bam ./results/remap_3512_gap_unmap_solved1.bam \
    --fasta ./results/reference_3512_gap_unmap_solved1.fasta \
    --out ./results/remap_3512_gap_unmap_solved1_coverage.png \
    --binsize 10

python $PATH_TO_MAIN \
    --bam ./results/remap_3567_gap_unmap_solved1.bam \
    --fasta ./results/reference_3567_gap_unmap_solved1.fasta \
    --out ./results/remap_3567_gap_unmap_solved1_coverage.png \
    --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/remap_3567_gap_unmap_solved2.bam \
    --fasta ./results/reference_3567_gap_unmap_solved2.fasta \
    --out ./results/remap_3567_gap_unmap_solved2_coverage.png \
    --binsize 10

python $PATH_TO_MAIN \
    --bam ./results/remap_3569_gap_unmap_solved1.bam \
    --fasta ./results/reference_3569_gap_unmap_solved1.fasta \
    --out ./results/remap_3569_gap_unmap_solved1_coverage.png \
    --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/remap_3569_gap_unmap_solved2.bam \
    --fasta ./results/reference_3569_gap_unmap_solved2.fasta \
    --out ./results/remap_3569_gap_unmap_solved2_coverage.png \
    --binsize 10
