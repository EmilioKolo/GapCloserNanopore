#!/bin/bash

# Extract reads that map to multiple locations in the reference assembly

# Make the script stop if any command fails
set -e

echo "Extracting double-mapping reads from BWA alignments..."

# Extract double-mapping reads with samtools and grep
samtools view -h ./input/runs_3511_bwa_sorted.bam \
    | awk '$2 ~ /256/ || $2 ~ /2048/ || NR<10' \
    | samtools view -bS - > ./results/multi_map/3511_bwa_multi_map.bam

samtools view -h ./input/runs_3512_bwa_sorted.bam \
    | awk '$2 ~ /256/ || $2 ~ /2048/ || NR<10' \
    | samtools view -bS - > ./results/multi_map/3512_bwa_multi_map.bam

samtools view -h ./input/runs_3567_bwa_sorted.bam \
    | awk '$2 ~ /256/ || $2 ~ /2048/ || NR<10' \
    | samtools view -bS - > ./results/multi_map/3567_bwa_multi_map.bam

samtools view -h ./input/runs_3569_bwa_sorted.bam \
    | awk '$2 ~ /256/ || $2 ~ /2048/ || NR<10' \
    | samtools view -bS - > ./results/multi_map/3569_bwa_multi_map.bam

echo "Double-mapping reads extracted and saved as BAM files in ./results/multi_map/"
echo "Proceeding to index BAM files and generate coverage plots..."

# Index the resulting BAM files
samtools index ./results/multi_map/3511_bwa_multi_map.bam
samtools index ./results/multi_map/3512_bwa_multi_map.bam
samtools index ./results/multi_map/3567_bwa_multi_map.bam
samtools index ./results/multi_map/3569_bwa_multi_map.bam

echo "BAM files indexed."
echo "Transforming BAM files to SAM format for easier inspection..."

# Transform BAM to SAM for easier inspection
samtools view -h -o ./results/multi_map/3511_bwa_multi_map.sam \
    ./results/multi_map/3511_bwa_multi_map.bam
samtools view -h -o ./results/multi_map/3512_bwa_multi_map.sam \
    ./results/multi_map/3512_bwa_multi_map.bam
samtools view -h -o ./results/multi_map/3567_bwa_multi_map.sam \
    ./results/multi_map/3567_bwa_multi_map.bam
samtools view -h -o ./results/multi_map/3569_bwa_multi_map.sam \
    ./results/multi_map/3569_bwa_multi_map.bam

echo "BAM files transformed to SAM format for easier inspection."
echo "Proceeding to generate coverage plots for double-mapping reads..."

# Coverage plots (done using CoveragePlotter: https://github.com/EmilioKolo/CoveragePlotter )
PATH_TO_MAIN=../CoveragePlotter/main.py
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3511_bwa_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3511_bwa_multi_map_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3512_bwa_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3512_bwa_multi_map_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3567_bwa_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3567_bwa_multi_map_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3569_bwa_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3569_bwa_multi_map_coverage.png --binsize 10

echo "Coverage plots generated and saved in ./results/multi_map/"
