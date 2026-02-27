#!/bin/bash

# Extract reads that map to multiple locations in the reference assembly

# Make the script stop if any command fails
set -e

# Extract double-mapping reads with samtools and grep
samtools view -h ./input/runs_3511_minimap_sorted.bam \
    | grep -E "^@|NH:i:[2-9]" \
    | samtools view -bS - > ./results/multi_map/3511_minimap_multi_map.bam

samtools view -h ./input/runs_3512_minimap_sorted.bam \
    | grep -E "^@|NH:i:[2-9]" \
    | samtools view -bS - > ./results/multi_map/3512_minimap_multi_map.bam

samtools view -h ./input/runs_3567_minimap_sorted.bam \
    | grep -E "^@|NH:i:[2-9]" \
    | samtools view -bS - > ./results/multi_map/3567_minimap_multi_map.bam

samtools view -h ./input/runs_3569_minimap_sorted.bam \
    | grep -E "^@|NH:i:[2-9]" \
    | samtools view -bS - > ./results/multi_map/3569_minimap_multi_map.bam

# Index the resulting BAM files
samtools index ./results/multi_map/3511_minimap_multi_map.bam
samtools index ./results/multi_map/3512_minimap_multi_map.bam
samtools index ./results/multi_map/3567_minimap_multi_map.bam
samtools index ./results/multi_map/3569_minimap_multi_map.bam

# Transform BAM to SAM for easier inspection
samtools view -h -o ./results/multi_map/3511_minimap_multi_map.sam \
    ./results/multi_map/3511_minimap_multi_map.bam
samtools view -h -o ./results/multi_map/3512_minimap_multi_map.sam \
    ./results/multi_map/3512_minimap_multi_map.bam
samtools view -h -o ./results/multi_map/3567_minimap_multi_map.sam \
    ./results/multi_map/3567_minimap_multi_map.bam
samtools view -h -o ./results/multi_map/3569_minimap_multi_map.sam \
    ./results/multi_map/3569_minimap_multi_map.bam

# Coverage plots (done using CoveragePlotter: https://github.com/EmilioKolo/CoveragePlotter )
PATH_TO_MAIN=../CoveragePlotter/main.py
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3511_minimap_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3511_minimap_multi_map_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3512_minimap_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3512_minimap_multi_map_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3567_minimap_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3567_minimap_multi_map_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ./results/multi_map/3569_minimap_multi_map.bam \
    --fasta ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --out ./results/multi_map/3569_minimap_multi_map_coverage.png --binsize 10
