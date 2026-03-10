#!/bin/bash

# Make output directories
mkdir -p ./results/unmapped_reads
mkdir -p ./results/unmapped_reads/fastqc

# Extract unmapped reads from the BAM file
samtools view -b -f 4 ./input/runs_3511_bwa_sorted.bam \
    > ./results/unmapped_reads/3511_bwa_sorted_unmapped.bam
samtools view -b -f 4 ./input/runs_3511_minimap_sorted.bam \
    > ./results/unmapped_reads/3511_minimap_sorted_unmapped.bam

bamtools filter -in ./results/unmapped_reads/3511_bwa_sorted_unmapped.bam \
    -out ./results/unmapped_reads/3511_bwa_sorted_unmapped_filtered.bam \
    -script <(echo '{"filters": [{"property": "length", "max": 100000}]}')
bamtools filter -in ./results/unmapped_reads/3511_minimap_sorted_unmapped.bam \
    -out ./results/unmapped_reads/3511_minimap_sorted_unmapped_filtered.bam \
    -script <(echo '{"filters": [{"property": "length", "max": 100000}]}')


samtools view -b -f 4 ./input/runs_3512_bwa_sorted.bam \
    > ./results/unmapped_reads/3512_bwa_sorted_unmapped.bam
samtools view -b -f 4 ./input/runs_3512_minimap_sorted.bam \
    > ./results/unmapped_reads/3512_minimap_sorted_unmapped.bam

samtools view -b -f 4 ./input/runs_3567_bwa_sorted.bam \
    > ./results/unmapped_reads/3567_bwa_sorted_unmapped.bam
samtools view -b -f 4 ./input/runs_3567_minimap_sorted.bam \
    > ./results/unmapped_reads/3567_minimap_sorted_unmapped.bam

samtools view -b -f 4 ./input/runs_3569_bwa_sorted.bam \
    > ./results/unmapped_reads/3569_bwa_sorted_unmapped.bam
samtools view -b -f 4 ./input/runs_3569_minimap_sorted.bam \
    > ./results/unmapped_reads/3569_minimap_sorted_unmapped.bam

# Run through Fastqc
fastqc -o ./results/unmapped_reads/fastqc \
    ./results/unmapped_reads/*.bam
