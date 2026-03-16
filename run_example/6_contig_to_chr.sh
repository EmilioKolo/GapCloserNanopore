#!/bin/bash

# Make output directory
mkdir -p ./results/low_qual_filt_reads/realignment

# Align de novo contigs to Cglab reference genome
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3511_bwa_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3511_bwa_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3511_bwa_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3511_bwa_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3511_bwa_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3511_bwa_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3511_bwa_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3511_bwa_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3511_bwa_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3511_bwa_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3511_bwa_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3511_bwa_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3511_bwa_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3511_bwa_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3511_minimap_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3511_minimap_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3511_minimap_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3511_minimap_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3511_minimap_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3511_minimap_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3511_minimap_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3511_minimap_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3511_minimap_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3511_minimap_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3511_minimap_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3511_minimap_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3511_minimap_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3511_minimap_unmap_gfa.sorted.bam

minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3512_bwa_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3512_bwa_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3512_bwa_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3512_bwa_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3512_bwa_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3512_bwa_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3512_bwa_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3512_bwa_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3512_bwa_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3512_bwa_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3512_bwa_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3512_bwa_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3512_bwa_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3512_bwa_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3512_minimap_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3512_minimap_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3512_minimap_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3512_minimap_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3512_minimap_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3512_minimap_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3512_minimap_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3512_minimap_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3512_minimap_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3512_minimap_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3512_minimap_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3512_minimap_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3512_minimap_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3512_minimap_unmap_gfa.sorted.bam

minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3567_bwa_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3567_bwa_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3567_bwa_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3567_bwa_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3567_bwa_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3567_bwa_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3567_bwa_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3567_bwa_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3567_bwa_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3567_bwa_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3567_bwa_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3567_bwa_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3567_bwa_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3567_bwa_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3567_minimap_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3567_minimap_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3567_minimap_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3567_minimap_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3567_minimap_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3567_minimap_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3567_minimap_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3567_minimap_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3567_minimap_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3567_minimap_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3567_minimap_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3567_minimap_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3567_minimap_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3567_minimap_unmap_gfa.sorted.bam

minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3569_bwa_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3569_bwa_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3569_bwa_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3569_bwa_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3569_bwa_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3569_bwa_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3569_bwa_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3569_bwa_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3569_bwa_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3569_bwa_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3569_bwa_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3569_bwa_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3569_bwa_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3569_bwa_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3569_minimap_no_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3569_minimap_no_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3569_minimap_no_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3569_minimap_no_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3569_minimap_no_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3569_minimap_no_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3569_minimap_no_unmap_gfa.sorted.bam
minimap2 -ax asm5 \
    ./input/C_glabrata_CBS138_ChrC.fasta \
    ./results/low_qual_filt_reads/runs_3569_minimap_unmap_gfa.fasta > \
    ./results/low_qual_filt_reads/realignment/runs_3569_minimap_unmap_gfa.sam
samtools view -S \
    -b ./results/low_qual_filt_reads/realignment/runs_3569_minimap_unmap_gfa.sam \
    > ./results/low_qual_filt_reads/realignment/runs_3569_minimap_unmap_gfa.bam
samtools sort \
    ./results/low_qual_filt_reads/realignment/runs_3569_minimap_unmap_gfa.bam \
    -o ./results/low_qual_filt_reads/realignment/runs_3569_minimap_unmap_gfa.sorted.bam
samtools index ./results/low_qual_filt_reads/realignment/runs_3569_minimap_unmap_gfa.sorted.bam
