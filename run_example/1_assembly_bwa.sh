#!/bin/bash

# To be run with the BWA alignments instead of Minimap2
# Cleanup for fastq files with empty reads and initial assembly with Raven

# Make the script stop if any command fails
set -e

# Activate Raven in conda (better to do it before running this script)
#conda init
#conda activate raven_assembler

# Cleanup fastq files

python remove_empty_reads.py \
    --input ./results/runs_3511_bwa_gap_no_unmap.fastq.gz \
    --output ./results/runs_3511_bwa_gap_no_unmap_clean.fastq
python remove_empty_reads.py \
    --input ./results/runs_3511_bwa_gap_unmap.fastq.gz \
    --output ./results/runs_3511_bwa_gap_unmap_clean.fastq

python remove_empty_reads.py \
    --input ./results/runs_3512_bwa_gap_no_unmap.fastq.gz \
    --output ./results/runs_3512_bwa_gap_no_unmap_clean.fastq
python remove_empty_reads.py \
    --input ./results/runs_3512_bwa_gap_unmap.fastq.gz \
    --output ./results/runs_3512_bwa_gap_unmap_clean.fastq

python remove_empty_reads.py \
    --input ./results/runs_3567_bwa_gap_no_unmap.fastq.gz \
    --output ./results/runs_3567_bwa_gap_no_unmap_clean.fastq
python remove_empty_reads.py \
    --input ./results/runs_3567_bwa_gap_unmap.fastq.gz \
    --output ./results/runs_3567_bwa_gap_unmap_clean.fastq

python remove_empty_reads.py \
    --input ./results/runs_3569_bwa_gap_no_unmap.fastq.gz \
    --output ./results/runs_3569_bwa_gap_no_unmap_clean.fastq
python remove_empty_reads.py \
    --input ./results/runs_3569_bwa_gap_unmap.fastq.gz \
    --output ./results/runs_3569_bwa_gap_unmap_clean.fastq

# Run raven

raven -t 6 -F ./results/runs_3511_bwa_gap_no_unmap.gfa \
    ./results/runs_3511_bwa_gap_no_unmap_clean.fastq \
    > ./results/runs_3511_bwa_gap_no_unmap.fasta
raven -t 6 -F ./results/runs_3511_bwa_gap_unmap.gfa \
    ./results/runs_3511_bwa_gap_unmap_clean.fastq \
    > ./results/runs_3511_bwa_gap_unmap.fasta

raven -t 6 -F ./results/runs_3512_bwa_gap_no_unmap.gfa \
    ./results/runs_3512_bwa_gap_no_unmap_clean.fastq \
    > ./results/runs_3512_bwa_gap_no_unmap.fasta
raven -t 6 -F ./results/runs_3512_bwa_gap_unmap.gfa \
    ./results/runs_3512_bwa_gap_unmap_clean.fastq \
    > ./results/runs_3512_bwa_gap_unmap.fasta

raven -t 6 -F ./results/runs_3567_bwa_gap_no_unmap.gfa \
    ./results/runs_3567_bwa_gap_no_unmap_clean.fastq \
    > ./results/runs_3567_bwa_gap_no_unmap.fasta
raven -t 6 -F ./results/runs_3567_bwa_gap_unmap.gfa \
    ./results/runs_3567_bwa_gap_unmap_clean.fastq \
    > ./results/runs_3567_bwa_gap_unmap.fasta

raven -t 6 -F ./results/runs_3569_bwa_gap_no_unmap.gfa \
    ./results/runs_3569_bwa_gap_no_unmap_clean.fastq \
    > ./results/runs_3569_bwa_gap_no_unmap.fasta
raven -t 6 -F ./results/runs_3569_bwa_gap_unmap.gfa \
    ./results/runs_3569_bwa_gap_unmap_clean.fastq \
    > ./results/runs_3569_bwa_gap_unmap.fasta
