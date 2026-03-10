#!/bin/bash

# To be run with the BWA alignments instead of Minimap2
# Filtering of reads that map to the gap region with and without unmapped reads

# Make the script stop if any command fails
set -e

# 3511
python extract_reads.py --input input/runs_3511_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3511_bwa_gap_no_unmap.fastq.gz
python extract_reads.py --input input/runs_3511_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3511_bwa_gap_unmap.fastq.gz \
    --include-unmapped

# 3512
python extract_reads.py --input input/runs_3512_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3512_bwa_gap_no_unmap.fastq.gz
python extract_reads.py --input input/runs_3512_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3512_bwa_gap_unmap.fastq.gz \
    --include-unmapped

# 3567
python extract_reads.py --input input/runs_3567_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3567_bwa_gap_no_unmap.fastq.gz
python extract_reads.py --input input/runs_3567_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3567_bwa_gap_unmap.fastq.gz \
    --include-unmapped

# 3569
python extract_reads.py --input input/runs_3569_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3569_bwa_gap_no_unmap.fastq.gz
python extract_reads.py --input input/runs_3569_bwa_sorted.bam \
    --region ChrC_C_glabrata_CBS138:95600-100700 \
    --region ChrC_C_glabrata_CBS138:101600-106700 \
    --output-folder results/ \
    --output-name runs_3569_bwa_gap_unmap.fastq.gz \
    --include-unmapped