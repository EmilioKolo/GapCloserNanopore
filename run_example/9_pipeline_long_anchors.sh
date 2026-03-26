#!/bin/bash

### IMPORTANT: Activate venv and conda flye environment before running this script

# Make the script stop if any command fails
set -e

RUN_NUM=3511
RUN_NUM=3512
RUN_NUM=3567
RUN_NUM=3569

OUTPUT_DIR=results/full_pipeline_mar_26/$RUN_NUM

mkdir -p ${OUTPUT_DIR}

# Extract the reads that map to the defined gap region
python extract_reads.py --input input/runs_${RUN_NUM}_minimap_sorted.bam \
    --region ChrC_C_glabrata_CBS138:90000-150000 \
    --output-folder ${OUTPUT_DIR}/ \
    --output-name ${RUN_NUM}_minimap_unmap.fastq.gz \
    --include-unmapped

# Create fasta file with selected region for the gap
python extract_fasta_regions.py \
    --input ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --region ChrC_C_glabrata_CBS138:80000-100662 \
    --region ChrC_C_glabrata_CBS138:101663-160000 \
    --output-folder input/ \
    --output-name gap_80k_to_160k.fasta

# Extract fastq.gz file
gunzip -c ${OUTPUT_DIR}/${RUN_NUM}_minimap_unmap.fastq.gz \
    > ${OUTPUT_DIR}/${RUN_NUM}_minimap_unmap.fastq

# Use seqkit sana to properly format file
seqkit sana -j 6 \
    ${OUTPUT_DIR}/${RUN_NUM}_minimap_unmap.fastq \
    -o ${OUTPUT_DIR}/${RUN_NUM}_minimap_unmap_sana.fastq
# Remove empty and low-quality reads
seqkit seq -Q 20 -M 40000 -m 10 -g \
    ${OUTPUT_DIR}/${RUN_NUM}_minimap_unmap_sana.fastq \
    | gzip > ${OUTPUT_DIR}/runs_${RUN_NUM}_minimap_unmap_clean.fastq.gz

mkdir -p ${OUTPUT_DIR}/flye_assembly
mkdir -p ${OUTPUT_DIR}/raven_assembly

# Assemble with Flye and Raven
raven -t 6 -F \
    ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap.gfa \
    ${OUTPUT_DIR}/runs_${RUN_NUM}_minimap_unmap_clean.fastq.gz \
    > ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap.fasta
flye \
    --read-error 0.03 \
    --nano-hq ${OUTPUT_DIR}/runs_${RUN_NUM}_minimap_unmap_clean.fastq.gz \
    --out-dir ${OUTPUT_DIR}/flye_assembly \
    --threads 6

# Turn Raven gfa file into fasta file
python gfa_to_fasta.py \
    --input ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap.gfa \
    --output-name runs_${RUN_NUM}_minimap_unmap_gfa.fasta \
    --output-folder ${OUTPUT_DIR}/raven_assembly/

# Map reads to contigs and see depth
minimap2 -t 6 -ax map-ont \
    ${OUTPUT_DIR}/flye_assembly/assembly.fasta \
    ${OUTPUT_DIR}/runs_${RUN_NUM}_minimap_unmap_clean.fastq.gz \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/flye_assembly/reads_to_contigs.bam
samtools index ${OUTPUT_DIR}/flye_assembly/reads_to_contigs.bam
samtools depth ${OUTPUT_DIR}/flye_assembly/reads_to_contigs.bam \
    > ${OUTPUT_DIR}/flye_assembly/reads_to_contigs.depth
minimap2 -t 6 -ax map-ont \
    ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap_gfa.fasta \
    ${OUTPUT_DIR}/runs_${RUN_NUM}_minimap_unmap_clean.fastq.gz \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/raven_assembly/reads_to_contigs.bam
samtools index ${OUTPUT_DIR}/raven_assembly/reads_to_contigs.bam
samtools depth ${OUTPUT_DIR}/raven_assembly/reads_to_contigs.bam \
    > ${OUTPUT_DIR}/raven_assembly/reads_to_contigs.depth

# Coverage plots (done using CoveragePlotter: https://github.com/EmilioKolo/CoveragePlotter )
PATH_TO_MAIN=../CoveragePlotter/main.py
python $PATH_TO_MAIN \
    --bam ${OUTPUT_DIR}/flye_assembly/reads_to_contigs.bam \
    --fasta ${OUTPUT_DIR}/flye_assembly/assembly.fasta \
    --out ${OUTPUT_DIR}/flye_assembly/reads_to_contigs_coverage.png --binsize 10
python $PATH_TO_MAIN \
    --bam ${OUTPUT_DIR}/raven_assembly/reads_to_contigs.bam \
    --fasta ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap_gfa.fasta \
    --out ${OUTPUT_DIR}/raven_assembly/reads_to_contigs_coverage.png --binsize 10

# Map assembled contigs to reference genome
minimap2 -ax asm5 -t 6 \
    input/C_glabrata_CBS138_current_chromosomes.fasta \
    ${OUTPUT_DIR}/flye_assembly/assembly.fasta \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/flye_assembly/contigs_to_reference.bam
samtools index ${OUTPUT_DIR}/flye_assembly/contigs_to_reference.bam
minimap2 -ax asm5 -t 6 \
    input/C_glabrata_CBS138_current_chromosomes.fasta \
    ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap_gfa.fasta \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/raven_assembly/contigs_to_reference.bam
samtools index ${OUTPUT_DIR}/raven_assembly/contigs_to_reference.bam

# Map anchor reference sequence to assembled contigs
minimap2 -ax asm5 -t 6 \
    ${OUTPUT_DIR}/flye_assembly/assembly.fasta \
    input/reference_anchors.fasta \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/flye_assembly/anchors_to_contigs.bam
samtools index ${OUTPUT_DIR}/flye_assembly/anchors_to_contigs.bam
minimap2 -ax asm5 -t 6 \
    ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap_gfa.fasta \
    input/reference_anchors.fasta \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/raven_assembly/anchors_to_contigs.bam
samtools index ${OUTPUT_DIR}/raven_assembly/anchors_to_contigs.bam

# Map ChrC reference sequence to assembled contigs
minimap2 -ax asm5 -t 6 \
    ${OUTPUT_DIR}/flye_assembly/assembly.fasta \
    input/C_glabrata_CBS138_ChrC.fasta \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/flye_assembly/chrC_to_contigs.bam
samtools index ${OUTPUT_DIR}/flye_assembly/chrC_to_contigs.bam
minimap2 -ax asm5 -t 6 \
    ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap_gfa.fasta \
    input/C_glabrata_CBS138_ChrC.fasta \
    | samtools view -bS - \
    | samtools sort -o ${OUTPUT_DIR}/raven_assembly/chrC_to_contigs.bam
samtools index ${OUTPUT_DIR}/raven_assembly/chrC_to_contigs.bam

# Align the gap region to the reference assembly
python align_fastas.py \
    --fasta1 ./input/gap_80k_to_160k.fasta \
    --fasta2 ${OUTPUT_DIR}/flye_assembly/assembly.fasta \
    --out-prefix ${OUTPUT_DIR}/flye_assembly/runs_${RUN_NUM}_minimap_unmap_gfa_align \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./input/gap_80k_to_160k.fasta \
    --fasta2 ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap_gfa.fasta \
    --out-prefix ${OUTPUT_DIR}/raven_assembly/runs_${RUN_NUM}_minimap_unmap_gfa_align \
    --mafft-threads 4
