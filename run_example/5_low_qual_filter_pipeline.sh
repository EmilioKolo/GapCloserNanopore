#!/bin/bash

# Make output directories
mkdir -p ./input/low_qual_filt_reads
mkdir -p ./results/low_qual_filt_reads

# Filter low quality reads from inputs
seqkit seq -Q 20 -M 40000 ./input/runs_3511.fastq.gz \
    | gzip > ./input/low_qual_filt_reads/runs_3511.fastq.gz
seqkit seq -Q 20 -M 40000 ./input/runs_3512.fastq.gz \
    | gzip > ./input/low_qual_filt_reads/runs_3512.fastq.gz
seqkit seq -Q 20 -M 40000 ./input/runs_3567.fastq.gz \
    | gzip > ./input/low_qual_filt_reads/runs_3567.fastq.gz
seqkit seq -Q 20 -M 40000 ./input/runs_3569.fastq.gz \
    | gzip > ./input/low_qual_filt_reads/runs_3569.fastq.gz

# Filter low quality reads from selected reads
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3511_bwa_no_unmap.fastq \
    ./results/runs_3511_bwa_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3511_bwa_unmap.fastq \
    ./results/runs_3511_bwa_gap_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3511_minimap_no_unmap.fastq \
    ./results/runs_3511_minimap_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3511_minimap_unmap.fastq \
    ./results/runs_3511_minimap_gap_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3512_bwa_no_unmap.fastq \
    ./results/runs_3512_bwa_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3512_bwa_unmap.fastq \
    ./results/runs_3512_bwa_gap_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3512_minimap_no_unmap.fastq \
    ./results/runs_3512_minimap_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3512_minimap_unmap.fastq \
    ./results/runs_3512_minimap_gap_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3567_bwa_no_unmap.fastq \
    ./results/runs_3567_bwa_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3567_bwa_unmap.fastq \
    ./results/runs_3567_bwa_gap_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3567_minimap_no_unmap.fastq \
    ./results/runs_3567_minimap_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3567_minimap_unmap.fastq \
    ./results/runs_3567_minimap_gap_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3569_bwa_no_unmap.fastq \
    ./results/runs_3569_bwa_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3569_bwa_unmap.fastq \
    ./results/runs_3569_bwa_gap_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3569_minimap_no_unmap.fastq \
    ./results/runs_3569_minimap_gap_no_unmap_clean.fastq
seqkit seq -Q 20 -M 40000 \
    -o ./results/low_qual_filt_reads/runs_3569_minimap_unmap.fastq \
    ./results/runs_3569_minimap_gap_unmap_clean.fastq

# Run raven
raven -t 6 -F ./results/low_qual_filt_reads/runs_3511_bwa_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3511_bwa_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3511_bwa_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3511_bwa_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3511_bwa_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3511_bwa_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3511_minimap_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3511_minimap_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3511_minimap_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3511_minimap_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3511_minimap_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3511_minimap_unmap.fasta

raven -t 6 -F ./results/low_qual_filt_reads/runs_3512_bwa_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3512_bwa_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3512_bwa_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3512_bwa_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3512_bwa_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3512_bwa_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3512_minimap_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3512_minimap_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3512_minimap_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3512_minimap_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3512_minimap_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3512_minimap_unmap.fasta

raven -t 6 -F ./results/low_qual_filt_reads/runs_3567_bwa_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3567_bwa_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3567_bwa_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3567_bwa_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3567_bwa_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3567_bwa_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3567_minimap_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3567_minimap_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3567_minimap_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3567_minimap_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3567_minimap_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3567_minimap_unmap.fasta

raven -t 6 -F ./results/low_qual_filt_reads/runs_3569_bwa_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3569_bwa_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3569_bwa_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3569_bwa_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3569_bwa_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3569_bwa_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3569_minimap_no_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3569_minimap_no_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3569_minimap_no_unmap.fasta
raven -t 6 -F ./results/low_qual_filt_reads/runs_3569_minimap_unmap.gfa \
    ./results/low_qual_filt_reads/runs_3569_minimap_unmap.fastq \
    > ./results/low_qual_filt_reads/runs_3569_minimap_unmap.fasta

# Turn gfa files into fasta files
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3511_bwa_no_unmap.gfa \
    --output-name runs_3511_bwa_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3511_bwa_unmap.gfa \
    --output-name runs_3511_bwa_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3511_minimap_no_unmap.gfa \
    --output-name runs_3511_minimap_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3511_minimap_unmap.gfa \
    --output-name runs_3511_minimap_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/

python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3512_bwa_no_unmap.gfa \
    --output-name runs_3512_bwa_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3512_bwa_unmap.gfa \
    --output-name runs_3512_bwa_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3512_minimap_no_unmap.gfa \
    --output-name runs_3512_minimap_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3512_minimap_unmap.gfa \
    --output-name runs_3512_minimap_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/

python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3567_bwa_no_unmap.gfa \
    --output-name runs_3567_bwa_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3567_bwa_unmap.gfa \
    --output-name runs_3567_bwa_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3567_minimap_no_unmap.gfa \
    --output-name runs_3567_minimap_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3567_minimap_unmap.gfa \
    --output-name runs_3567_minimap_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/

python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3569_bwa_no_unmap.gfa \
    --output-name runs_3569_bwa_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3569_bwa_unmap.gfa \
    --output-name runs_3569_bwa_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3569_minimap_no_unmap.gfa \
    --output-name runs_3569_minimap_no_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/
python gfa_to_fasta.py \
    --input ./results/low_qual_filt_reads/runs_3569_minimap_unmap.gfa \
    --output-name runs_3569_minimap_unmap_gfa.fasta \
    --output-folder results/low_qual_filt_reads/

# Make a alignments
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3511_bwa_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_bwa_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3511_bwa_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_bwa_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3511_minimap_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_minimap_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3511_minimap_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_minimap_unmap \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3512_bwa_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_bwa_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3512_bwa_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_bwa_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3512_minimap_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_minimap_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3512_minimap_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_minimap_unmap \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3567_bwa_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_bwa_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3567_bwa_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_bwa_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3567_minimap_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_minimap_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3567_minimap_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_minimap_unmap \
    --mafft-threads 4

python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3569_bwa_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_bwa_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3569_bwa_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_bwa_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3569_minimap_no_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_minimap_no_unmap \
    --mafft-threads 4
python align_fastas.py \
    --fasta1 ./results/gap_5k.fasta \
    --fasta2 ./results/low_qual_filt_reads/runs_3569_minimap_unmap_gfa.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_minimap_unmap \
    --mafft-threads 4

# Visualize the alignments

python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3511_bwa_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3511_bwa_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_bwa_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3511_bwa_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3511_bwa_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_bwa_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3511_minimap_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3511_minimap_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_minimap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3511_minimap_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3511_minimap_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3511_minimap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3512_bwa_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3512_bwa_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_bwa_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3512_bwa_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3512_bwa_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_bwa_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3512_minimap_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3512_minimap_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_minimap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3512_minimap_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3512_minimap_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3512_minimap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3567_bwa_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3567_bwa_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_bwa_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3567_bwa_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3567_bwa_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_bwa_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3567_minimap_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3567_minimap_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_minimap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3567_minimap_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3567_minimap_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3567_minimap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138

python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3569_bwa_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3569_bwa_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_bwa_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3569_bwa_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3569_bwa_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_bwa_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3569_minimap_no_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3569_minimap_no_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_minimap_no_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
python plot_alignments.py \
    --blast-tsv ./results/low_qual_filt_reads/alignments/3569_minimap_unmap.blast.tsv \
    --mafft-fasta ./results/low_qual_filt_reads/alignments/3569_minimap_unmap.mafft.fasta \
    --out-prefix ./results/low_qual_filt_reads/alignments/3569_minimap_unmap_plot \
    --anchor-id-substring ChrC_C_glabrata_CBS138
