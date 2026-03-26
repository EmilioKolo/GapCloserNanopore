#!/bin/bash

# Make the script stop if any command fails
set -e

python extract_fasta_regions.py \
    --input ./input/C_glabrata_CBS138_current_chromosomes.fasta \
    --region ChrC_C_glabrata_CBS138:95662-100662 \
    --region ChrC_C_glabrata_CBS138:101663-106663 \
    --output-folder input/ \
    --output-name reference_anchors.fasta
