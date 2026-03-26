#!/bin/bash

# Make the script stop if any command fails
set -e

# Extract the reads that map to the defined gap region
python extract_reads.py --input input/runs_3511_minimap_sorted.bam \
    --region ChrC_C_glabrata_CBS138:90000-150000 \
    --output-folder results/ \
    --output-name runs_3511_minimap_gap_unmap.fastq.gz \
    --include-unmapped
