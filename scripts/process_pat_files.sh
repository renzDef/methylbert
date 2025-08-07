#!/bin/bash

# Usage: ./process_pat_files.sh /path/to/pat_dir /path/to/output_dir
# Creates Bam files for all pat files in the given directory using the pat_to_sam script

set -euo pipefail

# Input arguments
PAT_DIR="$1"
OUTPUT_DIR="$2"
TMP_DIR="tmp_processing_dir"
REFERENCE="../data/reference/hg38.fa"
CPG_CITES="../data/reference/cpg_sites.pkl"
PAT_TO_SAM_SCRIPT="pat_to_sam.py"

# Create tmp and output dirs if they don't exist
mkdir -p "$TMP_DIR"
mkdir -p "$OUTPUT_DIR"

# Loop through all .pat.gz files in the PAT_DIR
for PAT_GZ_FILE in "$PAT_DIR"/*.pat.gz; do
    BASENAME=$(basename "$PAT_GZ_FILE" .pat.gz)
    BAM_FILE="$OUTPUT_DIR/$BASENAME.bam"

    # Step 1: Skip if BAM file already exists
    if [[ -f "$BAM_FILE" ]]; then
        echo "Skipping $BASENAME: BAM already exists."
        continue
    fi

    echo "Processing $BASENAME..."

    # Step 2: Unzip to tmp dir, keep original
    UNZIPPED_PAT="$TMP_DIR/$BASENAME.pat"
    gunzip -c "$PAT_GZ_FILE" > "$UNZIPPED_PAT"

    # Step 3: Run pat_to_sam.py to generate SAM file
    SAM_FILE="$TMP_DIR/$BASENAME.sam"
    python3 "$PAT_TO_SAM_SCRIPT" \
         "$UNZIPPED_PAT" \
         "$CPG_CITES" \
         "$REFERENCE" \
        -o "$SAM_FILE" \


    # Step 4: Convert SAM to BAM
    samtools view -b -o "$BAM_FILE" "$SAM_FILE"

    # Step 5: Clean up temp files
    rm -f "$UNZIPPED_PAT" "$SAM_FILE"
    echo "Finished $BASENAME."

done

# Optionally remove tmp dir if empty
rmdir "$TMP_DIR" 2>/dev/null || true
