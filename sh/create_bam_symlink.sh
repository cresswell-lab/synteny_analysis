#!/usr/bin/env bash

# elements=$(ls /nobackup/lab_cresswell/vgazziero/shallow_sarek/results/results/preprocessing/recalibrated)

# for e in elements; do
#     bams=$(ls e)
#     echo $bams 
# done

# This script creates soft links for BAM files in a specified directory.
 
# Define the source directory containing BAM files
SOURCE_DIR="/nobackup/lab_cresswell/vgazziero/synteny_prj/cut_and_tag_processing/run_and_tag_alignment/results/02_alignment/bowtie2/target/markdup"
 
# Define the target directory where soft links will be created
TARGET_DIR="/nobackup/lab_cresswell/vgazziero/synteny_prj/cut_and_tag_processing/cna_calling"
 
# Get list of BAM files in the source directory
mapfile -t BAMS < <(find "$SOURCE_DIR" -type f -name "*.bam")
 
# Show what is in BAMS
echo "BAM files found in source directory:"
for BAM in "${BAMS[@]}"; do
    echo "$BAM"
done
 
# Prompt user for confirmation before proceeding
# read -p "Press enter to continue or Ctrl+C to cancel..."
 
# Loop through each BAM file and create a soft link in the target directory
for BAM in "${BAMS[@]}"; do
    # Extract the base name of the BAM file
    BASENAME=$(basename "$BAM")
 
    # Create the soft link in the target directory
    ln -s "$BAM" "$TARGET_DIR/$BASENAME"
 
    echo "Created soft link for $BASENAME in $TARGET_DIR"
done
 
# Notify user that the process is complete
echo "Soft links created successfully in $TARGET_DIR."

# Get list of BAM files in the source directory
mapfile -t BAIS < <(find "$SOURCE_DIR" -type f -name "*.bai")
 
# Show what is in BAMS
echo "BAI files found in source directory:"
for BAI in "${BAIS[@]}"; do
    echo "$BAI"
done
 
# Prompt user for confirmation before proceeding
# read -p "Press enter to continue or Ctrl+C to cancel..."
 
# Loop through each BAM file and create a soft link in the target directory
for BAI in "${BAIS[@]}"; do
    # Extract the base name of the BAM file
    BASENAME=$(basename "$BAI")
 
    # Create the soft link in the target directory
    ln -s "$BAI" "$TARGET_DIR/$BASENAME"
 
    echo "Created soft link for $BASENAME in $TARGET_DIR"
done
 
# Notify user that the process is complete
echo "Soft links created successfully in $TARGET_DIR."

