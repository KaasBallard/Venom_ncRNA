#!/bin/bash

# Define variables
    # Define variable for the SAM files
    sam_files=(
        "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_SAMs/RVG_5S_trimmed_aligned.sam"
        "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_SAMs/RVG_6S_trimmed_aligned.sam"
        "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_SAMs/RVG_7S_trimmed_aligned.sam"
        "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_SAMs/LVG_2_trimmed_aligned.sam"
        "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_SAMs/LVG_4_trimmed_aligned.sam"
        "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_SAMs/LVG_9_trimmed_aligned.sam"
    )

    # Define a new path for the BAM files 
    bam_directory="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_BAMs"

    # Create the directory if it doesn't exist
    mkdir -p "$bam_directory"

# Loop through the directory and make the BAMs
for sam_file in "${sam_files[@]}"; do

    # Extract the basename to give the BAMs names
    filename=$(basename "${sam_file}")
    output_bam="${bam_directory}/${filename%.sam}.bam"

    # Convert SAM to BAM using Samtools
    samtools view -bS "$sam_file" -o "$output_bam"

    # Index the BAM
    samtools index "$output_bam"

    # Echo that the BAM file is created
    echo "SAM to BAM conversion and indexing complete for: $output_bam"
done