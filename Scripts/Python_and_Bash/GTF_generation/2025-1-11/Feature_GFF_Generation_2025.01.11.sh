#!/bin/bash

<<'Step3-GettingGFFs'
The point of this file is to record the commands used to generate the different GFFs for each sample containing the different 
target types.

Previous step: mature_miRNA_extraction_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/ShortStack_out_to_miRanda_data/mature_miRNA_extraction_2024-5-13.sh)
            and
            BCFtools_liftoff_get_fasta_and_gff_for_all_samples_2025.01.11.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BCFtools/BCFtools_liftoff_get_fasta_and_gff_for_all_samples_2025.01.11.sh)

Next step: BEDtools_getfasta_for_three_features_2025.01.12.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BEDtools/Genome_Chunking_for_miRanda/2025-01-12/BEDtools_getfasta_for_three_features_2025.01.12.sh)
Step3-GettingGFFs

# Find all of the GFF files and create an array from them
# Define the file path to the directory containing them
gff_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/gff_files"

# Define the output directory
output_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/gff_files/Sub_features"

# Make the directory if it doesn't exist yet
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"

# Define the types of features
features=("three_prime_utr" "five_prime_utr" "CDS")

# Find GFF files
sample_gffs=()
while IFS= read -r gff; do
    sample_gffs+=("$gff")
done < <(find "$gff_dir" -type f -name "*.gff")

# Add the reference genome GFF
ref_gff="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18.gff"

# If else statement that adds the reference to the list 
if [ -f "$ref_gff" ]; then
    sample_gffs+=("$ref_gff")
    echo "Added additional GFF: $ref_gff"
else
    echo "Additional GFF not found: $ref_gff"
fi

# Ensure files are found
if [ ${#sample_gffs[@]} -eq 0 ]; then
    echo "No GFF files found in $gff_dir or additional file."
    exit 1
fi

# Run a loop that finds each GFF and 
for gff in "${sample_gffs[@]}"; do
    echo "Processing GFF: $gff"

    # Get base filename
    base_file_name=$(basename "$gff" .gff)

    # Loop through each feature and execute the corresponding awk command
    for feature in "${features[@]}"; do

        # Set output file
        output_file="${output_dir}/${base_file_name}_${feature}.gff"

        # Run awk
        if awk -v feature="$feature" '$3 == feature' "$gff" > "$output_file"; then
            echo "Extracted $feature to $output_file"
        else
            echo "Error processing $feature for $gff"
        fi
    done
done