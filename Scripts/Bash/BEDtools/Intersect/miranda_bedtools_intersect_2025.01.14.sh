#!/bin/bash

<<'Step8-Bedtools_intersect'
This script is meant to run the bedtools intesect command in a loop so that the miRanda 
formatted outputs in BED format can be compared to the genome.

Previous step: miranda_output_conversion_BED_2025.01.13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miranda_formating/miranda_output_conversion_BED_2025.01.13.ipynb)
The next step in the pipeline is: 3UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/3UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb)
                                  5UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/5UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb)
                                  exons_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/exons_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb)
Step8-Bedtools_intersect

# Define log file for capturing both stdout and stderr
log_file="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BEDtools/Intersect/miranda_bedtools_intersect_2025.01.14.log"
exec > >(tee -a "$log_file") 2>&1

# Activate the needed conda/mamba environment
source /home/administrator/miniforge3/bin/activate bed-sam-bcf-tools

# Set an output directory and make it
output_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2025-01-12/bedtools_intersect"
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"

# Define the directory for the input BED files that will be intersected with bedtools
input_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2025-01-12"
echo "The following BED containing directory set"
printf '%s\n' "${input_dir}"

# Set an array for the BED in that dir
input_beds=()
while IFS= read -r bed; do
    input_beds+=("$bed")
done < <(find "$input_dir" -type f -name "*.bed")
echo "The following BED files will be processed:"
printf '%s\n' "${input_beds[@]}"


# Set the directory that the genomes are contained
genome_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/gff_files/Sub_features"
echo "The following GFF containing directory set:"
printf '%s\n' "${genome_dir}"


# Process each BED file
for bed in "${input_beds[@]}"; do
    echo "Processing BED: $bed"

    # Get the base filename
    base_name=$(basename "$bed" _miranda_miRNA_targets.bed)

    # Find matching GFF files
    matching_genomes=$(find "$genome_dir" -type f -name "${base_name}.gff")
    if [ -z "$matching_genomes" ]; then
        echo "No matching genomic GFF files found for $base_name"
        continue
    fi

    echo "Matching GFF files:"
    printf '%s\n' "${matching_genomes}"

    # Process each matching GFF file
    for gff in $matching_genomes; do
        
        # Set output file name
        output_bed="${output_dir}/${base_name}_bedtools_intersect.bed"

        echo "Running bedtools intersect for $bed and $gff -> $output_bed"
        if ! bedtools intersect -wa -wb -a "$gff" -b "$bed" > "$output_bed"; then
            echo "Error: bedtools intersect failed for $bed and $gff"
            continue
        else
            echo "Success: bedtools intersect succeeded for $bed and $gff"
        fi
    done
done