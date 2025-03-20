#!/bin/bash
<<'Step4-Get_three_five_utr_and_CDS_FASTAs'
The point of this file is to record the commands used to generate the FASTA files for miRanda to be fed.

Previous step: Feature_GFF_Generation_2025.01.11.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/GTF_generation/2025-1-11/Feature_GFF_Generation_2025.01.11.sh)
Next step: miranda_mature_2025.01.12.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miRanda_runs/2025-01-12/miranda_mature_2025.01.12.sh).
Step4-Get_three_five_utr_and_CDS_FASTAs

# Activate the needed conda/mamba environment
source /home/administrator/miniforge3/bin/activate bed-sam-bcf-tools

# Define the reference genome FASTA file
ref_genome_fasta="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos.fasta"

# Define the output directory
output_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/feature_files"
[ ! -d "$output_dir" ] && mkdir -p "$output_dir"

# Set the directories for GFF and FASTA files
gff_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/gff_files/Sub_features"
fasta_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files"

# Add reference genome to FASTA list
input_fasta_files=()
while IFS= read -r fasta; do
    input_fasta_files+=("$fasta")
done < <(find "$fasta_dir" -type f -name "*.fasta")

# If else statement that adds the reference to the list 
if [ -f "$ref_genome_fasta" ]; then
    input_fasta_files+=("$ref_genome_fasta")
    echo "Added reference FASTA: $ref_genome_fasta"
else
    echo "Reference FASTA not found: $ref_genome_fasta"
fi

echo "The following FASTA files will be processed:"
printf '%s\n' "${input_fasta_files[@]}"

# Process each FASTA file and match it with corresponding GFF files
for input_fasta in "${input_fasta_files[@]}"; do

    # Get base filename
    base_name=$(basename "$input_fasta" .fasta)
    
    # If the base name doesn't match exactly, manually adjust it
    if [[ "$base_name" == "Cvv_2017_genome_with_myos" ]]; then
        base_name="Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18"
    fi

    echo "Processing FASTA: $input_fasta"

    # Find matching GFF files
    matching_gffs=$(find "$gff_dir" -type f -name "${base_name}_*.gff")
    if [ -z "$matching_gffs" ]; then
        echo "No matching GFF files found for $base_name"
        continue
    fi

    # Process each matching GFF
    while IFS= read -r gff; do
        gff_suffix=$(basename "$gff" | sed -E "s/^${base_name}_//;s/\.gff$//")

        # Set the output file name
        output_fasta="${output_dir}/${base_name}_${gff_suffix}.fasta"
        echo "Generating FASTA for $gff -> $output_fasta"

        # Run bedtools
        bedtools getfasta -fi "$input_fasta" \
            -bed "$gff" \
            -fo "$output_fasta"

        # Check if the output file was created
        if [ ! -f "$output_fasta" ]; then
            echo "Error: FASTA file: $output_fasta not created"
        else
            echo "Success: FASTA file: $output_fasta created"
        fi
    done <<< "$matching_gffs"
done

echo "FASTA generation completed."

<<Name_change_note
Note that the names for:
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/temp/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18_CDS.fasta
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/temp/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18_five_prime_utr.fasta
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/temp/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18_three_prime_utr.fasta

/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/gff_files/Sub_features/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18_CDS.gff
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/gff_files/Sub_features/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18_five_prime_utr.gff
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/gff_files/Sub_features/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18_three_prime_utr.gff


where changed to:
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/feature_files/Crotalus_viridis_reference_CDS.fasta
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/feature_files/Crotalus_viridis_reference_five_prime_utr.fasta
/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/feature_files/Crotalus_viridis_reference_three_prime_utr.fasta



Because, for some reason miRanda experiences a error that looks like this:
*** stack smashing detected ***: terminated
miranda_mature_2025.01.12.sh: line 39: 24376 Aborted                 (core dumped) miranda "$mirna_file" "$genome" -out "$output_file"

Basically the buffer overflows because the file name is to long. Incredible.
Name_change_note
