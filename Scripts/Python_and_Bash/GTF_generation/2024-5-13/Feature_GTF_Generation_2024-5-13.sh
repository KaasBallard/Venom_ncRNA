#!/bin/bash

<<Step3
The point of this file is to record the commands used to generate the different CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf's.
The previous steps of this pipeline is: mature_miRNA_extraction_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/ShortStack_out_to_miRanda_data/mature_miRNA_extraction_2024-5-13.sh)
The next steps of the pipeline is: BEDtools_getfasta_for_three_features_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BEDtools/Genome_Chunking_for_miRanda/2024-5-13/BEDtools_getfasta_for_three_features_2024-5-13.sh)
Step3

# Define the file path
input_file="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP.gtf"

# Define the output directory
output_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/"

# Define the types of features
features=("three_prime_utr" "five_prime_utr" "CDS")

# Loop through each feature and execute the corresponding awk command
for feature in "${features[@]}"; do
    awk -v feature="$feature" '$3 == feature' "$input_file" > "${output_dir}CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_${feature}.gtf"
done
