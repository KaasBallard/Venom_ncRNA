#!/bin/bash

<<Step8
This script is meant to run the bedtools intesect command in a loop so that the miRanda 
formatted outputs (in .bed format) can be compared to the genome.
The previous steps in the pipeline is: miranda_output_conversion_BED_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miranda_formating/miranda_output_conversion_BED_2024-5-13.ipynb)
The next step in the pipeline is: 3UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/3UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb)
                                  5UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/5UTR_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb)
                                  exons_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/exons_miranda.tab_bedtools_intersect_fusion_2024-5-13.ipynb)
Step8

# Create array for bedtools
inputs=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2024-5-13/Cvv_2017_genome_with_myos_CDS_miranda_miRNA_targets.bed"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2024-5-13/Cvv_2017_genome_with_myos_five_prime_utr_miranda_miRNA_targets.bed"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2024-5-13/Cvv_2017_genome_with_myos_three_prime_utr_miranda_miRNA_targets.bed"
    )

genomes=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_CDS.gtf"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_five_prime_utr.gtf"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_three_prime_utr.gtf"
    )

output_prefix="bedtools_intersect"

for ((i=0; i<${#inputs[@]}; i++)); do
    # Initialize variables for the loop
    input="${inputs[$i]}"
    genome="${genomes[$i]}"

    input_dir=$(dirname "$input")
    
    out_file="${input_dir}/${output_prefix}_$(basename "${input%.*}").bed"

    # Using the -wa and -wb flags so that the entries are written right next to each other
    bedtools intersect -wa -wb -a "$genome" -b "$input" > "$out_file"

    echo "The intersection between the genome: $genome and $input has been completed"
done