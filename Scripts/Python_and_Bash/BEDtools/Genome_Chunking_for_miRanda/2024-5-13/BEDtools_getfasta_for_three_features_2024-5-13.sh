#!/bin/bash
<<Step4
The point of this file is to record the commands used to generate the Cvv_2017_genome_with_myos_3UTR.fasta that I then fed to miranda.
The previous step of this pipeline is: Feature_GTF_Generation_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/GTF_generation/2024-5-13/Feature_GTF_Generation_2024-5-13.sh)
The next steps of the pipeline is: miranda_mature_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miRanda_runs/2024-5-13/miranda_mature_2024-5-13.sh).
Step4

# Define the input fasta file
genome_fasta="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos.fasta"

# Define the output directory
output_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/"

# Define the bed files
gtf_files=(
    "CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_CDS.gtf"
    "CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_three_prime_utr.gtf" 
    "CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_five_prime_utr.gtf"
    )

# Define the output fasta files
output_fastas=(
    "Cvv_2017_genome_with_myos_CDS.fasta" 
    "Cvv_2017_genome_with_myos_three_prime_utr.fasta" 
    "Cvv_2017_genome_with_myos_five_prime_utr.fasta"
    )

# Loop through each bed file and execute the corresponding bedtools command
for ((i=0; i<${#gtf_files[@]}; i++)); do
    bedtools getfasta -fo "${output_dir}${output_fastas[$i]}" -fi "$genome_fasta" -bed "${output_dir}${gtf_files[$i]}"
done
