#!/bin/bash

# Define variables
reference_genome="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos.fasta"
genome_index="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos.fasta"
fastqs=(
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/All/LVG_2_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/All/LVG_4_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/All/LVG_9_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/All/RVG_5S_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/All/RVG_6S_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/All/RVG_7S_trimmed.fq
)
output_directory="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_SAMs/"
mkdir -p "$output_directory"

# Create Bowtie index if it doesn't exist
if [ ! -f "${genome_index}.1.ebwt" ]; then
    bowtie-build $reference_genome $genome_index
fi

# Loop through FASTQs and align with Bowtie
for fastq in "${fastqs[@]}"; do
    # Define output files
    output_sam="$output_directory/$(basename "${fastq%.fq}")_aligned.sam"

    # Run Bowtie for single-end alignment
    bowtie -q $genome_index $fastq -S $output_sam
    
    # Let the user know when each alignment is done
    echo "Alignment of $fastq completed! Output: $output_sam"
done
