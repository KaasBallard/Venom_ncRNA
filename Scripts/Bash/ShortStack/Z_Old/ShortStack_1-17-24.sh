#!/bin/bash

# Define variables:

    # Define genomefile path:
    genome_file="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos.fasta"

    # Define readfile paths:
    read_file_1="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/short_RNA_reads_fastq/RVG_5S_trim.fq" 
    read_file_2="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/short_RNA_reads_fastq/RVG_6S_trim.fq"
    read_file_3="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/short_RNA_reads_fastq/RVG_7S_trim.fq"  
    read_file_4="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/short_RNA_reads_fastq/RVG_12S_trim.fq" 

    # Output directory:
    output_directory="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/ShortStack/1-17-24_ShortStack_results_readfiles"

    # Number of threads:
    threads=32

# Run the program:
ShortStack --genomefile "${genome_file}" \
            --readfile "${read_file_1}" "${read_file_2}" "${read_file_3}" "${read_file_4}" \
            --outdir "${output_directory}" \
            --mmap u \
            --dn_mirna \
            --threads "${threads}" 