#!/bin/bash

# Define variables:

    # Define genomefile path:
    genome_file="/home/administrator/Desktop/ExtraSSD2/Kaas/Venom_ncRNA_project/Data/Genome_files/Cvv_2017_genome_with_myos.fasta"

    # Define readfile paths:
    read_files=(
        "/home/administrator/Desktop/ExtraSSD2/Kaas/Venom_ncRNA_project/Data/short_RNA_reads_fastq/RVG_5S_trim.fq"
        "/home/administrator/Desktop/ExtraSSD2/Kaas/Venom_ncRNA_project/Data/short_RNA_reads_fastq/RVG_6S_trim.fq"
        "/home/administrator/Desktop/ExtraSSD2/Kaas/Venom_ncRNA_project/Data/short_RNA_reads_fastq/RVG_7S_trim.fq"
        "/home/administrator/Desktop/ExtraSSD2/Kaas/Venom_ncRNA_project/Data/short_RNA_reads_fastq/RVG_12S_trim.fq"
    )

    # Output directory:
    output_location="/home/administrator/Desktop/ExtraSSD2/Kaas/Venom_ncRNA_project/ShortStack_results/9-1-23_ShortStack_results_readfiles"

    # Number of threads:
    threads=4


# This section runs the the program for each readfile:
    for read_file in "${read_files[@]}"
    do
        # Extract output location so that the loop doesn't overwrite direcotry and creates new ones for each output instead:
        output_directory="${output_location}/$(basename "${read_file}")"

        # Code for running the program:
        ShortStack --genomefile "${genome_file}" \
            --readfile "${read_file}" \
            --outdir "${output_directory}" \
            --mmap u \
            --dn_mirna \
            --threads "${threads}" 

        echo "Search complete for read $readfile"
    done