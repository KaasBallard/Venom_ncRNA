#!/bin/bash

<<Step1
This script is how I ran the second step of the pipeline. Here I simply ran ShortStack for every
.fasta file in the Clean_sRNA_reads directory. The .fasta files were generated from the clean data
by unziping all of the .fa.gz files and running gzip -dk *.gz. For some reason ShortStack didn't
like the .fa files because they start with ">@", so I removed the "@" with the following sed command:
sed 's/^>@/>/' .fa > .fasta.
I did this for all of the .fa files in the directory. 
The previous step in the pipeline is: 3_and_5_prime_UTR_Calculations.py (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/3_and_5_prime_UTR_Calculations.py)
The next step in the pipeline is: mature_miRNA_extraction_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/ShortStack_out_to_miRanda_data/mature_miRNA_extraction_2024-4-9.sh).
Step1

# Define variables:

    # Define genomefile path:
    genome_file="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos.fasta"

    # Create a variable that has the path for the directory that contains all of the necessary readfiles
    read_file_directory="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTA"

    # Output directory:
    output_directory="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/ShortStack/ShortStack_results_2024-5-13"

    # Number of threads:
    threads=32

    mirbase="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/miRBase/miRBase_database.fa"

# Run the program:
ShortStack --genomefile "${genome_file}" \
            --readfile "${read_file_directory}"/*.fasta \
            --outdir "${output_directory}" \
            --mmap u \
            --known_miRNAs "${mirbase}" \
            --dn_mirna \
            --threads "${threads}" 