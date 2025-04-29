#!/bin/bash
<<Step0.1
This script is how I ran trimgalore to remove adaptor sequences from the raw smRNASeq files.
The previous step in the pipeline is: fastqc_raws_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/FastQC/fastqc_raws_2024-5-13.sh)
The next step in the pipeline is: fastqc_trimmed_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/FastQC/fastqc_trimmed_2024-5-13.sh).
Step0.1

# Define variables
    # Fastq files to be trimmed
    fastqs=(
        /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All/LVG_2.fq
        /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All/LVG_4.fq
        /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All/LVG_9.fq
        /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All/RVG_5S.fq
        /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All/RVG_6S.fq
        /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All/RVG_7S.fq
        /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All/RVG_12S.fq
    )

    # Set output directory
    output_directory="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ"
    # Make said directory
    mkdir -p "$output_directory"

# Loop through FASTQs and trim
for fastq in "${fastqs[@]}"; do
    # Define output files
    output_fastq="$output_directory/$(basename "${fastq%.fq}")_trimmed.fq"

    # Run trimgalore
    trim_galore --length 17 --output_dir "$output_fastq" "$fastq"
    # Let the user know when each trimming is done
    echo "Trimming of $fastq completed! Output file $output_fastq complete!"
done
