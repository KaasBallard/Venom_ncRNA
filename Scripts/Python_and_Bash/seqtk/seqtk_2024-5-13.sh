#/bin/bash
<<Step0.3
This script is how I turned fastqs to fasta files after trimming.
The previous step in the pipeline is: fastqc_trimmed_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/FastQC/fastqc_trimmed_2024-5-13.sh)
The next step in the pipeline is: ShortStack_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/ShortStack/ShortStack_2024-5-13.sh).
Step0.3

# This program converts fastqs to fasta.
# Var definitions
fastqs=(
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All/LVG_2_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All/LVG_4_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All/LVG_9_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All/RVG_5S_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All/RVG_6S_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All/RVG_7S_trimmed.fq
    /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All/RVG_12S_trimmed.fq
)

output_location="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTA"

# Run seqtk:
for fastq in "${fastqs[@]}"; do
    # Loop that extracts basename of the file iteratively:
    fasta="${output_location}/$(basename "${fastq%.fq}.fasta")"
    echo "Name conversion complete for fastqs file $fastq"

    # Run seqtk:
    seqtk seq -a "${fastq}" > "${fasta}"
    echo "fastqs to fasta conversion complete for $fastq"
done
