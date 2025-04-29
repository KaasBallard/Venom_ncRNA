#/bin/bash

<<Step0.2
This is the first step I did for this project. This script is for recording how I checked the data quality for the various read files provided by Novogene.
The previous step in the pipeline is: trimgalore_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/trimgalore/trimgalore_2024-5-13.sh)
The next step in the pipeline is: seqtk_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/seqtk/seqtk_2024-5-13.sh).
Step0.2

# Create an array contain the directories I want to loop through
directories=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_reads/FASTQ/All"
)

# Set the number of threads
threads=32

# Run FastQC in a loop
for dir in "${directories[@]}"; do
    echo "Checking read quality for sRNA reads in the directory: $dir"

    # Run fastqc
    fastqc "$dir"/*.fq -t "$threads"
done


<<Note
I moved the resulting fastqc outputs to the directory: /home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/fastqc/Trimmed_Quality_2024-5-13
Note