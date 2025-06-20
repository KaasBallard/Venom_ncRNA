#/bin/bash

<<Step0
This is the first step I did for this project. This script is for recording how I checked the data quality for the various read files provided by Novogene.
The next step in the pipeline is: trimgalore_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/trimgalore/trimgalore_2024-5-13.sh).
Step0

# Create an array contain the directories I want to loop through
directories=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Raw_sRNA_reads/All"
)

# Set the number of threads
threads=32

# Run FastQC in a loop
for dir in "${directories[@]}"; do
    echo "Checking read quality for sRNA reads in the directory: $dir"

    # Run fastqc
    fastqc "$dir"/*.fq.gz -t "$threads"
done


<<Note
I moved the resulting fastqc outputs to the directory: /home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/fastqc/Raw_Quality_2024-5-13
Note