#/bin/bash

<<Step6
This script is how I created tabular outputs from miRanda for all 3'UTR, 5'UTR, and CDS.
The previous step in the pipeline is: miranda_mature_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miRanda_runs/2024-3-28/miranda_mature_2024-5-13.sh)
The next step in the pipeline is: miranda_output_conversion_BED_2024-5-13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miranda_formating/miranda_output_conversion_BED_2024-5-13.ipynb)
Step6

<<Source
I found this site: https://bioinformaticsworkbook.org/dataAnalysis/SmallRNA/Miranda_miRNA_Target_Prediction.html#gsc.tab=0 and used it to help me create this.
Source

# Files to be reformated
miranda_outs=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2024-5-13/Cvv_2017_genome_with_myos_CDS_miranda_miRNA_targets.out"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2024-5-13/Cvv_2017_genome_with_myos_five_prime_utr_miranda_miRNA_targets.out"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2024-5-13/Cvv_2017_genome_with_myos_three_prime_utr_miranda_miRNA_targets.out"
    )
    
# For loop that reformates the files
for file in "${miranda_outs[@]}"; do
    cat "$file" | grep ">>" | sort -k5,5nr | sed 's/>>//g' | cat <(echo "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions" | tr "," "\t") - > "${file%.out}.tab"
done