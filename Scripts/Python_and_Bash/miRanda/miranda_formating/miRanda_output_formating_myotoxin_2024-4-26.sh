#/bin/bash

<<myotoxin
This script is how I created tabular outputs from miRanda for the miRanda myotoxin run.
myotoxin

<<Source
I found this site: https://bioinformaticsworkbook.org/dataAnalysis/SmallRNA/Miranda_miRNA_Target_Prediction.html#gsc.tab=0 and used it to help me create this.
Source

# Files to be reformated
miranda_outs=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_myotoxin_2024-4-26/myotoxin_lncRNA_miranda_miRNA_targets.out"
    )
    
# For loop that reformates the files
for file in "${miranda_outs[@]}"; do
    cat "$file" | grep ">>" | sort -k5,5nr | sed 's/>>//g' | cat <(echo "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions" | tr "," "\t") - > "${file%.out}.tab"
done