#/bin/bash

<<Source
I found this site: https://bioinformaticsworkbook.org/dataAnalysis/SmallRNA/Miranda_miRNA_Target_Prediction.html#gsc.tab=0 and used it to help me create this.
Source

# For loop for the files
miranda_outs=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_5S_trim_mir"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_6S_trim_mir"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_7S_trim_mir"
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_12S_trim_mir"
)

for file in "${miranda_outs[@]}"; do
    cat "$file" | grep ">>" | sort -k5,5nr | sed 's/>>//g' | cat <(echo "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions" | tr "," "\t") - > "${file%.out}_miranda_output.tab"
done