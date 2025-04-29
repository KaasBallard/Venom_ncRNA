#/bin/bash

<<'Step6-Convert_miRanda_outs_to_tab'
This script is how I created tabular outputs from miRanda for all 3'UTR, 5'UTR, and CDS for all of the samples.
Previous step: miranda_2025.01.12.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miRanda_runs/2025-01-12/miranda_2025.01.12.sh)
Next step: miranda_output_conversion_BED_2025.01.13.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miranda_formating/miranda_output_conversion_BED_2025.01.13.ipynb)
Step6-Convert_miRanda_outs_to_tab

<<Source
I found this site: https://bioinformaticsworkbook.org/dataAnalysis/SmallRNA/Miranda_miRNA_Target_Prediction.html#gsc.tab=0 and used it to help me create this.
Source

# Set the file path containing the miRanda outs
miranda_dir="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_2025-01-12"

# Find miRanda outputs and add them to an array
miranda_outs=()
while IFS= read -r out; do
    miranda_outs+=("$out")
done < <(find "$miranda_dir" -type f -name "*.out")

# Check to make sure the array is not empty
if [ ${#miranda_outs[@]} -eq 0 ]; then
    echo "Error: No miRanda output files found in $miranda_outs"
    exit 1
else
    echo "Success: miRanda output files found. The following will be converted to tabular format:"
    printf '%s\n' "${miranda_outs[@]}"
fi

  
# For loop that reformates the files
for file in "${miranda_outs[@]}"; do
    
    # Create path for new tabular file
    base_name=$(basename "$file" .out)
    tab_file="${miranda_dir}/${base_name}.tab"

    # Check if the output file exists
    if [ ! -f "$tab_file" ]; then
        echo "$tab_file does not exist, converting $file to tabular output."

        # Format tabular file
        cat "$file" | grep ">>" | sort -k5,5nr | sed 's/>>//g' | cat <(echo "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions" | tr "," "\t") - > "$tab_file"

        echo "Success: Tabular output created: $tab_file"
    else
        echo "miRanda tabular output already exists. Skipping."
    fi
done