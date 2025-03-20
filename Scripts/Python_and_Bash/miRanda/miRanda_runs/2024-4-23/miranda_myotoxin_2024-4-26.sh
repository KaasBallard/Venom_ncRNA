#/bin/bash

<<myotoxin
This script is how I ran miRanda against the myotoxin sequences to see if any thing finds. There are no next steps as this is an asside.
myotoxin

# Define variables:

    # mirna_file or the RNA sequence to be compared to genomes:
    mirna_file="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/miRanda_mirna_inputs_from_shortstack/2024-4-9_Run/Post_clean/mature_mir.fasta"

    # genomes or the genomic DNA/RNA sequence, alignment procedure uses complementarity and not sequence identity:
    genomes=(
        "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/myotoxin_lncRNA/myotoxin_lncRNA.fasta"
    )

    # Name of location for the output file:
    output_directory="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_myotoxin_2024-4-26"

# Run the program:
    # Actually run the program:
    for genome in "${genomes[@]}"; do
        # Extract the genome name without the path for the directory
        genome_name=$(basename "$genome" .fasta)

        # Create new file name for the output
        output_file="${output_directory}/${genome_name}_miranda_miRNA_targets.out"

        # Run miRanda
        miranda "$mirna_file" "$genome" -out "$output_file"

        # Report to me that the output has been complete
        echo "miRanda run for $genome_name complete."
    done