#/bin/bash

<<Step2
This is the 3rd step of the pipeline. This script removes all of the hairpin and the immature
sequences from the mir.fasta file generated by ShortStack.
The previous step in the pipeline is: ShortStack_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/ShortStack/ShortStack_2024-5-13.sh)
The next step in the pipeline is: Feature_GTF_Generation_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/GTF_generation/2024-5-13/Feature_GTF_Generation_2024-5-13.sh)
Step2

# Run the command to extract mature miRNA sequences
shortstack_output="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/miRanda_mirna_inputs_from_shortstack/2024-5-13_Run/Pre_clean/mir.fasta"
   
out_file="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/miRanda_mirna_inputs_from_shortstack/2024-5-13_Run/Post_clean/mature_mir.fasta"

grep -A 1 ".mature" $shortstack_output | grep -v "^--$" > "$out_file"

# Run the command to extract the hairpin sequences:
out_file_haripin="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/miRanda_mirna_inputs_from_shortstack/2024-5-13_Run/Post_clean/haripin_mir.fasta"
awk '
  BEGIN {print_sequence = 0}
  /^>/ { 
    if ($0 ~ /.mature/ || $0 ~ /.star/) 
      print_sequence = 0; 
    else 
      print_sequence = 1;
  }
  print_sequence { print }
' "$shortstack_output" > "$out_file_haripin"