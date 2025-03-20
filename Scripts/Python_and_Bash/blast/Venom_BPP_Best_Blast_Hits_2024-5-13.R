#### Step 0.5-Finding_BPP ####
# The point of this file is to find the best blast hits for BPP and print them.
# The previous step in the pipeline is: blast_BPP_alignment_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/blast/blast_BPP_alignment_2024-5-13.sh)
# After this I went into AliView and for all of the top sequences I found in the out put file. From AliView I took the longest sequence and used the start and stop of the BPP transcript model
# and as the coding region. I then took anything before that as the 5'UTR and anything after that as the 3'UTR.
# I then and added the 3'UTR and 5'UTR seqences to the gtf files for the 3UTR and 5UTR (it was already in the CDS fasta). I then calculated the start and end of the 3'UTR and 5'UTR using the next step.
# The next step in the pipeline is: (3_and_5_prime_UTR_Calculations.py) /home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/Python/3_and_5_prime_UTR_Calculations.py
library(tidyverse)

# Filter Blast data for the best hits
blast_df <- read.table('/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/blast/Venom_BPP/Venom_BPP_Blast_Hits_2024-5-13.txt') %>% 
  rename(
    Q_id = V1,
    S_id = V2,
    percent_id = V3,
    alignment_length = V4,
    mismatches = V5,
    gap_openings = V6,
    q_start = V7,
    q_end = V8,
    s_start = V9,
    s_end = V10,
    e_value = V11,
    bit_score = V12
  ) %>% 
  slice_max(bit_score) %>% 
  slice_min(e_value)  %>% 
  slice_max(percent_id)

# Save the column to file
write.csv(blast_df$S_id, file = '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/blast/Venom_BPP/BPP_Best_Blast_Hits_2024-5-13.csv', quote = F, row.names = F)
