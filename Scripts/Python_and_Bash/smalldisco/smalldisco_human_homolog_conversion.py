#!/usr/bin/env python3

'''The point of this file is to convert all the sirna_name's in the smalldisco outputs to human homologs found in the converted names file.'''

# Import pandas and os
import pandas as pd
import os

# Create a dataframe for the .txt
converter_df = pd.read_csv('/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_GTF_to_converted_names_05.22.23.txt', sep='\t')

# All of the bed files to be converted
bed_files = [
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/smalldisco/sirna_11-20-23/RVG_5S/RVG_5S_trim.bed',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/smalldisco/sirna_11-20-23/RVG_6S/RVG_6S_trim.bed',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/smalldisco/sirna_11-20-23/RVG_7S/RVG_7S_trim.bed',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/smalldisco/sirna_11-20-23/RVG_12S/RVG_12S_trim.bed'
]

# Loop through files to be converted
for bed in bed_files:
    # Data fram for each file that gets looped through:
    sirna_df = pd.read_csv(bed, sep='\t', comment='#', header=None, names=['chrom', 'start', 'end', 'sirna_name', 'score', 'strand', 'num_reads'], skiprows=1)

    # Function to extract the matching entry from the converter_df
    def Get_Matching_siRNA (sirna_name):
        # Extract the unique part by removing the prefix and suffix
        shared_entry = sirna_name.replace('sirna_', '').replace('_01', '')
    
        # Find the matching entry in the converter_df
        match = converter_df[converter_df['gtf_gene'].str.contains(shared_entry)]['converted_id_no_dups']
    
        # Return the first matching entry, or None if not found
        return match.iloc[0] if not match.empty else None
    
    # Apply the function to replace sirna_name with converted_id_no_dups
    sirna_df.insert(4, 'converted_id_no_dups', sirna_df['sirna_name'].apply(Get_Matching_siRNA))

    # Save the converted .bed to a new file in the same directory.
    out_path = os.path.join(os.path.dirname(bed), os.path.basename(bed).replace('.bed', '_converted.bed'))
    sirna_df.to_csv(out_path, sep='\t', header=True, index=False)
