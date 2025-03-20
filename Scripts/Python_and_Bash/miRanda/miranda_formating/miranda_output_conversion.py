#!/usr/bin/env python3

import pandas as pd
import os

'''
The point of this file is to convert miranda's tab formated outputs (converted from the miRanda_output_formating_top50K.sh and miRanda_output_formating.sh) into .bed files.
'''

# Create a array of for the tab files to be converted to a .bed
tab_files = [
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_5S_trim_mir_miranda_top50k_output.tab',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_5S_trim_mir_miranda_output.tab',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_6S_trim_mir_miranda_top50k_output.tab',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_6S_trim_mir_miranda_output.tab',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_7S_trim_mir_miranda_top50k_output.tab',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_7S_trim_mir_miranda_output.tab',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_12S_trim_mir_miranda_top50k_output.tab',
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/miRanda/miRanda_mature_11-18-23/RVG_12S_trim_mir_miranda_output.tab'
]

# Create a loop that converts all of the files.
for tab_file in tab_files:
    # Create a data frame to read the miranda tab-formatted output
    mir_df = pd.read_csv(tab_file, sep='\t')

    # Create variables for storing chrom info
    mir_df['chrom'] = mir_df['Seq2'].apply(lambda x: x.split(':')[-2].split('.')[-1])
    mir_df['chromStart'] = mir_df['Seq2'].apply(lambda x: int(x.split(':')[-1].split('-')[0]))
    mir_df['chromEnd'] = mir_df['Seq2'].apply(lambda x: int(x.split(':')[-1].split('-')[1].split('(')[0]))

    # Print intermediate results for troubleshooting
    #print(mir_df[['Seq1', 'chrom', 'chromStart', 'chromEnd']])

    # Create a new data frame with only the updated names
    bed_mir_df = mir_df[['chrom', 'chromStart', 'chromEnd']]

    # Create an output path to the directory that these new files should be put in.
    out_path = os.path.join(os.path.dirname(tab_file), os.path.basename(tab_file).replace('.tab', '.bed'))
    # Save the data frame to a new file
    bed_mir_df.to_csv(out_path, sep='\t', header=True, index=False)