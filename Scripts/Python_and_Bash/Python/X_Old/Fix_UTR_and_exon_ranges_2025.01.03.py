# %% [markdown]
# Step 1
# Date: 2024/12/19
# The purpose of this script is to take the GFF created by gffread with the following command:
# 
# gffread -E /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP.gtf -o /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Crotalus_viridis_annotation_with_BPP_and_myotoxin_2024.12.18.gff
# 
# and fix re-add the three_prime_utr and five_prime_utr sequences back.
# 
# Next step: BCFtools_get_fasta_file_three_prime_utr_2024.12.20.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BCFtools/BCFtools_get_fasta_file_three_prime_utr_2024.12.20.sh)
# Previous step: None

# %%
# Import needed packages
import pandas as pd
import polars as pl
import os

# %%
# Laod file path into memory
gtf = '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP.gtf'

# Load GTF into memory
gtf_df = pd.read_csv(
    gtf, sep = '\t', comment = '#', header = None, names= [
        'seqid', 'source', 'type2', 'start', 'end', 'score', 'strand', 'phase', 'attributes2'
    ]
)

# Laod the file path of the GFF into memory
gff = '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Crotalus_viridis_annotation_with_BPP_and_myotoxin_2024.12.18.gff'

# Load the GFF into memory
gff_df = pd.read_csv(
    gff, sep='\t', comment='#', header=None, names = [
        'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'
    ]
)

# %%
# Get common columns for the two data frames
shared_cols = gff_df.columns.intersection(gtf_df.columns).tolist()
print(shared_cols)

# %%
# Join the data frames and filter out any removing everthing unneaded from the GFF section
gxf_df = (pd.merge(
    gff_df, gtf_df,
    on = shared_cols,
    how = 'outer'
)
)
gxf_df

# %%
# Fix the ranges for the three_prime_utr and five_prime_utr and get proper feature values into the type column so that 3' and 5' UTRs aren't missing anymore
gxf_df2 = (
    gxf_df
        .query('type2 != "gene"') # Filter out instances of gene IDs in that row to reduce data complexity
        # Copy they three_prime_utr and five_prime_utrs into the main type column
        # .assign(
        #     type = lambda x: x['type2'].where(x['type2'].isin(['three_prime_utr', 'five_prime_utr']), x['type'])
        # )
        # .drop(columns = ['type2', 'attributes2']) # Remove columns I don't need
        .drop_duplicates() # Equivalent to distinct() in R
        # Add a new column to clarify what genes are what
        .assign(
            gene_id = lambda x: x['attributes2'].str.extract(r'gene_id "([^"]+)"')
        )
)
gxf_df2

# Create a dictionary to map gene_id to attributes where type2 is 'exon'
exon_attributes = gxf_df2[gxf_df2['type2'] == 'exon'].set_index('gene_id')['attributes'].to_dict()

# Replace missing values in 'attributes' for 'three_prime_utr' or 'five_prime_utr' rows when 'type2' is 'exon'
gxf_df2['attributes'] = gxf_df2.apply(
    lambda row: exon_attributes.get(row['gene_id'], row['attributes']) if pd.isna(row['attributes']) else row['attributes'],
    axis=1
)

# %%
gxf_df2

# %%
# Separate CDS from everthing else
cds_df = (
    gxf_df2
        .query("type == 'exon'")
        .drop_duplicates()
        .query("gene_id == 'myotoxin1'") # Filter out other gene_ids as a test if the function works
)

# Seperate 3' UTRs from everthing else
three_utr_df = gxf_df2.query("type2 == 'three_prime_utr'").drop_duplicates().query("gene_id == 'myotoxin1'") # Filter out other gene_ids as a test if the function works

# Seperate 5' UTRs from everthing else
five_utr_df = gxf_df2.query("type2 == 'five_prime_utr'").drop_duplicates().query("gene_id == 'myotoxin1'") # Filter out other gene_ids as a test if the function works

# Get 3' and 5' UTRs
three_five_utr_df = gxf_df2.query("type2 == 'three_prime_utr' or type2 == 'five_prime_utr'").drop_duplicates().query("gene_id == 'myotoxin1'")

# %%
# Define a function that takes a data frame and re-calculates ranges for the exon
def range_recalc(df1, df2):

    # Iterate through rows, checking if the data ranges overlap
    # df1 is the exon df
    for i, rows1 in df1.iterrows():

        # df2 is the 3'UTR df
        for j, rows2 in df2.iterrows():
            
            # Check if the row is a three_prime_utr row, and if so do this stuff:
            if rows2['type2'] == 'three_prime_utr':

                # Check if the sequence IDs are the same and only continue if true
                if (
                    rows1['seqid'] == rows2['seqid'] 
                    and rows1['type'] == 'exon' 
                    and rows1['gene_id'] == rows2['gene_id']
                    and rows1['source'] == rows2['source']
                ):

                    print(f"Matching seqid found: {rows1['seqid']}")

                    # Check if what strand the sequence is on
                    # If the strand is sense do that
                    if rows1['strand'] == '+':
                        # Check if the ranges overlap
                        if rows1['start'] < rows2['start'] and rows1['end'] == rows2['end']:
                            # Adjust the end point for row1 (the exon)
                            df1.at[i, 'end'] = rows2['start'] - 1
                            print(f"Adjusted end for row {i}: {df1.at[i, 'end']}")
                        else:
                            print(f"Ranges do not overlap for row {i}")
                            continue

                    # If the strand is anti-sense do this
                    elif rows1['strand'] == '-':
                        # Check if the ranges overlap
                        if rows1['start'] == rows2['start'] and rows1['end'] > rows2['end']:
                            # Adjust the start point for row1 (the exon)
                            df1.at[i, 'start'] = rows2['end'] + 1 
                            print(f"Adjusted start for row {i}: {df1.at[i, 'start']}")
                        else:
                            print(f"Ranges do not overlap for row {i}")
                            continue 
                    else:
                        continue
                else:
                    continue


            # Check if the row is a five_prime_utr row, and if so do this stuff:        
            elif rows2['type2'] == 'five_prime_utr':

                # Check if the sequence IDs are the same and only continue if true
                if (
                    rows1['seqid'] == rows2['seqid'] 
                    and rows1['type'] == 'exon' 
                    and rows1['gene_id'] == rows2['gene_id']
                    and rows1['source'] == rows2['source']
                ):

                    print(f"Matching seqid found: {rows1['seqid']}")

                    # Check if what strand the sequence is on
                    # If the strand is sense do that
                    if rows1['strand'] == '+':
                        # Check if the ranges overlap
                        if rows1['start'] == rows2['start'] and rows1['end'] > rows2['end']:
                            # Adjust the start point for row1 (the exon)
                            df1.at[i, 'start'] = rows2['end'] + 1 
                            print(f"Adjusted start for row {i}: {df1.at[i, 'start']}")
                        else:
                            print(f"Ranges do not overlap for row {i}")
                            continue 
                    
                    # If the strand is anti-sense do this
                    elif rows1['strand'] == '-':
                        # Check if the ranges overlap
                        if rows1['start'] < rows2['start'] and rows1['end'] == rows2['end']:
                            # Adjust the end point for rows1 (the exon)
                            df1.at[i, 'end'] == rows2['start'] - 1
                            print(f"Adjusted end for row {i}: {df1.at[i, 'end']}")
                        else:
                            print(f"Ranges do not overlap for row {i}")
                            continue
                    else:
                        continue
                else:
                    continue

            else:
                continue

    return df1

# %%
# Filter out everything but myotoxin from gxf_df2
gxf_myo_df = gxf_df2.query("gene_id == 'myotoxin1'")
gxf_myo_df

# Seperate 3' UTRs from everthing else
three_utr_df = gxf_myo_df.query("type2 == 'three_prime_utr'").drop_duplicates() 
# Seperate 5' UTRs from everthing else
five_utr_df = gxf_myo_df.query("type2 == 'five_prime_utr'").drop_duplicates()

# Get 3' and 5' UTRs
three_five_utr_df = gxf_myo_df.query("type2 == 'three_prime_utr' or type2 == 'five_prime_utr'").drop_duplicates()

# %%
# Use the function to re-calculate the start or end for the CDS
three_five_utr_recalc_df = range_recalc(gxf_myo_df, three_five_utr_df)
# print(type(three_five_utr_recalc_df))

# %%
three_five_utr_df

# %%
three_five_utr_recalc_df

# %%
# Now that I know it works, I can apply the function to the entire data frame and give it a second data frame that contains both 3' and 5' UTRs
# Seperate three_prime_utr and five_prime_utrs
utrs_df = (
    gxf_df2
        .query("type2 in ['three_prime_utr', 'five_prime_utr']")
        .drop_duplicates()
)


# %%
# Apply the function to the main data frame and the above data frame
recscaled_cds_df = range_recalc(gxf_df2, utrs_df)
recscaled_cds_df

# %%
# Format the recalcualted CDS data frame
formated_recscaled_cds_df = (
    recscaled_cds_df
        .get(['seqid', 'source', 'type2', 'start', 'end', 'score', 'strand', 'phase', 'attributes']) # Reorder and remove columns
        .rename(columns = {'type2': 'type'}) # Rename columns
)
formated_recscaled_cds_df


# %%
# Save the new gff file
formated_recscaled_cds_df.to_csv(
    '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Crotalus_viridis_annotation_with_BPP_and_myotoxin_with_three_and_five_prime_utrs_2024.12.18.gff',
    sep = '\t',
    index = False,
    header = False
)


