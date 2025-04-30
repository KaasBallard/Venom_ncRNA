# Last Edited: 2025/04/29

# Load packages ----
library(tidyverse)
library(arrow)

# Read data ----
# miRNA data
mirna_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Conversion table for expressed genes
conversion_table <- 'Data/Conversion_Table/Cvv_GTF_to_converted_names_2024.2.18.txt'

# Read the miRNA data
mirna_df <- read_parquet(file = mirna_data) |> 
  filter(
    str_detect(genes, 'Venom_'),
    # Filter by targeting strength
    total.score >= 155,
    total.energy <= -7
  ) |> 
  select(sample.id, miRNA.cluster, genes, venom.family, mRNA.counts, feature.type) |> 
  distinct() |> 
  mutate(
    # Create a column to say if a miRNA targets a given gene
    targeted = 'yes'
  ) |> 
  # Remove sample ids and mRNA counts as they are no longer necessary
  # Also remove the miRNA cluster column as it is no longer needed
  select(
    -sample.id, -mRNA.counts, -miRNA.cluster
  ) |> 
  distinct()

# Read in the conversion table to find all annotated venom genes
annotated_genes_df <- read_tsv(file = conversion_table) |> 
  # Convert Venom_ADAM28_1 to Venom_SVMP12
  mutate(
    converted_id_no_dups = ifelse(converted_id_no_dups == "Venom_ADAM28_1", 'Venom_SVMP12', converted_id_no_dups),
    converted_id_no_dups = ifelse(converted_id_no_dups == 'Venom_ADAM28_2', 'Venom_ADAM28', converted_id_no_dups)
  ) |> 
  # Get venoms
  filter(str_starts(converted_id_no_dups, 'Venom_')) |> 
  select(genes = converted_id_no_dups) |> 
  # Create an annotated column
  mutate(
    annotated = 'yes'
  )

# Fuse the miRNA-mRNA data to the annotated venoms ----
# Fuse the data frames
known_genes_df <- full_join(
  mirna_df,
  annotated_genes_df,
  by = 'genes',

) |> 
  filter(!str_detect(genes, 'ADAM')) |> 
  # Fix entries for PLAK, as it was the only one that wasn't targeted
  mutate(
    venom.family = ifelse(genes == 'Venom_PLA2K', 'PLA2', venom.family),
    targeted = ifelse(genes == 'Venom_PLA2K', 'no', targeted)
  )
# Only PLA2K doesn't have any strongly predicted binders

# Targeting by location ----
# Find for genes targeted by type
target_by_type_df <- known_genes_df |>
  pivot_wider(
    names_from = feature.type,
    names_prefix = 'targeted.in.',
    values_from = targeted,
    values_fill = 'no'
  ) |> 
  select(-targeted.in.NA) |> 
  arrange(venom.family)

# Write the file
write_csv(target_by_type_df, file = 'Tables/miRNAs_and_Targets/Filtered_venom_gene_targeting_by_region_2025.04.29.csv')

# Calculate the target counts
target_counts_df <- known_genes_df |> 
  filter(!is.na(feature.type)) |> 
  group_by(feature.type) |> 
  summarise(
    # Create a column containing the total venom genes
    total.genes = nrow(known_genes_df |> distinct(genes)),
    # Count the number of genes target by feature type
    total.by.target = n()
  ) |> 
  mutate(
    # Create a percent column
    percent.targeted = (total.by.target / total.genes) * 100,
    # Create a fraction column
    fraction.targeted = paste(as.character(total.by.target), as.character(total.genes), sep = '/')
  ) |> 
  select(feature.type, percent.targeted, fraction.targeted)
write_csv(target_counts_df, file = 'Tables/miRNAs_and_Targets/Filtered_venom_gene_targeted_by_miRNAs_fraction_2025.04.29.csv')
