# Last Edited: 2025/02/27

# Set up and Read Data ----

## Load packages ----

# Load in packages
library(compositions)
library(cowplot)
library(tidyverse)
library(scales)
# library(matrixStats)
library(ggrepel)
# library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(readxl)
library(ggpmisc)
library(ggplot2)
library(gganimate)
library(arrow)
# Source my functions
source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')

## Set data ----
# Expression data for venom genes and their miRNAs
mRNA_protein_data <- 'Data/Merged/mRNA_Protein_Combined_Data_2025.01.22.parquet'
file.exists(mRNA_protein_data)

## Read reference data ----

# Read reference data in as a DataFrame
mRNA_protein_df <- read_parquet(file = mRNA_protein_data) %>% 
  dplyr::filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_')
  )


# Set up DataFrame for BRMS modeling ----

## Reference data ----

# Apparently zero values are not treated properly by the clr transform, so I am going to create a pseudo-count to fix it.
# First I need to find the lowest value in each set of data though
# Find the lowest non-zero intensity value
min_nonzero_intensity <- mRNA_protein_df %>%
  filter(intensity > 0) %>%
  summarise(min_value = min(intensity)) %>%
  pull(min_value)
print(min_nonzero_intensity)

# Find the lowest non-zero RNA.VST value
min_nonzero_rna_counts <- mRNA_protein_df %>%
  filter(mRNA.counts > 0) %>%
  summarise(min_value = min(mRNA.counts)) %>%
  pull(min_value)
print(min_nonzero_rna_counts)

# Initialize pseudo-count for mRNA
mRNA_pseudo_count <- 1e-7

# Initialize pseudo-count for protein
protein_pseduo_count <- 1e-7

# Apparently zero values are not treated properly by the clr transform, so I have already filtered out proteins that were not been observed.
# Calculate CLR
clr_df <- mRNA_protein_df %>%
  # Get mRNA and Protein data seperate
  select(
    sample.id, genes, venom.family, mRNA.counts, intensity, protein.observed
  ) %>% 
  distinct() %>% 
  # Replace zeros with the small constant before CLR
  mutate(intensity = ifelse(intensity == 0, protein_pseduo_count, intensity)) %>%
  mutate(mRNA.counts = ifelse(mRNA.counts == 0, mRNA_pseudo_count, mRNA.counts)) %>% 
  group_by(sample.id) %>%
  mutate(Scaled.Protein.Expression = intensity / sum(intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.mRNA.Expression = mRNA.counts / sum(mRNA.counts)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  # Perform CLR transformation
  mutate(protein.clr = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(mRNA.clr = as.numeric(compositions::clr(Scaled.mRNA.Expression))) %>% 
  select(-Scaled.mRNA.Expression, -Scaled.Protein.Expression, -mRNA.counts, -intensity)
glimpse(clr_df)

# Save the above as a data frame to be read in another script
write_parquet(clr_df, sink = 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_venom_2025.01.30.parquet')


# Do a second transformation excluding proteins that were not observed
clr_df2 <- mRNA_protein_df %>% 
  # Filter unexpressed proteins
  dplyr::filter(
    protein.observed == 'Yes'
  ) %>%
  # Get mRNA and Protein data seperate
  select(
    sample.id, genes, venom.family, mRNA.counts, intensity, protein.observed
  ) %>% 
  distinct() %>%
  # Replace zeros with the small constant before CLR
  mutate(intensity = ifelse(intensity == 0, protein_pseduo_count, intensity)) %>%
  mutate(mRNA.counts = ifelse(mRNA.counts == 0, mRNA_pseudo_count, mRNA.counts)) %>% 
  group_by(sample.id) %>%
  mutate(Scaled.Protein.Expression = intensity / sum(intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.mRNA.Expression = mRNA.counts / sum(mRNA.counts)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  # Perform CLR transformation
  mutate(protein.clr = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(mRNA.clr = as.numeric(compositions::clr(Scaled.mRNA.Expression))) %>% 
  select(-Scaled.mRNA.Expression, -Scaled.Protein.Expression, -mRNA.counts, -intensity)
glimpse(clr_df2)

# Save
write_parquet(clr_df2, sink = 'Data/CLR_transformed_data/Filtered_CLR_mRNA_vs_protein_expression_venom_2025.02.27.parquet')
