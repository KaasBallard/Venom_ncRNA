#### Set up ####

# Last updated 2024/8/31

# Load in packages
library(tidyverse)
library(readxl)
library(Biostrings)
library(DESeq2)

# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd('C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/')

# Define the conversion table path
conversion_table <- 'Data/Conversion_Table/Cvv_GTF_to_converted_names_2024.2.18.txt'

# Define the path to the library we sent Anthony
protein_library <- 'Data/Protein/CroVir_rnd1.all.maker.proteins.final_withMyo_library_sent_to_Anthony.fasta'

#### Format the Conversion Table ####

# Read a second name conversion df in so that I can fuse and convert easily.
conversion_df2 <- read.table(file = conversion_table, header = T) %>% 
  rename_with(~ gsub('_', '.', .), contains('_'))

#### Convert FASTA names to a Data Table ####

# Read the FASTA file
fasta_data <- readAAStringSet(protein_library)

# Extract the headers from the FASTA file
fasta_names <- names(fasta_data)
head(print(fasta_names))

# Creata a data table of all the names from the FASTA file
fasta_df <- tibble(crovir.protein = fasta_names) %>% 
  mutate(crovir.transcript = gsub('-protein-', '-transcript-', crovir.protein)) # Convert transcripts to proteins

# Fuse the conversion_df2 with the fasta_df
library_check_df <- left_join(
  conversion_df2,
  fasta_df,
  by = c('crovir.transcript')
) %>% 
  dplyr::rename('Genes' = 'converted.id.no.dups')


# Filter for rows that have NA values in the crovir.protein column and print them
na_genes <- library_check_df %>% 
  filter(is.na(crovir.protein))

# Variable to hold the na_genes
not_sent_genes <- na_genes$Genes

# Print the result
print(not_sent_genes)


# Check what venom genes were not sent
venom_library_check_df <- library_check_df %>% 
  dplyr::filter(str_detect(Genes, 'Venom_'))

na_venom_genes <- venom_library_check_df %>% 
  filter(is.na(crovir.protein))

not_sent_venom_genes <- na_venom_genes$Genes

print(not_sent_venom_genes)
# Looks like the only venom genes that weren't sent were SVSP10 and SVSP11