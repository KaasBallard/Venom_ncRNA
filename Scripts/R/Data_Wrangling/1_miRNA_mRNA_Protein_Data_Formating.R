# Last updated: 2025/01/20

# Set up ----

## Set up environment ----

# Load in packages
library(tidyverse)
library(readxl)
library(Biostrings)
library(DESeq2)
library(arrow)


## Set up data ----

# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd('C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/')

# Define file paths for different data sets
metadata <- 'Data/metadata/DESeq_Sample_metadata_2024.6.4.csv'
rnaseq_data <- 'Data/mRNA/cvv_12SnakeVenomRNAseq_rawCounts_09.16.22.txt'
protein_data <- 'Data/Protein/2023_Aug9_QEHF_Castoe.xlsx'
conversion_table <- 'Data/Conversion_Table/Cvv_GTF_to_converted_names_2024.2.18.txt' # nolint: line_length_linter.
miRNA_reference_data <- 'Data/miRNA/Proccessed/miRNA_counts_for_reference/reference_target_and_count_data.2025.01.22.parquet'
miRNA_number_per_sample_data <- 'Data/miRNA/Proccessed/Number_of_miRNAs_per_sample/all_miRNA_numbers_per_sample-gene-feature_type.2025.01.22.parquet'
no_feature_type_sample_data <- 'Data/miRNA/Proccessed/Number_of_miRNAs_per_sample/all_miRNA_numbers_per_sample-gene.2025.01.22.parquet'
filt_miRNA_number_per_sample <- 'Data/miRNA/Proccessed/Number_of_miRNAs_per_sample/filtered_miRNA_numbers_per_sample-gene-feature_type.2025.01.22.parquet'
filt_no_feature_type_sample_data <- 'Data/miRNA/Proccessed/Number_of_miRNAs_per_sample/filtered_miRNA_numbers_per_sample-gene.2025.01.22.parquet'
protein_library <- 'Data/Protein/CroVir_rnd1.all.maker.proteins.final_withMyo_library_sent_to_Anthony.fasta' # Define the path to the library we sent Anthony


# Conversion table ----

# Read a second name conversion df in so that I can fuse and convert easily.
conversion_df <- read.table(file = conversion_table, header = TRUE) %>%
  rename_with(~ gsub('_', '.', .), contains('_')) %>%
  # Change the name of the venom adam genes so that they are correct
  mutate(
    converted.id.no.dups = ifelse(converted.id.no.dups == "Venom_ADAM28_1", 'Venom_SVMP12', converted.id.no.dups),
    converted.id.no.dups = ifelse(converted.id.no.dups == 'Venom_ADAM28_2', 'Venom_ADAM28', converted.id.no.dups)
  )
glimpse(conversion_df)

# Read the miRNA sample data (number of miRNAs per sample version) ----

## Read the miRNA sample data ----

# Read the miRNA sample data for the number of miRNAs per sample
miRNA_sample_df <- read_parquet(file = miRNA_number_per_sample_data) %>%
  rename_with(~ gsub('_', '.', .), contains('_')) %>%
  # Change the name of the venom adam genes so that they are correct
  mutate(
    genes = ifelse(genes == "Venom_ADAM28_1", 'Venom_SVMP12', genes),
    genes = ifelse(genes == 'Venom_ADAM28_2', 'Venom_ADAM28', genes)
  )
glimpse(miRNA_sample_df)

### Read the miRNA sample data (no feature type version) ----
no_feature_miRNA_sample_df  <- read_parquet(file = no_feature_type_sample_data) %>%
  rename_with(~ gsub('_', '.', .), contains('_')) %>%
  # Change the name of the venom adam genes so that they are correct
  mutate(
    genes = ifelse(genes == "Venom_ADAM28_1", 'Venom_SVMP12', genes),
    genes = ifelse(genes == 'Venom_ADAM28_2', 'Venom_ADAM28', genes)
  )
glimpse(no_feature_miRNA_sample_df)

## Read the filtered miRNA sample data ----

# Read the filtered miRNA sample data for the number of miRNAs per sample
filtered_miRNA_sample_df <- read_parquet(file = filt_miRNA_number_per_sample) %>%
  rename_with(~ gsub('_', '.', .), contains('_')) %>%
  # Change the name of the venom adam genes so that they are correct
  mutate(
    genes = ifelse(genes == "Venom_ADAM28_1", 'Venom_SVMP12', genes),
    genes = ifelse(genes == 'Venom_ADAM28_2', 'Venom_ADAM28', genes)
  )
glimpse(filtered_miRNA_sample_df)

### Read the filtered miRNA sample data (no feature type version) ----

# Read the filtered miRNA sample data for the number of miRNAs per sample
no_ft_filtered_miRNA_sample_df <- read_parquet(file = filt_no_feature_type_sample_data) %>%
  rename_with(~ gsub('_', '.', .), contains('_')) %>%
  # Change the name of the venom adam genes so that they are correct
  mutate(
    genes = ifelse(genes == "Venom_ADAM28_1", 'Venom_SVMP12', genes),
    genes = ifelse(genes == 'Venom_ADAM28_2', 'Venom_ADAM28', genes)
  )
glimpse(no_ft_filtered_miRNA_sample_df)



# Read and format miRNA data reference data ----

## Read miRNA reference data ----

# Read the miRNA data frame I created from the shortstack counts.txt, miRanda tab output, and bedtools intersect between the miRanda tab output and the genome gtf.
miRNA_df <- read_parquet(file = miRNA_reference_data) %>%
  select(
    -assembler
  ) %>%
  rename_with(~ gsub('_', '.', .), contains('_')) %>%
  # Change the name of the venom adam genes so that they are correct
  mutate(
    genes = ifelse(genes == "Venom_ADAM28_1", 'Venom_SVMP12', genes),
    genes = ifelse(genes == 'Venom_ADAM28_2', 'Venom_ADAM28', genes)
  )
glimpse(miRNA_df)



### Calculate reads per million for miRNAs ----

# Load new data frame with count data only
miRNA_df <- miRNA_df %>%
  group_by(sample.id) %>%
  mutate(
    # Create a column containing the total counts
    total.counts = sum(miRNA.counts),
    # Normalize to RPM
    miRNA.rpm = (miRNA.counts / total.counts) * 1e6
  ) %>%
  ungroup() %>%
  select(-total.counts)
glimpse(miRNA_df)



# Read and format mRNA Data ----

# Read raw mRNA data counts in
mRNA_df <- read.table(rnaseq_data, header = TRUE)
glimpse(mRNA_df)

# Change the Crovir_Transcript_ID entry to what it is in the conversion_table
mRNA_df[1, "Crovir_Transcript_ID"] <- 'crovir-transcript-myotoxin' # For myotoxin
mRNA_df[18188, 'Crovir_Transcript_ID'] <- 'crovir-transcript-17319'

# Filter out uneccessary columns
mRNA_df <- mRNA_df %>%
  select(
    'Geneid', 'Crovir_Transcript_ID',
    '..STAR_mapped.RVG_5Aligned.sortedByCoord.out.bam', '..STAR_mapped.RVG_6Aligned.sortedByCoord.out.bam', '..STAR_mapped.RVG_7Aligned.sortedByCoord.out.bam', '..STAR_mapped.RVG_12Aligned.sortedByCoord.out.bam',
    '..STAR_mapped.LVG_2Aligned.sortedByCoord.out.bam', '..STAR_mapped.LVG_4Aligned.sortedByCoord.out.bam', '..STAR_mapped.LVG_9Aligned.sortedByCoord.out.bam'
  ) %>%
  dplyr::rename(
    'gtf.gene' = 'Geneid', 'crovir.transcript' = 'Crovir_Transcript_ID',
    'CV1087_viridis' = '..STAR_mapped.RVG_5Aligned.sortedByCoord.out.bam', 'CV0987_lutosus' = '..STAR_mapped.RVG_6Aligned.sortedByCoord.out.bam',
    'CV0985_concolor' = '..STAR_mapped.RVG_7Aligned.sortedByCoord.out.bam', 'CV1082_viridis' = '..STAR_mapped.RVG_12Aligned.sortedByCoord.out.bam',
    'CV1081_viridis' = '..STAR_mapped.LVG_2Aligned.sortedByCoord.out.bam', 'CV0857_viridis' = '..STAR_mapped.LVG_4Aligned.sortedByCoord.out.bam',
    'CV1086_viridis' = '..STAR_mapped.LVG_9Aligned.sortedByCoord.out.bam'
  )
glimpse(mRNA_df)

# Create a variable to hold shared column names between the conversion table and the above
mRNA_conversion_cols <- intersect(names(mRNA_df), names(conversion_df))

# Join the mRNA data with the conversion table
mRNA_df <- left_join(
  mRNA_df,
  conversion_df,
  by = mRNA_conversion_cols
) %>%
  select(
    converted.id.no.dups, converted.id, gtf.gene, gtf.gene.trimmed, contains('CV')
  ) %>%
  dplyr::rename('genes' = 'converted.id.no.dups')
glimpse(mRNA_df)



# DESeq2 ----

## Fuse miRNA and mRNA data to do VST conversion through DESeq2 ----

# Format mRNA data to be bound to the miRNA data
mRNA_deseq_df <- mRNA_df %>%
  select(
    genes, CV1081_viridis, CV0857_viridis, CV1086_viridis, CV1082_viridis, CV1087_viridis, CV0987_lutosus, CV0985_concolor
  ) %>%
  filter(!is.na(genes)) %>% # Remove empty values
  column_to_rownames(var = 'genes')
glimpse(mRNA_deseq_df)

# Format miRNA data to be bound to mRNA data
miRNA_deseq_df <- miRNA_df %>%
  dplyr::select(
    miRNA.cluster.original,
    sample.id,
    miRNA.counts
  ) %>%
  distinct() %>%
  # Pivot as the current format of the columns is by sample.id counts, not sample.ids as columns with counts in them
  pivot_wider(
    names_from = sample.id,
    values_from = miRNA.counts,
    values_fill = 0
  ) %>%
  column_to_rownames(var = 'miRNA.cluster.original')
glimpse(miRNA_deseq_df)


## Run DESeq2 ----

# Concatenated miRNA and RNA data
counts_matrix <- bind_rows(mRNA_deseq_df, miRNA_deseq_df)
rm(mRNA_deseq_df, miRNA_deseq_df)
gc()

# Read sample data in
sample_data <- read.csv(metadata, row.names = 1) %>% select(species, sub_species, sex)
sample_data$species <- factor(sample_data$species) # Turn the sample data into factors
sample_data$sub_species <- factor(sample_data$sub_species) # Turn the sample data into factors
sample_data$sex <- factor(sample_data$sex)

# Create the deseq data set
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_data,
  design = ~ species
)


## Get normalized counts for miRNAs and mRNAs ----

# Do the vst, rlog, and normTransform transformations
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
ntd <- normTransform(dds)

# Extract the vst transformed values
vst_both <- as.data.frame(assay(vsd))
rlog_both <- as.data.frame(assay(rld))
ntd_both <- as.data.frame(assay(ntd))

# Turn rows into a column and pivot
vst_both <- rownames_to_column(vst_both, var = "genes") %>%
  # Pivot the sample columns so that Sample.IDs are in a single column and expression is given it's own column
  pivot_longer(
    cols = matches('CV'),
    names_to = 'sample.id',
    values_to = 'vst'
  ) %>%
  distinct()
rlog_both <- rownames_to_column(rlog_both, var = 'genes') %>%
  # Pivot the sample columns so that Sample.IDs are in a single column and expression is given it's own column
  pivot_longer(
    cols = matches('CV'),
    names_to = 'sample.id',
    values_to = 'rlog'
  ) %>%
  distinct()
ntd_both <- rownames_to_column(ntd_both, var = 'genes') %>%
  # Pivot the sample columns so that Sample.IDs are in a single column and expression is given it's own column
  pivot_longer(
    cols = matches('CV'),
    names_to = 'sample.id',
    values_to = 'ntd'
  ) %>%
  distinct()



### Fuse the different count types ----

# Create a DataFrame containing all the count types
norm_counts_both <- full_join(
  vst_both,
  rlog_both,
  by = c('genes', 'sample.id')
) %>%
  full_join(
    ntd_both,
    by = c('genes', 'sample.id')
  )
glimpse(norm_counts_both)


### Get normalized counts for miRNA and mRNA ----

# Create a data frame that contains vst values for only mRNA
norm_counts_mrna <- norm_counts_both %>%
  dplyr::filter(
    !str_detect(genes, 'Cluster_')
  ) %>%
  dplyr::rename(
    mRNA.vst = vst, mRNA.rlog = rlog, mRNA.ntd = ntd
  ) %>%
  distinct()
glimpse(norm_counts_mrna)

# Create a data frame that contains vst values for only mirna
norm_counts_mirna <- norm_counts_both %>%
  dplyr::filter(
    str_detect(genes, 'Cluster_')
  ) %>%
  dplyr::rename(miRNA.cluster.original = genes, miRNA.vst = vst, miRNA.rlog = rlog, miRNA.ntd = ntd)
glimpse(norm_counts_mirna)
rm(norm_counts_both)


### Fuse the counts back to their source data ----

# Fuse mRNA counts back to the main data
mRNA_df <- mRNA_df %>%
  pivot_longer(
    cols = contains('CV'),
    names_to = 'sample.id',
    values_to = 'mRNA.counts'
  ) %>%
  full_join(
    norm_counts_mrna,
    by = c('genes', 'sample.id')
  ) %>%
  mutate(
    venom.family = case_when(
      grepl('Venom_SVMP', genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('Venom_VEGF', genes) ~ 'VEGF',
      grepl('Venom_ohanin', genes) ~ 'ohanin',
      grepl('Venom_vQC', genes) ~ 'vQC',
      grepl('Venom_SVSP', genes) ~ 'SVSP',
      grepl('Venom_PLA2', genes) ~ 'PLA2',
      grepl('Venom_ADAM', genes) ~ 'ADAM',
      grepl('Venom_CRISP', genes) ~ 'CRISP',
      grepl('Venom_CTL', genes) ~ 'CTL',
      grepl('Venom_EXO', genes) ~ 'EXO',
      grepl('Venom_LAAO', genes) ~ 'LAAO',
      grepl('Venom_myotoxin', genes) ~ 'myotoxin',
      grepl('Venom_BPP', genes) ~ 'BPP',
      TRUE ~ 'others'
    )
  )
glimpse(mRNA_df)

# Fuse miRNA normalized values to the main data
miRNA_df <- full_join(
  miRNA_df,
  norm_counts_mirna,
  by = c('sample.id', 'miRNA.cluster.original'),
  relationship = 'many-to-many'
) %>%
  mutate(
    venom.family = case_when(
      grepl('Venom_SVMP', genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('Venom_VEGF', genes) ~ 'VEGF',
      grepl('Venom_ohanin', genes) ~ 'ohanin',
      grepl('Venom_vQC', genes) ~ 'vQC',
      grepl('Venom_SVSP', genes) ~ 'SVSP',
      grepl('Venom_PLA2', genes) ~ 'PLA2',
      grepl('Venom_ADAM', genes) ~ 'ADAM',
      grepl('Venom_CRISP', genes) ~ 'CRISP',
      grepl('Venom_CTL', genes) ~ 'CTL',
      grepl('Venom_EXO', genes) ~ 'EXO',
      grepl('Venom_LAAO', genes) ~ 'LAAO',
      grepl('Venom_myotoxin', genes) ~ 'myotoxin',
      grepl('Venom_BPP', genes) ~ 'BPP',
      TRUE ~ 'others'
    )
  )
glimpse(miRNA_df)


## Save the miRNA expression and target data ----

# Save the miRNA_df
write_parquet(miRNA_df, 'Data/miRNA/Proccessed/miRNA_counts_for_reference/miRNA_counts-rpm-vst_for_reference.2025.01.22.parquet')

# Remove some columns to decrease the size of the data
miRNA_reduced_df <- miRNA_df %>%
  select(
    miRNA.cluster.original, miRNA.cluster, best.miRNA.ortholog, base.miRNA.name, miRNA.sequence.chrom, miRNA.start, miRNA.end, miRNA.length, miRNA.strandedness,
    sample.id, miRNA.counts, miRNA.rpm, miRNA.vst, miRNA.rlog, miRNA.ntd
  ) %>%
  distinct()
glimpse(miRNA_reduced_df)
write_parquet(miRNA_reduced_df, 'Data/miRNA/Proccessed/miRNA_counts_for_reference/miRNA_counts-rpm-vst-rlog-norm_counts.2025.01.22.parquet')



## Save the mRNA expression data ----

# Format the data, removing uneccessary columns
mRNA_reduced_df <- mRNA_df %>%
  select(
    -converted.id, -gtf.gene, -gtf.gene.trimmed
  ) %>%
  distinct()
glimpse(mRNA_reduced_df)

# Write the data to file
write_parquet(mRNA_reduced_df, 'Data/mRNA/mRNA_raw_counts-vst-rlog-norm_counts.2025.01.22.parquet')



# Read and format protein data ----

# Read protein data as a data frame
protein_df <- read_excel(protein_data, sheet = 1)

# Create a data frame to map Anthony's protein names to standard names
# Convert anthony names to standard names (version 1)
protein_name_conversion <- data.frame(
  stringsAsFactors = FALSE,
  Anthony_name = c("C_o_cerberus_AZ_Kingman",
                   "C_o_cerberus_AZ_NW_Kingman",
                   "C_o_concolor_UT_Duchesne",
                   "C_o_lutosus_528_UT_Beaver_Co",
                   "C_o_lutosus_UT_Beaver_Co",
                   "C_o_tigris_AZ_Pima_Co",
                   "C_v_viridis_CO_Hwy93",
                   "C_v_viridis_CO_Otero_Co",
                   "C_v_viridis_CO_Weld_Co",
                   "C_v_viridis_MT",
                   "C_v_viridis_NM_6",
                   "C_v_viridis_NM_Luna_Co",
                   "C_v_viridis_NM_Socorro",
                   "C_v_viridis_CO_San_Miguel_Co"),
  New_name = c("NONE",
               "CV1090_cerberus_Other_M",
               "RVG_7S_CV0985_concolor_Other_F",
               "RVG_6S_CV0987_lutosus_Other_M",
               "NONE",
               "NONE",
               "CV1096_viridis_North_F",
               "CV1084_viridis_Mid_M",
               "RVG_5S_CV1087_viridis_North_F",
               "CV0857_viridis_North_M",
               "NONE",
               "CV1086_viridis_South_M",
               "CV1089_viridis_South_M",
               "CV1081_viridis_Mid_M")
)

# Process protein data: filter, rename, and remove certain columns
pivoted_protein_df <- protein_df %>%
  select(-'Indistinguishable Proteins') %>% # remove indistinguishable proteins column # 51
  pivot_longer(
    cols = 9:50,
    names_sep = '\\s',
    names_to = c("Sample.ID", "Feature"),
    values_to = "Intensity"
  ) %>% #9:50
  left_join(protein_name_conversion, by = c('Sample.ID' = 'Anthony_name')) %>%
  filter(!New_name == 'NONE') %>% # remove samples for which we have no study samples
  select(-Sample.ID) %>%
  dplyr::rename(Sample.ID = New_name) %>%
  filter(!is.na(Protein)) # It seems that CV1082 doesn't have any data for this part.
glimpse(pivoted_protein_df)
# Forget protein_df
rm(protein_df)

# Samples that I have
protein_samples <- c("RVG_7S_CV0985_concolor_Other_F", 'RVG_6S_CV0987_lutosus_Other_M', 'RVG_5S_CV1087_viridis_North_F',
                     'CV1081_viridis_Mid_M', 'CV0857_viridis_North_M', 'CV1086_viridis_South_M')

# Filter out rows that don't have individuals I need and change the name of the Protein column and it's entries so it can be fused
sample_filtered_protein_df <- pivoted_protein_df %>%
  filter(Sample.ID %in% protein_samples) %>%
  dplyr::rename('crovir.transcript' = 'Protein') %>%
  mutate('crovir.transcript' = gsub('-protein-', '-transcript-', crovir.transcript))
glimpse(sample_filtered_protein_df)
rm(pivoted_protein_df)

# Change the conversion_df to have a better name for that column
conversion_df_for_protein <- conversion_df %>%
  dplyr::rename('Genes' = 'converted.id.no.dups') %>%
  filter(!grepl('fgenesh', gtf.gene))

# List of possible Sample.ID entries
study_sample_ids <- c(
  'RVG.7S.CV0985.concolor.Other.F',
  'RVG.6S.CV0987.lutosus.Other.M',
  'RVG.5S.CV1087.viridis.North.F',
  'LVG.2.CV1081.viridis.Mid.M',
  'LVG.4.CV0857.viridis.North.M',
  'LVG.9.CV1086.viridis.South.M'
)

# Create a function to duplicate rows with NA in Sample.ID for each possible Sample.ID
expand_na_rows <- function(df, sample_ids) {
  df %>%
    filter(is.na(Sample.ID)) %>%
    tidyr::uncount(length(sample_ids)) %>%
    mutate(Sample.ID = rep(sample_ids, length.out = n()))
}

# This is better because it fuses the data together without getting rid of anything
protein_df_with_ids <- full_join(
  sample_filtered_protein_df,
  conversion_df_for_protein,
  by = 'crovir.transcript'
) %>%
  rename_with(~ gsub(' ', '.', .), contains(' ')) %>%
  mutate(Sample.ID = case_when(
    Sample.ID == 'RVG_7S_CV0985_concolor_Other_F' ~ 'RVG.7S.CV0985.concolor.Other.F',
    Sample.ID == 'RVG_6S_CV0987_lutosus_Other_M' ~ 'RVG.6S.CV0987.lutosus.Other.M',
    Sample.ID ==  'RVG_5S_CV1087_viridis_North_F' ~ 'RVG.5S.CV1087.viridis.North.F',
    Sample.ID == 'CV1081_viridis_Mid_M' ~ 'LVG.2.CV1081.viridis.Mid.M',
    Sample.ID == 'CV0857_viridis_North_M' ~ 'LVG.4.CV0857.viridis.North.M',
    Sample.ID == 'CV1086_viridis_South_M' ~ 'LVG.9.CV1086.viridis.South.M',
    TRUE ~ Sample.ID
  ))
glimpse(protein_df_with_ids)

# Expand rows with NA in Sample.ID
expanded_rows <- expand_na_rows(protein_df_with_ids, study_sample_ids)

# Bind the expanded rows back to the original dataframe
protein_df_with_ids <- protein_df_with_ids %>%
  filter(!is.na(Sample.ID)) %>%
  bind_rows(expanded_rows) %>%
  mutate(
    Protein.Observed = ifelse(is.na(Intensity), 'No', 'Yes'),
    Feature = ifelse(is.na(Feature), 'Intensity', Feature),
    Intensity = ifelse(is.na(Intensity), 0, Intensity)
  ) %>% # This mutate function call adds the "Intensity" entry to the Feature column and makes the Intensity column have zeros in it if there is an NA those columns. This allows me to treat un-observed proteins like they were observed, just with zero intensity
  distinct()

# Add the some of the code from the "6_Protein_Library_Check.R" script

# Read a second name conversion df in so that I can fuse and convert easily.
conversion_df2 <- read.table(file = conversion_table, header = TRUE) %>%
  rename_with(~ gsub('_', '.', .), contains('_'))

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
rm(conversion_df2)
gc()

# Filter for rows that have NA values in the crovir.protein column and print them
na_genes <- library_check_df %>%
  filter(is.na(crovir.protein))

# Variable to hold the na_genes
not_sent_genes <- na_genes$Genes

# Print the result
print(not_sent_genes)

# Add another column to the protein_df_with_ids so that I can filter out anything that wasn't in the sent library
protein_df_with_ids <- protein_df_with_ids %>%
  mutate(
    In.Library = ifelse(Genes %in% not_sent_genes, 'No', 'Yes')
  ) %>%
  filter(Feature == 'Intensity') %>%
  select(-contains('Peptide'), -contains('Probability'), -contains('Count'), -Feature) %>%
  dplyr::rename(
    intensity = Intensity,
    sample.id = Sample.ID,
    in.library = In.Library,
    protein.observed = Protein.Observed,
    genes = Genes
  ) %>%
  mutate(
    sample.id = case_when(
      sample.id == 'RVG.7S.CV0985.concolor.Other.F' ~ 'CV0985_concolor',
      sample.id == 'RVG.6S.CV0987.lutosus.Other.M' ~ 'CV0987_lutosus',
      sample.id == 'RVG.5S.CV1087.viridis.North.F' ~ 'CV1087_viridis',
      sample.id == 'LVG.2.CV1081.viridis.Mid.M' ~ 'CV1081_viridis',
      sample.id == 'LVG.4.CV0857.viridis.North.M' ~ 'CV0857_viridis',
      sample.id == 'LVG.9.CV1086.viridis.South.M' ~ 'CV1086_viridis',
      TRUE ~ sample.id
    )
  )
glimpse(protein_df_with_ids)

# Forget conversion_df_for_protein
rm(conversion_df_for_protein, sample_filtered_protein_df)

## Save the protein data ----

# Remove unnecessary columns
protein_df_with_ids2 <- protein_df_with_ids %>%
  select(
    genes, sample.id, intensity, protein.observed, 'in.library'
  ) %>%
  distinct()
glimpse(protein_df_with_ids2)

# Save the data
write_parquet(protein_df_with_ids2, 'Data/Protein/protein_intensity_data.2025.01.22.parquet')


# Fuse mRNA, miRNA, and protein data ----

## Fuse mRNA and Protein data ----

# Create a vector of shared columns
protein_mRNA_cols <- intersect(names(mRNA_df), names(protein_df_with_ids))
protein_mRNA_cols

# Fuse data
protein_mRNA_df <- full_join(
  mRNA_df,
  protein_df_with_ids,
  by = protein_mRNA_cols
) %>%
  select(-Protein.Length) %>%
  distinct()
glimpse(protein_mRNA_df)

# # Check what can't be joined
# protein_mRNA_df2 <- anti_join(
#   mRNA_df,
#   protein_df_with_ids,
#   by = protein_mRNA_cols
# ) %>%
#   filter(str_detect(genes, 'Venom'))
# glimpse(protein_mRNA_df2)
# # It's just the sample that doesn't exist in the protein data
rm(mRNA_df, protein_df_with_ids)
gc()

### Format and save protein and mRNA data ----
protein_mRNA_df2 <- protein_mRNA_df %>%
  select(
    -contains('gtf.gene'), -converted.id, -crovir.transcript
  )
write_parquet(protein_mRNA_df2, 'Data/Merged/mRNA_Protein_Combined_Data_2025.01.22.parquet')
rm(protein_mRNA_df2)
gc()


## Fuse miRNA data to mRNA-Protein data ----

# Remove the sample without protein
protein_mRNA_df3 <- protein_mRNA_df %>%
  filter(!(sample.id == 'CV1082_viridis'))
glimpse(protein_mRNA_df3)

# Get shared columns
miRNA_mRNA_prot_cols <- intersect(names(protein_mRNA_df3), names(miRNA_df))
miRNA_mRNA_prot_cols

# Fuse miRNA data to mRNA-Protein data
protein_mRNA_miRNA_df <- full_join(
  protein_mRNA_df3,
  miRNA_df %>% 
    filter(sample.id != 'CV1082_viridis'), # Remove the sample that doesn't have any protein data
  by = miRNA_mRNA_prot_cols
) %>%
  # Select columns to change order and only get useful info
  select(
    sample.id, genes, venom.family, miRNA.target.chrom, miRNA.target.start, miRNA.target.end, positions, miRNA.target.length, strand, mRNA.counts, mRNA.vst, mRNA.rlog, mRNA.ntd, intensity, in.library, protein.observed,
    total.score, total.energy, max.score, max.energy, feature.type,
    miRNA.cluster.original, miRNA.cluster, best.miRNA.ortholog, miRNA.sequence.chrom, miRNA.start, miRNA.end, miRNA.strandedness, miRNA.length, miRNA.counts, miRNA.rpm, miRNA.vst, miRNA.rlog, miRNA.ntd,
    miRNA.name.probability, blast.percent.identity, E.value, bit.score, miRNA.sequence, hairpin.sequence
  ) %>%
  distinct()
glimpse(protein_mRNA_miRNA_df)
# Save the file
write_parquet(protein_mRNA_miRNA_df, 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet')

# Create a version of the data that only has venom genes and their non-venom paralogous
venom_protein_mRNA_miRNA_df <- protein_mRNA_miRNA_df %>%
  filter(str_starts(genes, 'Venom_|PLA2|ADAM|PRSS|BDEF')) %>%
  distinct()
# Save the files
write_parquet(venom_protein_mRNA_miRNA_df, 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_Venom_2025.01.22.parquet')

# Find what is missing from the data
missing_df <- protein_mRNA_miRNA_df %>%
  filter(if_any(everything(), is.na))
glimpse(missing_df)
View(missing_df)

# Add the mRNA-Protein data to the miRNA sample data ----

## Fuse miRNA sample data to mRNA-Protein data ----
# Get shared columns
miRNA_mRNA_prot_intersect <- intersect(names(miRNA_sample_df), names(protein_mRNA_df3))
miRNA_mRNA_prot_intersect

# Join the data
miRNA_sample_protein_mRNA_df <- full_join(
  miRNA_sample_df,
  protein_mRNA_df3,
  by = miRNA_mRNA_prot_intersect,
  # relationship = 'many-to-many'
) %>%
  select(
    -contains('gtf'), -converted.id, -crovir.transcript
  ) %>%
  distinct()
glimpse(miRNA_sample_protein_mRNA_df)

# Save the data
write_parquet(miRNA_sample_protein_mRNA_df, 'Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet')

## Fuse miRNA sample data to mRNA-Protein data (no feature type version) ----
# Get shared columns
no_feature_type_intersect <- intersect(names(no_feature_miRNA_sample_df), names(protein_mRNA_df3))
no_feature_type_intersect

# Join the data
no_ft_miRNA_protein_mRNA_df <- full_join(
  no_feature_miRNA_sample_df,
  protein_mRNA_df3,
  by = no_feature_type_intersect,
  # relationship = 'many-to-many'
) %>%
  select(
    -contains('gtf'), -converted.id, -crovir.transcript
  ) %>%
  distinct()
glimpse(no_ft_miRNA_protein_mRNA_df)

# Save the data
write_parquet(no_ft_miRNA_protein_mRNA_df, 'Data/Merged/No_feature_type_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet')


## Fuse filtered miRNA sample data to mRNA-Protein data ----
# Get shared columns
miRNA_mRNA_prot_intersect2 <- intersect(names(filtered_miRNA_sample_df), names(protein_mRNA_df3))
miRNA_mRNA_prot_intersect2

# Join the data
filt_miRNA_protein_mRNA_df <- full_join(
  filtered_miRNA_sample_df,
  protein_mRNA_df3,
  by = miRNA_mRNA_prot_intersect2,
  # relationship = 'many-to-many'
) %>%
  select(
    -contains('gtf'), -converted.id, -crovir.transcript
  ) %>%
  distinct()
glimpse(filt_miRNA_protein_mRNA_df)

# Save the data
write_parquet(filt_miRNA_protein_mRNA_df, 'Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet')


## Fuse filtered miRNA sample data to mRNA-Protein data ----
# Get shared columns
no_feature_type_intersect2 <- intersect(names(no_ft_filtered_miRNA_sample_df), names(protein_mRNA_df3))
no_feature_type_intersect2

# Join the data
no_ft_filt_miRNA_protein_mRNA_df <- full_join(
  no_ft_filtered_miRNA_sample_df,
  protein_mRNA_df3,
  by = no_feature_type_intersect2,
  # relationship = 'many-to-many'
) %>%
  select(
    -contains('gtf'), -converted.id, -crovir.transcript
  ) %>%
  distinct()
glimpse(no_ft_filt_miRNA_protein_mRNA_df)

# Save the data
write_parquet(no_ft_filt_miRNA_protein_mRNA_df, 'Data/Merged/No_feature_type_filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet')
