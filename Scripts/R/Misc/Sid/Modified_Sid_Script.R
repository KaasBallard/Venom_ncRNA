#### Set up ####

# Load in packages
library(cowplot)
library(tidyverse)
library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(readxl)
library(ggpmisc)

# Set working directory
setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')

# Define file paths for different data sets
rnaseq_data <- 'Data/mRNA/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt'
protein_data <- 'Data/Protein/2023_Aug9_QEHF_Castoe.xlsx'
conversion_table <- 'Data/Conversion_Table/Cvv_GTF_to_converted_names_05.22.23.txt'
miRNA_data <- 'Data/Merged/Merged_Column_and_Row_Filtered.tsv'


#### Protein Data ####

# Read protein data from Excel file
protein_df <- read_excel(protein_data, sheet = 1) # version 1 old anthony

# Create a dataframe to map Anthony's protein names to standard names
# Convert anthony names to standard names (version 1)
name_conversion <- data.frame(
  stringsAsFactors = FALSE,
  Anthony_name = c("C_o_cerberus_AZ_Kingman",
                   "C_o_cerberus_AZ_NW_Kingman","C_o_concolor_UT_Duchesne",
                   "C_o_lutosus_528_UT_Beaver_Co","C_o_lutosus_UT_Beaver_Co",
                   "C_o_tigris_AZ_Pima_Co","C_v_viridis_CO_Hwy93",
                   "C_v_viridis_CO_Otero_Co","C_v_viridis_CO_Weld_Co",
                   "C_v_viridis_MT","C_v_viridis_NM_6","C_v_viridis_NM_Luna_Co",
                   "C_v_viridis_NM_Socorro","C_v_viridis_CO_San_Miguel_Co"),
  New_name = c("NONE",
               "CV1090_cerberus_Other_M","CV0985_concolor_Other_F","CV0987_lutosus_Other_F",
               "NONE","NONE","CV1096_viridis_North_F",
               "CV1084_viridis_Mid_M","CV1087_viridis_North_F",
               "CV0857_viridis_North_M","NONE","CV1086_viridis_South_M",
               "CV1089_viridis_South_M","CV1081_viridis_Mid_M")
)

# Process protein data: filter, rename, and remove certain columns
filtered_protein_df <- protein_df %>%
  select(-51) %>% # remove indistinguishable proteins column # 51
  pivot_longer(cols = 9:50, names_sep = '\\s', names_to = c("SampleID", "Feature"), values_to = "Value") %>% #9:50
  left_join(name_conversion, by = c('SampleID'='Anthony_name')) %>%
  filter(!New_name == 'NONE') %>% # remove samples for which we have no study samples
  select(-SampleID) %>%
  rename(SampleID = New_name) %>%
  filter(!is.na(Protein))
  # It seems that CV1082 doesn't have any data for this part. I also removed it for the RNAseq data.

# Read and preprocess protein conversion table
protein_to_orthos <- read.table(conversion_table, header=T) 
protein_to_orthos <- protein_to_orthos %>%
  filter(!grepl('fgenesh', gtf_gene)) %>% # remove crappy fgenesh annotations
  select(crovir_transcript, converted_id_no_dups) %>%
  mutate(crovir_transcript = gsub('transcript', 'protein', crovir_transcript)) %>%
  mutate(crovir_transcript = ifelse(converted_id_no_dups == 'Venom_myotoxin', 'crovir-protein-myotoxin', crovir_transcript))

# Merge protein data with conversion table, filter for Venoms
filtered_protein_df <- filtered_protein_df %>%
  left_join(protein_to_orthos, by = c('Protein' = 'crovir_transcript')) %>%
  select(-Protein) %>%
  rename(Protein = converted_id_no_dups) %>%
  filter(grepl('Venom', Protein)) # keep only Venoms


#### RNA Data ####

# RNAseq VST normalized count matrix
rna_df <- read.table(rnaseq_data, header=T)

# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
rna_df <- rna_df %>%
  select(matches('LVG|RVG')) %>%
  filter(grepl('Venom', row.names(.)))


column_numbers <- 1:13
target_names <- paste0('VG_', column_numbers)

for (i in column_numbers) {
  target_columns <- grep(paste0('_', i, '$'), colnames(rna_df))
  
  if (length(target_columns) > 0) {
    if (length(target_columns) == 1) {
      # If only one column is found, copy it to become VG_i
      rna_df[target_names[i]] <- rna_df[, target_columns]
    } else {
      # Calculate row mean for multiple columns
      rna_df[target_names[i]] <- rowMeans(rna_df[, target_columns], na.rm = TRUE)
    }
  }
}


#### RNA Data Filtering ####

rna_df_filtered <- rna_df %>%
  select(-contains('LVG'),-contains('RVG'))

new_names <- data.frame(
  stringsAsFactors = FALSE,
  VG_name = c("VG_1","VG_2","VG_3","VG_4","VG_5","VG_6",
              "VG_7","VG_8","VG_9","VG_11","VG_12",
              "VG_13"),
  New_name = c("CV1084_viridis_Mid_M","CV1081_viridis_Mid_M",
               "CV1089_viridis_South_M",
               "CV0857_viridis_North_M","CV1087_viridis_North_F",
               "CV0987_lutosus_Other_F","CV0985_concolor_Other_F",
               "CV1090_cerberus_Other_M",
               "CV1086_viridis_South_M","CV1095_viridis_North_M",
               "CV1082_viridis_South_M","CV1096_viridis_North_F")
)


rna_df_filtered <- rna_df_filtered %>%
  filter(str_detect(row.names(.), 'Venom')) %>%
  rownames_to_column(var = "Toxin") %>%
  mutate(Toxin = gsub("CRISP_", 'CRISP', Toxin)) %>%
  mutate(Toxin = gsub("CTL_", 'CTL', Toxin)) %>%
  filter(str_detect(Toxin, 'ADAM|CRISP3|CRISP4|CTL1|CTL4|CTL5|CTL6|EXO|LAAO1|LAAO2|VEGF2|SVMP11', negate = T)) %>% # remove lowly expressed genes in all samples
  pivot_longer(cols = VG_1:VG_13,
               names_to = "VG_Sample",
               values_to = "Average_Exp") %>%
  left_join(new_names, by = c('VG_Sample'='VG_name')) %>%
  select(-VG_Sample) %>%
  rename(VG_Sample = New_name)

rna_df_filtered <- rna_df_filtered %>%
  filter(VG_Sample %in% c('CV1087_viridis_North_F', 'CV0985_concolor_Other_F', 'CV0987_lutosus_Other_F'))

# Add GEX to prot
protein_GEX <- filtered_protein_df %>%
  left_join(rna_df_filtered, by = c('SampleID'='VG_Sample', 'Protein'='Toxin')) %>%
  separate(SampleID, into = c('CV_ID','Lineage', 'Population','Sex'), sep = '_', remove = F) %>%
  mutate(Venom_family = case_when(grepl('SVMP', Protein) ~ 'SVMP',
                                  grepl('SVSP', Protein) ~ 'SVSP',
                                  grepl('PLA2', Protein) ~ 'PLA2',
                                  grepl('ADAM', Protein) ~ 'ADAM',
                                  grepl('CRISP', Protein) ~ 'CRISP',
                                  #grepl('CTL', Protein) ~ 'CTL',
                                  grepl('EXO', Protein) ~ 'EXO',
                                  grepl('LAAO', Protein) ~ 'LAAO',
                                  grepl('myotoxin', Protein) ~ 'myotoxin',
                                  grepl('BPP', Protein) ~ 'BPP',
                                  TRUE ~ 'others'))


#### miRNA Data ####

# Read the miRNA dataframe I created from the shortstack counts.txt, miRanda tab output, and bedtools intersect between the miRanda tab output and the genome gtf.
miRNA_df <- read.table(file = miRNA_data, sep = '\t', header = T)

# Read a second name conversion df in so that I can fuse and convert easily.
conversion_df <- read.table(file = conversion_table, header = T)

#### miRNA Data Curation ####

# Filter out unnecessary columns and get rid of bad fgenesh anotations
filtered_miRNA_df <- miRNA_df %>%
  select(-Feature.type, -miRNA.Sequence.for.Counts..Hairpin., -Unnamed..0, -X, -Genome.Chrom, -Genome.Start, -Genome.End) %>%
  rename(miRNA.Target.Strandedness = Genome.Strandedness, gtf_gene = gene_id) %>% # Why the hell do I have to list the new column name first? This is in my opinion, terrible design. R is bad.
  filter(!grepl('fgenesh', gtf_gene)) %>%
  distinct()

#### Fuse miRNA data frame to the other dataframes ####

# Left_join the conversion_df to the filtered miRNA_df.
miRNA_with_ids_df <- left_join(filtered_miRNA_df, conversion_df, by = "gtf_gene") %>% select(-converted_id)





#### Gene expression and Protein Amount ####

# According to Anthony, best to use unscaled Intensity
gene_expression_figure <- protein_GEX %>%
  filter(Feature == 'Intensity') %>%
  filter(Average_Exp > 10) %>%
  ggplot() +
  geom_point(aes(x=Average_Exp, y=Value, color = Venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    aes(x = Average_Exp, y = Value),
    method = "lm",
    se = FALSE,  # Do not display confidence intervals
    color = "black",
    linetype = "dashed",
    # formula = y ~ poly(x, 3, raw = TRUE) # 3rd degree polynomial
    formula = y ~ x # linear
  ) +
  stat_poly_eq(aes(x = Average_Exp, y = Value),
               formula = y ~ x) +
  scale_color_brewer(palette = 'Dark2') +
  xlab('Gene expression') +
  ylab('Peak intensity') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 3))
