# Last Edited: 2025/03/25

# Load Packages ----
library(tidyverse)
library(arrow)
# Source my functions
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')


# Set Up ----

# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read both in as data frames
# Let's keep the analysis limited to the 5UTR
mi_df <- read_parquet(file = miRNA_mRNA_protein_data) %>% 
  dplyr::filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    !sample.id == "CV1082_viridis",
    feature.type == 'five_prime_utr'
  ) %>% 
  select(
    miRNA.cluster, genes, total.energy, total.score
  ) %>% 
  mutate(
    miRNA.cluster = str_replace(miRNA.cluster, 'cvi-', '')
  ) %>% 
  distinct()


# Create color scheme for the venom genes
SVMP_color <- '#4A70B5'
ADAM_color <- '#9A70B5'
SVSP_color <- '#F0B830' 
PLA2_color <- '#7570B3'
miRNA_color <- '#8B0AA5'
VEGF_color <- '#74ADD1'
ohanin_color <- '#3A489C'
myotoxin_color <- '#B2182B'
vQC_color <- '#80BC50'
CRISP_color <- '#E7298A'
CTL_color <- '#F67E17'
EXO_color <- '#005824'
LAAO_color <- '#B35806'
BPP_color <- '#1B9E77'
miRNA_color <- '#BEBEBE' # This is the same as typing 'grey'


# Data Formating ----

## miRNA - Venom gene relationships ----
### Edge table ----
# Create a smaller data frame for only the relavent data for Cytoscape for venom genes and PLA2G2
edge_venom_miRNA_cytoscape_df <- mi_df %>% 
  filter(str_detect(genes, 'Venom_|PLA2G2E.1')) %>% # IMPORTANT: THIS IS HOW YOU DETECT MULTIPLE STINGS!!!!!
  select(
    miRNA.cluster, genes, total.score, total.energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.cluster,
    Gene_ID = genes,
    Binding_Score = total.score,
    Binding_Energy = total.energy
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(edge_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/5UTR/Edge_Cytoscape_miRNA_Venom_interaction_Data.2025.03.25.tsv', sep = '\t', quote = FALSE, row.names = F)

### Note table ----
# Create a node table based on the above edge table
node_venom_miRNA_cytoscape_df <- edge_venom_miRNA_cytoscape_df %>% 
  select(-contains('Binding')) %>% 
  pivot_longer(
    cols = c('miRNA_ID', 'Gene_ID'),
    names_to = 'Locus_Type',
    values_to = 'Locus'
  ) %>% 
  mutate(
    nodeColor = case_when(
      grepl('SVMP', Locus) ~ SVMP_color,
      grepl('SVSP', Locus) ~ SVSP_color,
      grepl('PLA2', Locus) ~ PLA2_color,
      grepl('VEGF', Locus) ~ VEGF_color,
      grepl('ohanin', Locus) ~ ohanin_color,
      grepl('vQC', Locus) ~ vQC_color,
      grepl('CRISP', Locus) ~ CRISP_color,
      grepl('CTL', Locus) ~ CTL_color,
      grepl('EXO', Locus) ~ EXO_color,
      grepl('LAAO', Locus) ~ LAAO_color,
      grepl('myotoxin', Locus) ~ myotoxin_color,
      grepl('BPP', Locus) ~ BPP_color,
      grepl('miRNA', Locus_Type) ~ miRNA_color
    ),
    nodeShape = case_when(
      grepl('Gene_ID', Locus_Type) ~ 'Ellipse',
      grepl('miRNA_ID', Locus_Type) ~ 'Round Rectangle'
    )
  )
# Save the file
write.table(node_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/5UTR/Node_Cytoscape_miRNA_Venom_interaction_Data.2025.03.25.tsv', sep = '\t', quote = FALSE, row.names = F)


## miRNA - Venom gene relationships (Filtered) ----
### Edge table ----
# Create a smaller data frame for only the relavent data for Cytoscape for venom genes and PLA2G2
filtered_edge_venom_miRNA_cytoscape_df <- mi_df %>% 
  filter(
    str_detect(genes, 'Venom_|PLA2G2E.1'),
    # Filter out miRNAs and relationships that don't meet the binding energy and score criteria
    total.score >= 155,
    total.energy <= -7
  ) %>% # IMPORTANT: THIS IS HOW YOU DETECT MULTIPLE STINGS!!!!!
  select(
    miRNA.cluster, genes, total.score, total.energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.cluster,
    Gene_ID = genes,
    Binding_Score = total.score,
    Binding_Energy = total.energy
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(filtered_edge_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/5UTR/Filtered_Edge_Cytoscape_miRNA_Venom_interaction_Data.2025.03.25.tsv', sep = '\t', quote = FALSE, row.names = F)

### Note table ----
# Create a node table based on the above edge table
filtered_node_venom_miRNA_cytoscape_df <- filtered_edge_venom_miRNA_cytoscape_df %>% 
  select(-contains('Binding')) %>% 
  pivot_longer(
    cols = c('miRNA_ID', 'Gene_ID'),
    names_to = 'Locus_Type',
    values_to = 'Locus'
  ) %>% 
  mutate(
    nodeColor = case_when(
      grepl('SVMP', Locus) ~ SVMP_color,
      grepl('SVSP', Locus) ~ SVSP_color,
      grepl('PLA2', Locus) ~ PLA2_color,
      grepl('VEGF', Locus) ~ VEGF_color,
      grepl('ohanin', Locus) ~ ohanin_color,
      grepl('vQC', Locus) ~ vQC_color,
      grepl('CRISP', Locus) ~ CRISP_color,
      grepl('CTL', Locus) ~ CTL_color,
      grepl('EXO', Locus) ~ EXO_color,
      grepl('LAAO', Locus) ~ LAAO_color,
      grepl('myotoxin', Locus) ~ myotoxin_color,
      grepl('BPP', Locus) ~ BPP_color,
      grepl('miRNA', Locus_Type) ~ miRNA_color
    ),
    nodeShape = case_when(
      grepl('Gene_ID', Locus_Type) ~ 'Ellipse',
      grepl('miRNA_ID', Locus_Type) ~ 'Round Rectangle'
    )
  )
# Save the file
write.table(filtered_node_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/5UTR/Filtered_Node_Cytoscape_miRNA_Venom_interaction_Data.2025.03.25.tsv', sep = '\t', quote = FALSE, row.names = F)


