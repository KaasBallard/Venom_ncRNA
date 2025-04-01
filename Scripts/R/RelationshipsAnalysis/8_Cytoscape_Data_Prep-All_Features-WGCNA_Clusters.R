# Last Edited: 2025/03/23

# Load Packages ----
library(tidyverse)
library(arrow)
# Source my functions
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')


# Set Up ----

# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read in the cluster data
wgcn_clusters <- 'Data/WGCNA/WGCNA_Loci_Clusters_2025.03.26.parquet'

# Read both in as data frames
mi_df <- read_parquet(file = miRNA_mRNA_protein_data) %>% 
  dplyr::filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    !sample.id == "CV1082_viridis",
    # Filter out Venom_ADAM28 and CTL6 because they were not in the WGCNA data
    !str_detect(genes, 'Venom_ADAM28|Venom_CTL6')
  ) %>% 
  select(
    miRNA.cluster, genes, total.energy, total.score, feature.type
  ) %>% 
  mutate(
    miRNA.cluster = str_replace(miRNA.cluster, 'cvi-', '')
  ) %>% 
  distinct()

# Read in the wgcn clusters
wgcna_df <- read_parquet(file = wgcn_clusters) %>% 
  # Change the name of the genes col
  rename(Locus = genes) %>% 
  mutate(
    # Remove Venom_ from the venom gene names
    Locus = str_replace(Locus, 'Venom_', ''),
    # Remove 'cvi-' from infront of the miRNA names
    Locus = str_replace(Locus, 'cvi-', ''),
    # Remove the space between the cluster and number
    Locus = str_replace_all(Locus, '_', ' ')
  )


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
  filter(str_detect(genes, 'Venom_')) %>% # I am not doing PLA2GE.1 because it wasn't in the WGCNA data
  select(
    miRNA.cluster, genes, feature.type
  ) %>% 
  rename(
    miRNA_ID = miRNA.cluster,
    Gene_ID = genes,
    Target = feature.type
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(edge_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/WGCNA/Edge_Cytoscape_miRNA_Venom_interaction_Data.2025.03.26.tsv', sep = '\t', quote = FALSE, row.names = F)

### Note table ----
# Create a node table based on the above edge table
node_venom_miRNA_cytoscape_df <- edge_venom_miRNA_cytoscape_df %>% 
  select(-Target) %>% 
  pivot_longer(
    cols = c('miRNA_ID', 'Gene_ID'),
    names_to = 'Locus_Type',
    values_to = 'Locus'
  ) %>% 
  distinct() %>% 
  full_join(
    wgcna_df,
    by = c('Locus')
  ) %>% 
  mutate(
    ClusterColor = case_when(
      str_detect(module, 'turquoise') ~ '#40E0D0',
      str_detect(module, 'blue') ~ '#0000FF',
      str_detect(module, 'brown') ~ '#A52A2A',
      str_detect(module, 'grey') ~ '#C0C0C0'
    ),
    NodeShape = case_when(
      grepl('Gene_ID', Locus_Type) ~ 'Ellipse',
      grepl('miRNA_ID', Locus_Type) ~ 'Round Rectangle'
    )
  )
# Save the file
write.table(node_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/WGCNA/Node_Cytoscape_miRNA_Venom_interaction_Data.2025.03.26.tsv', sep = '\t', quote = FALSE, row.names = F)

