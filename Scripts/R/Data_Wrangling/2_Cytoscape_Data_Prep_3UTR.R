# Last Edited: 2024/10/22

#### Load Packages ####
library(tidyverse)
# Source my functions
source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')


#### Set Up ####

# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")

# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.08.31.tsv'


# Read both in as data frames
miRNA_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T)

# # I am going to exclude the following genes because they are not well annotated
# excluded_genes = c(
#   'maker-scaffold-mi1-augustus-gene-59.13_crovir-transcript-12940',
#   'maker-scaffold-mi1-augustus-gene-59.20_crovir-transcript-12947',
#   'maker-scaffold-mi2-augustus-gene-22.17_crovir-transcript-739',
#   'XP_016876419',
#   'XP_011528471'
# )

# Let's keep the analysis limited to the 3UTR
# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  select(-miRNA.Cluster) %>%
  dplyr::rename(
    'Genes' = 'Converted.Gene.IDs',
      'miRNA.Cluster' = 'Putative.miRNA.Name'
  ) %>%
  select(miRNA.Cluster, everything()) %>% # Move the new miRNA.Clusters to the front
  filter(!str_detect(Genes, 'maker-scaffold|augustus|XP_')) %>% # This should get rid of all of the weirdly annotated genes
  # filter(!(Genes %in% excluded_genes)) %>% 
  filter(Origin == 'three_prime_utr') %>% # Filter out everything not targeting the 3' UTR
  select(-Origin) %>% # Remove to save memory
  mutate(Total.Energy = abs(Total.Energy)) %>% 
  mutate(miRNA.Cluster = gsub('cvi-', '', miRNA.Cluster)) %>% 
  # Change the name of the venom adam genes so that they are correct
  mutate(
    Genes = ifelse(Genes == "Venom_ADAM28_1", 'Venom_SVMP12', Genes),
    Genes = ifelse(Genes == 'Venom_ADAM28_2', 'ADAM28', Genes)
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
EXO_color <- '#49FFFF'
LAAO_color <- '#B35806'
BPP_color <- '#1B9E77'
miRNA_color <- '#BEBEBE' # This is the same as typing 'grey'




#### Data Formating ####

# Create a smaller data frame for only the relavent data for Cytoscape for venom genes and PLA2G2
edge_venom_miRNA_cytoscape_df <- mi_df %>% 
  filter(str_detect(Genes, 'Venom|PLA2G2E.1')) %>% # IMPORTANT: THIS IS HOW YOU DETECT MULTIPLE STINGS!!!!!
  select(
    miRNA.Cluster, Genes, Total.Score, Total.Energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.Cluster,
    Gene_ID = Genes,
    Binding_Score = Total.Score,
    Binding_Energy = Total.Energy
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(edge_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/3UTR/Edge_Cytoscape_miRNA_Venom_interaction_Data.2024.10.22.tsv', sep = '\t', quote = FALSE, row.names = F)

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
write.table(node_venom_miRNA_cytoscape_df, file = 'Data/Cytoscape/3UTR/Node_Cytoscape_miRNA_Venom_interaction_Data.2024.10.22.tsv', sep = '\t', quote = FALSE, row.names = F)

# Create a smaller data frame for only PLA2, SVMP, and SVSP
classic_families_cytoscape_df <- mi_df %>% 
  filter(str_detect(Genes, 'SVMP|SVSP|PLA')) %>% # IMPORTANT: THIS IS HOW YOU DETECT MULTIPLE STINGS!!!!!
  select(
    miRNA.Cluster, Genes, Total.Score, Total.Energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.Cluster,
    Gene_ID = Genes,
    Binding_Score = Total.Score,
    Binding_Energy = Total.Energy
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(classic_families_cytoscape_df, file = 'Data/Cytoscape/3UTR/Edge_Cytoscape_miRNA_Classic_Families_interaction_Data.2024.10.22.tsv', sep = '\t', quote = FALSE, row.names = F)


# Create data for cytoscape interaction network for SVMPs
SVMP_miRNA_cytoscape_df <- mi_df %>% 
  filter(str_detect(Genes, 'SVMP')) %>% 
  select(
    miRNA.Cluster, Genes, Total.Score, Total.Energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.Cluster,
    Gene_ID = Genes,
    Binding_Score = Total.Score,
    Binding_Energy = Total.Energy
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(SVMP_miRNA_cytoscape_df, file = 'Data/Cytoscape/3UTR/Edge_Cytoscape_miRNA_SVMP_interaction_Data.2024.10.22.tsv', sep = '\t', quote = FALSE, row.names = F)

# Create data for cytoscape interaction network for SVSPs
SVSP_miRNA_cytoscape_df <- mi_df %>% 
  filter(str_detect(Genes, 'SVSP')) %>% 
  select(
    miRNA.Cluster, Genes, Total.Score, Total.Energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.Cluster,
    Gene_ID = Genes,
    Binding_Score = Total.Score,
    Binding_Energy = Total.Energy
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(SVSP_miRNA_cytoscape_df, file = 'Data/Cytoscape/3UTR/Edge_Cytoscape_miRNA_SVSP_interaction_Data.2024.10.22.tsv', sep = '\t', quote = FALSE, row.names = F)


# Create data for cytoscape interaction network for PLA2s
PLA2_miRNA_cytoscape_df <- mi_df %>% 
  filter(str_detect(Genes, 'PLA2')) %>% 
  select(
    miRNA.Cluster, Genes, Total.Score, Total.Energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.Cluster,
    Gene_ID = Genes,
    Binding_Score = Total.Score,
    Binding_Energy = Total.Energy
  ) %>% 
  mutate(Gene_ID = str_replace_all(Gene_ID, 'Venom_', '')) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(PLA2_miRNA_cytoscape_df, file = 'Data/Cytoscape/3UTR/Edge_Cytoscape_miRNA_PLA2_interaction_Data.2024.10.22.tsv', sep = '\t', quote = FALSE, row.names = F)


# Create a smaller data frame for only the relavent data for Cytoscape for non-venom genes
nonvenom_miRNA_cytoscape_df <- mi_df %>% 
  filter(!str_detect(Genes, 'Venom')) %>% 
  select(
    miRNA.Cluster, Genes, Total.Score, Total.Energy
  ) %>% 
  rename(
    miRNA_ID = miRNA.Cluster,
    Gene_ID = Genes,
    Binding_Score = Total.Score,
    Binding_Energy = Total.Energy
  ) %>% 
  mutate(miRNA_ID = str_replace_all(miRNA_ID, '_', ' '))
# Save the file
write.table(nonvenom_miRNA_cytoscape_df, file = 'Data/Cytoscape/3UTR/Edge_Cytoscape_miRNA_Non-Venom_interaction_Data.2024.10.22.tsv', sep = '\t', quote = FALSE, row.names = F)

# Sample edge data
edges <- data.frame(
  source = c("node1", "node2", "node3", "node4"),
  target = c("node2", "node3", "node4", "node1"),
  stringsAsFactors = FALSE
)

# Extract unique nodes
unique_nodes <- unique(c(edges$source, edges$target))

# Create node data with node shapes
nodes <- data.frame(
  id = unique_nodes,
  type = c("enzyme", "structural", "enzyme", "structural"), # Example types
  nodeShape = c("diamond", "hexagon", "diamond", "hexagon"), # Corresponding shapes
  stringsAsFactors = FALSE
)



