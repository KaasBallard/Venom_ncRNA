# Last Edited: 2024/12/5

#### Set up and Read Data ####

# Load in packages
library(compositions)
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
library(ggplot2)
library(pheatmap)
# Source my functions
source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')

# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")

# Create variable for the fused dataset.
# If this is in Dropbox the file needs to be on the machine
miRNA_mRNA_protein_data <- 'Data/Merged/miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.12.5.tsv'

# Read both in as data frames
miRNA_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T) %>% 
  # Get rid of weirdly annotated genes
  filter(!str_detect(Genes, 'maker-scaffold|augustus|XP_|noncoding')) %>% 
  distinct() %>% 
  # filter(str_detect(Origin, 'three_prime_utr|five_prime_utr')) %>% # Filter out CDS targeting miRNAs
  dplyr::select(-contains('Blast'), -E.value, -Bit.Score, -miRNA.Identity.Type) %>%  # Remove columns to save memory
  dplyr::rename(
    'miRNA.Cluster.Original' = 'miRNA.Cluster',
    'miRNA.Cluster' = 'Putative.miRNA.Name'
  )

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  # Remove target info to reduce RAM use
  dplyr::select(
    -contains('Length'), -contains('Start'), -contains('End'), -contains('Strand'), -miRNA.Sequence,
    -contains('Max'), -miRNA.Cluster.Original, -Positions, -Origin
  ) %>% 
  distinct()
rm(miRNA_mRNA_protein_df)
gc()



# Create color scheme for the venom genes
SVMP_color <- '#4A70B5'
ADAM_color <- '#9A70B5'
SVSP_color <- '#F0B830' 
PLA2_color <- '#7570B3'
VEGF_color <- '#74ADD1'
ohanin_color <- '#3A489C'
myotoxin_color <- '#B2182B'
vQC_color <- '#80BC50'
CRISP_color <- '#E7298A'
CTL_color <- '#F67E17'
EXO_color <- '#005824'
LAAO_color <- '#B35806'
BPP_color <- '#1B9E77'
other_color <- '#666666'
miRNA_color <- 'grey14'
novel_miRNA_color <- 'grey'
three_prime_color <- 'black'
five_prime_color <- '#1B9E77'
cds_color <- '#4A70B5'

# Create an order for the Genes
venom_gene_order <- c(
  'BPP',
  'myotoxin',
  'ohanin',
  'PLA2B1', 'PLA2K', 'PLA2C1',
  'SVSP1', 'SVSP2', 'SVSP3', 'SVSP10', 'SVSP6', 'SVSP11', 'SVSP7', 'SVSP8', 'SVSP9',
  'SVMP1', 'SVMP2', 'SVMP3', 'SVMP4', 'SVMP5', 'SVMP6', 'SVMP7', 'SVMP8', 'SVMP9', 'SVMP10', 'SVMP11', 'SVMP12',
  'CRISP1', 'CRISP2', 'CRISP3', 'CRISP4',
  'CTL1', 'CTL2', 'CTL3', 'CTL4', 'CTL5', 'CTL6',
  'EXO1', 'EXO2', 'EXO3',
  'LAAO1', 'LAAO2', 'LAAO3',
  'VEGF1', 'VEGF2',
  'vQC1', 'vQC2'
) 

# Create color scheme for the venom genes
venom_colors <- c(
  SVMP = SVMP_color,
  ADAM = ADAM_color,
  SVSP = SVSP_color,
  PLA2 = PLA2_color,
  miRNA = miRNA_color,
  VEGF = VEGF_color,
  Ohanin = ohanin_color,
  Myotoxin = myotoxin_color,
  vQC = vQC_color,
  CRISP = CRISP_color,
  CTL = CTL_color,
  EXO = EXO_color,
  LAAO = LAAO_color,
  BPP = BPP_color,
  others = other_color
)

# Create color scheme for the 
loci_colors <- c(
  SVMP = SVMP_color,
  ADAM = ADAM_color,
  SVSP = SVSP_color,
  PLA2 = PLA2_color,
  miRNA = miRNA_color,
  VEGF = VEGF_color,
  Ohanin = ohanin_color,
  Myotoxin = myotoxin_color,
  vQC = vQC_color, 
  CRISP = CRISP_color,
  CTL = CTL_color,
  EXO = EXO_color,
  LAAO = LAAO_color,
  BPP = BPP_color,
  others = other_color,
  novel_miRNA = novel_miRNA_color,
  miRNA = miRNA_color
)



#### Format data ####

# Create data frame for formating
reduced_mi_df <- mi_df %>% 
  # Filter out mRNAs that weren't in the library
  # filter(
  #   In.Library == 'Yes',
  #   !str_starts(Genes, 'ADAM')
  # ) %>%
  dplyr::select(
    Sample.ID,
    Genes,
    Venom.Family,
    mRNA.VST,
    Intensity,
    miRNA.Cluster,
    miRNA.VST,
    miRNA.RPM,
  ) %>% 
  distinct() %>% 
  mutate(
    Species = case_when(
      grepl('.viridis', Sample.ID) ~ 'C.viridis',
      grepl('.concolor', Sample.ID) ~ 'C.concolor',
      grepl('.lutosus', Sample.ID) ~ 'C. lutosus'
    )
  )
glimpse(reduced_mi_df)
gc()


# # Check the number of targeted Genes
# target_genes <- binding_mrna_mirna_prot_df %>% filter(Total.Energy < -7, Total.Score > 155) %>% distinct(Genes)
# glimpse(target_genes)



#### Does the number of miRNAs that target a gene correlate with mRNA expression level ####

# Format data to get number of miRNAs per gene
number_of_mirnas_per_gene_df <- reduced_mi_df %>% 
  distinct(miRNA.Cluster, Genes) %>% 
  group_by(Genes) %>% 
  # dplyr::summarise(
  #   number.of.miRNAs = n()
  # ) %>% 
  dplyr::summarise(
    number.of.miRNAs = sum(!is.na(miRNA.Cluster))
  )

# Remove miRNA expression levels and fuse number of miRNAs to the main data for plot
mirna_count_gene_exp_df <- reduced_mi_df %>% 
  select(-miRNA.Cluster, -miRNA.VST, -miRNA.RPM) %>% 
  distinct() %>% 
  left_join(number_of_mirnas_per_gene_df, by = c('Genes')) %>% 
  # Filter out lowly expressed non-venom genes
  filter(!(Venom.Family == 'others' & Intensity == 0))
  
# Forget and clear unnecessary data
rm(number_of_mirnas_per_gene_df)
gc()



# Plot miRNA counts against mRNA expression data
mirna_count_vs_mrna_expression <- ggplot(mirna_count_gene_exp_df, aes(x = number.of.miRNAs, y = mRNA.VST, color = Venom.Family)) +
  geom_point(
    aes(shape = Species),
    size = 2.5
  ) + 
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = T,
    color = 'black',
    linetype = 'solid',
    linewidth = 0.8
  ) +
  stat_poly_eq(
    inherit.aes = F,
    aes(
      x = number.of.miRNAs, y = mRNA.VST,
      label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")
    ),
    formula = y ~ x
  ) +
  scale_x_continuous() +
  scale_y_continuous() +
  scale_color_manual(values = venom_colors) +
  labs(
    title = 'The number of miRNAs per gene correlates with mRNA expression level',
    x = 'Number of miRNA per gene',
    y = 'mRNA expression level',
    color = 'Gene Family'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', hjust = 0),
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5)
  )
mirna_count_vs_mrna_expression
ggsave(
  filename = 'Figures/Expression_Plots/Qualitative_Analysis/Dotplots/All_expressed_proteins_miRNA_numbers_vs_mRNA_expression_plot_2024.12.6.pdf', 
  plot = mirna_count_vs_mrna_expression,
  height = 15, width = 20,
  create.dir = T
)
  


#### Does the expression of an miRNA scale with the number of genes it targets ####


# Create data frame for formating
reduced_mi_df2 <- mi_df %>% 
  # Filter out mRNAs that weren't in the library
  # filter(
  #   In.Library == 'Yes',
  #   !str_starts(Genes, 'ADAM')
  # ) %>%
  filter(
    Total.Score >= 155,
    Total.Energy <= -7
  ) %>% 
  dplyr::select(
    Sample.ID,
    Genes,
    Venom.Family,
    mRNA.VST,
    Intensity,
    miRNA.Cluster,
    miRNA.VST,
    miRNA.RPM,
  ) %>% 
  distinct() %>% 
  mutate(
    Species = case_when(
      grepl('.viridis', Sample.ID) ~ 'C.viridis',
      grepl('.concolor', Sample.ID) ~ 'C.concolor',
      grepl('.lutosus', Sample.ID) ~ 'C. lutosus'
    )
  )
glimpse(reduced_mi_df)


# Get the number of genes per miRNA
genes_per_miRNA_df <- reduced_mi_df2 %>% 
  distinct(
    miRNA.Cluster, Genes
  ) %>% 
  group_by(miRNA.Cluster) %>% 
  dplyr::summarise(
    number.of.Genes = n()
  ) %>% 
  filter(!is.na(miRNA.Cluster))
glimpse(genes_per_miRNA_df)

# Fuse genes per miRNA to the main data with miRNA expression levels
mirna_count_mirna_exp_df <- reduced_mi_df2 %>% 
  # Filter out lowly expressed non-venom genes
  filter(!(Venom.Family == 'others' & Intensity == 0)) %>% 
  select(-Genes, -Venom.Family, -contains('mRNA'), -Intensity) %>% 
  distinct() %>% 
  left_join(genes_per_miRNA_df, by = c('miRNA.Cluster')) %>% 
  filter(!is.na(miRNA.Cluster))
glimpse(mirna_count_mirna_exp_df)

# Forget and clear unnecessary data
rm(genes_per_miRNA_df)
gc()

# Plot miRNA counts against mRNA expression data
mirna_count_vs_mirna_expression <- ggplot(mirna_count_mirna_exp_df, aes(x = number.of.Genes, y = miRNA.VST)) +
  geom_point(
    aes(shape = Species),
    size = 2.5
  ) + 
  geom_smooth(
    method = 'lm',
    formula = y ~ x,
    se = T,
    color = 'black',
    linetype = 'solid',
    linewidth = 0.8
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(
    title = 'The number of genes per miRNA correlates with miRNA expression level',
    x = 'Number of genes per miRNA',
    y = 'miRNA expression level'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', hjust = 0),
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5)
  )
mirna_count_vs_mirna_expression

  