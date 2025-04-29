# Last Edited: 2024/12/7

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
miRNA_mRNA_protein_data <- 'Data/Merged/Venom_miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.12.5.tsv'

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

# Remove columns that are unnecessary and filter out non-venom genes
base_data_df <- mi_df %>% 
  # Filter anything out that isn't a venom gene
  filter(
    # Total.Score > 155, 
    # Total.Energy < 7,
    !is.na(miRNA.Cluster),
    str_detect(Genes, 'Venom_'),
    Protein.Observed == 'Yes'
  ) %>% 
  select(
    Sample.ID, Genes, Venom.Family, mRNA.Counts, Intensity, miRNA.Cluster, miRNA.RPM
  ) %>% 
  distinct()

# Add up all the miRNA expression per gene
total_mirna_df <- base_data_df %>% 
  group_by(Genes, Sample.ID) %>%
  summarise(
    Total.miRNA.Exp = sum(miRNA.RPM, na.rm = T), 
  )

# Join the miRNA total expression back in
prot_mrna_mirna_df <- left_join(
  base_data_df,
  total_mirna_df
) %>% 
  select(-miRNA.Cluster, -miRNA.RPM) %>% 
  distinct() %>% 
  # Perform clr transformation
  group_by(Sample.ID) %>%
  mutate(Scaled.Protein.Expression = Intensity / sum(Intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.mRNA.Expression = mRNA.Counts / sum(mRNA.Counts)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  # Perform CLR transformation
  mutate(CLR.Scaled.Protein.Expression = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(CLR.Scaled.mRNA.Expression = as.numeric(compositions::clr(Scaled.mRNA.Expression)))




#### Create total miRNA expression plot ####

# Create expression graph
mirna_prot_mrna_gene_plot <- ggplot(prot_mrna_mirna_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Total.miRNA.Exp, size = Total.miRNA.Exp)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(
    low = 'blue', high = 'red'
  ) +
  scale_size_continuous(range = c(1, 6)) +
  scale_x_continuous(breaks = seq(-15, 15, by = 1.5)) +  # Set x-axis to continuous scale
  scale_y_continuous(breaks = seq(-15, 15, by = 1.5)) +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs protein expression',
    x = 'mRNA expression (CLR)',
    y = 'Protein expression (CLR)',
    color = 'miRNA expression (total reads/million)',
    size = 'miRNA expression (total reads/million)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12, hjust = 0.5),
    legend.position = "bottom",
    legend.title.position = 'top'
  ) +
  facet_wrap( ~ Genes)
mirna_prot_mrna_gene_plot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/miRNA_total_expression_mRNA_vs_protein_expression_by_genes_2024.12.7.pdf', plot = mirna_prot_mrna_gene_plot, height = 15, width = 25, create.dir = T)


# Create expression graph for families
mirna_prot_mrna_family_plot <- ggplot(prot_mrna_mirna_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Total.miRNA.Exp, size = Total.miRNA.Exp)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(
    low = 'blue', high = 'red'
  ) +
  scale_size_continuous(range = c(1, 6)) +
  scale_x_continuous(breaks = seq(-15, 15, by = 1.5)) +  # Set x-axis to continuous scale
  scale_y_continuous(breaks = seq(-15, 15, by = 1.5)) +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs protein expression',
    x = 'mRNA expression (CLR)',
    y = 'Protein expression (CLR)',
    color = 'miRNA expression (total reads/million)',
    size = 'miRNA expression (total reads/million)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12, hjust = 0.5),
    legend.position = "bottom",
    legend.title.position = 'top'
  ) +
  facet_wrap( ~ Venom.Family)
mirna_prot_mrna_family_plot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/miRNA_total_expression_mRNA_vs_protein_expression_by_family_2024.12.7.pdf', plot = mirna_prot_mrna_family_plot, height = 15, width = 25, create.dir = T)

