# Last Edited: 2025/01/29

# Set up and Read Data ----

## Read packages ----
library(tidyverse)
library(arrow)
library(viridis)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(scales)

## Read data ----

# Create variables for the data
miRNA_number_data <- 'Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filtered_miRNA_number_data <- 'Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'

# Read data
miRNA_num_df <- read_parquet(file = miRNA_number_data) %>% 
  filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    !str_detect(genes, 'ADAM28')
  ) %>% 
  select(
    sample.id, miRNA.target.chrom, genes, venom.family, feature.type, number.of.miRNAs
  ) %>% 
  distinct() %>% 
  mutate(
    genes = str_remove(genes, '^Venom_'),
    sample.id = case_when(
      grepl('CV0857_viridis', sample.id) ~ 'CV0857',
      grepl('CV0985_concolor', sample.id) ~ 'CV0985',
      grepl('CV0987_lutosus', sample.id) ~ 'CV0987',
      grepl('CV1081_viridis', sample.id) ~ 'CV1081',
      grepl('CV1086_viridis', sample.id) ~ 'CV1086',
      grepl('CV1087_viridis', sample.id) ~ 'CV1087'
    )
  ) # Remove 'Venom_' from genes
glimpse(miRNA_num_df)

# Read the filtered data
filt_miRNA_num_df <- read_parquet(file = filtered_miRNA_number_data) %>% 
  filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    !str_detect(genes, 'ADAM28'),
    !is.na(feature.type)
  ) %>% 
  select(
    sample.id, miRNA.target.chrom, genes, venom.family, feature.type, number.of.miRNAs
  ) %>% 
  distinct() %>% 
  mutate(
    genes = str_remove(genes, '^Venom_'),
    sample.id = case_when(
      grepl('CV0857_viridis', sample.id) ~ 'CV0857',
      grepl('CV0985_concolor', sample.id) ~ 'CV0985',
      grepl('CV0987_lutosus', sample.id) ~ 'CV0987',
      grepl('CV1081_viridis', sample.id) ~ 'CV1081',
      grepl('CV1086_viridis', sample.id) ~ 'CV1086',
      grepl('CV1087_viridis', sample.id) ~ 'CV1087'
    )
  ) # Remove 'Venom_' from genes
glimpse(filt_miRNA_num_df)

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
other_color <- '#666666'
three_prime_color <- 'black'
five_prime_color <- '#0072b2'
cds_color <- '#d55e00'
viridis_color <- '#2D8A5C'
lutosus_color <- '#8C4720'
concolor_color <- '#E0A229'

# Create an order for the genes
venom_gene_order <- c(
  'BPP',
  'myotoxin',
  'ohanin',
  'PLA2B1', 'PLA2K', 'PLA2C1', 'PLA2A1',
  'SVSP1', 'SVSP2', 'SVSP3', 'SVSP4', 'SVSP10', 'SVSP5', 'SVSP6', 'SVSP11', 'SVSP7', 'SVSP8', 'SVSP9',
  'SVMP1', 'SVMP2', 'SVMP3', 'SVMP4', 'SVMP5', 'SVMP6', 'SVMP7', 'SVMP8', 'SVMP9', 'SVMP10', 'SVMP11', 'SVMP12',
  'CRISP1', 'CRISP2', 'CRISP3', 'CRISP4',
  'CTL1', 'CTL2', 'CTL3', 'CTL4', 'CTL5', 'CTL6',
  'EXO1', 'EXO2', 'EXO3',
  'LAAO1', 'LAAO2', 'LAAO3',
  'VEGF1', 'VEGF2',
  'vQC1', 'vQC2'
) 

# Set the order of venom genes
miRNA_num_df$genes <- factor(miRNA_num_df$genes, levels = venom_gene_order)

# Set the order of venom genes
filt_miRNA_num_df$genes <- factor(filt_miRNA_num_df$genes, levels = venom_gene_order)


# Make figures ----

# Create bar plot for the unfiltered data
miRNA_number_bar_plot <- ggplot(miRNA_num_df, aes(x = genes, y = number.of.miRNAs, fill = sample.id)) +
  geom_bar(position = 'fill', stat = 'identity', color = 'black') +
  facet_wrap( ~ feature.type) +
  labs(
    x = 'Number of miRNAs per gene',
    y = 'Genes',
    fill = 'Sample ID'
  ) +
  scale_fill_brewer(palette = 'Spectral') +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 13)
  ) +
  geom_text(aes(label = number.of.miRNAs), position = position_fill(vjust = 0.5), size = 2)
miRNA_number_bar_plot
ggsave('Figures/Bar_Graphs/miRNA_Numbers/Number_of_miRNAs_per_gene_2025.02.02.pdf', plot = miRNA_number_bar_plot, create.dir = TRUE, width = 16, height = 10)


# Create bar plot for the unfiltered data
filt_miRNA_number_bar_plot <- ggplot(filt_miRNA_num_df, aes(x = genes, y = number.of.miRNAs, fill = sample.id)) +
  geom_bar(position = 'fill', stat = 'identity', color = 'black') +
  facet_wrap( ~ feature.type) +
  labs(
    x = 'Number of miRNAs per gene',
    y = 'Genes',
    fill = 'Sample ID'
  ) +
  scale_fill_brewer(palette = 'Spectral') +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 13)
  ) +
  geom_text(aes(label = number.of.miRNAs), position = position_fill(vjust = 0.5), size = 2)
filt_miRNA_number_bar_plot
ggsave('Figures/Bar_Graphs/miRNA_Numbers/Filtered_number_of_miRNAs_per_gene_2025.02.02.pdf', plot = filt_miRNA_number_bar_plot, create.dir = TRUE, width = 16, height = 10)

