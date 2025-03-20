# Last Edited: 2024/11/27

#### Set up and Read Data ####

# Load in packages
library(compositions)
library(cowplot)
library(tidyverse)
library(scales)
# library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(readxl)
library(ggpmisc)
library(ggplot2)
library(gganimate)
library(gifski)
library(modelr)
library(scales)
# Load parallel package for the loop later on
library(parallel)

# Load libraries that are important for the tree model
library(brms)
library(visreg)
library(ggplot2)
library(ape)
library(plyr)
library(ggtree)
library(phytools)
library(MuMIn)
library(ggtree)
library(tidytree)
library(tidybayes)
library(MuMIn)
library(broom.mixed)
library(lme4)
# library(cmdstanr)
# Source my functions
source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')

# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('C:/Users/kaasb/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA')


# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/Long_miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.11.20.tsv'

# Read both in as data frames
miRNA_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T)

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  dplyr::rename(
    'miRNA.Cluster.Original' = 'miRNA.Cluster',
    'miRNA.Cluster' = 'Putative.miRNA.Name'
  ) %>%
  dplyr::select(miRNA.Cluster, everything()) %>% # Move the new miRNA.Clusters to the front
  filter(!str_detect(Genes, 'maker-scaffold|augustus|XP_')) %>%  # This should get rid of all of the weirdly annotated genes
  # filter(str_detect(Origin, 'three_prime_utr')) %>% # Filter out non-3'UTRs
  select(-'miRNA.Cluster.Original') %>% 
  # Change the name of the venom adam genes so that they are correct
  mutate(
    Genes = ifelse(Genes == "Venom_ADAM28_1", 'Venom_SVMP12', Genes),
    Genes = ifelse(Genes == 'Venom_ADAM28_2', 'ADAM28', Genes)
  ) %>% 
  mutate(
    Venom.Family = case_when(
      grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('VEGF', Genes) ~ 'VEGF',
      grepl('ohanin', Genes) ~ 'Ohanin',
      grepl('vQC', Genes) ~ 'vQC',
      grepl('SVSP', Genes) ~ 'SVSP',
      grepl('PLA2', Genes) ~ 'PLA2',
      grepl('CRISP', Genes) ~ 'CRISP',
      grepl('CTL', Genes) ~ 'CTL',
      grepl('EXO', Genes) ~ 'EXO',
      grepl('LAAO', Genes) ~ 'LAAO',
      grepl('myotoxin', Genes) ~ 'Myotoxin',
      grepl('BPP', Genes) ~ 'BPP',
      TRUE ~ 'others'
    )
  ) %>% 
  mutate(
    Species = case_when(
      grepl('.viridis.', Sample.ID) ~ 'C.viridis',
      grepl('.lutosus.', Sample.ID) ~ 'C.lutosus',
      grepl('.concolor.', Sample.ID) ~ 'C.concolor'
    )
  )
glimpse(mi_df)
rm(miRNA_mRNA_protein_df)


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


#### miRNA vs mRNA (Colored by protein expression) ####

# Get only neccessary columns
m_p_mi_df <- mi_df %>% 
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  select(
    Sample.ID, Species, miRNA.Cluster, Genes, miRNA.VST, mRNA.VST, Intensity, Venom.Family
  ) %>% 
  distinct()

# Reset rownames
rownames(m_p_mi_df) <- NULL


# Create plot
mirna_mrna_protein_color_plot <- ggplot(
  data = m_p_mi_df,
  mapping = aes(x = miRNA.VST, y = mRNA.VST)
) +
  geom_point(
    aes(
      shape = Species, 
      color = log(Intensity + 1),
      size = log(Intensity + 1)
    )
  ) +
  scale_color_gradient2(
    low = 'blue',
    mid = 'purple',
    high = 'red',
    midpoint = (min(log(m_p_mi_df$Intensity + 1)) + max(log(m_p_mi_df$Intensity + 1))) / 2 
  ) +
  scale_size_continuous(range = c(1, 3)) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miRNA vs mRNA expression',
    x = 'miRNA expression',
    y = 'mRNA expression',
    color = 'Protein expression (log)',
    size = 'Protein expression (log)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12),
    legend.position = "right"
  ) +
  facet_wrap( ~ Genes)
mirna_mrna_protein_color_plot
# Well, that doesn't really show anything
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/miRNA_mRNA_protein_analysis_2024.12.4.pdf', plot = mirna_mrna_protein_color_plot, height = 15, width = 20, create.dir = T)


#### miRNA vs mRNA (Colored by protein expression) and filtered by score and energy ####

# Get data and format it
be_bs_m_mi_p_df <- mi_df %>% 
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  select(
    Sample.ID, Species, miRNA.Cluster, Total.Score, Total.Energy, Genes, miRNA.VST, mRNA.VST, Intensity, Venom.Family
  ) %>% 
  # Filter by binding energy and binding score
  filter(
    Total.Energy < -7,
    Total.Score > 155
  ) %>% 
  distinct()

# Reset rownames
rownames(be_bs_m_mi_p_df) <- NULL


# Create plot
be_bs_m_mi_p_plot <- ggplot(
  data = be_bs_m_mi_p_df,
  mapping = aes(x = miRNA.VST, y = mRNA.VST)
) +
  geom_point(
    aes(
      shape = Species, 
      color = log(Intensity + 1),
      size = log(Intensity + 1)
    )
  ) +
  scale_color_gradient2(
    low = 'blue',
    mid = 'purple',
    high = 'red',
    midpoint = (min(log(be_bs_m_mi_p_df$Intensity + 1)) + max(log(be_bs_m_mi_p_df$Intensity + 1))) / 2 
  ) +
  scale_size_continuous(range = c(1, 3)) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miRNA vs mRNA expression',
    x = 'miRNA expression',
    y = 'mRNA expression',
    color = 'Protein expression (log)',
    size = 'Protein expression (log)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12),
    legend.position = "right"
  ) +
  facet_wrap( ~ Genes)
be_bs_m_mi_p_plot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/Filtered_miRNA_mRNA_protein_analysis_2024.12.4.pdf', plot = be_bs_m_mi_p_plot, height = 15, width = 20, create.dir = T)




#### miRNA vs mRNA colored by protein, miRNAs added together ####


# Get miRNA total data
total_mirna_df <- mi_df %>% 
  # Filter by binding energy and binding score
  filter(
    Total.Energy < -7,
    Total.Score > 155
  ) %>%
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  select(
    Sample.ID, Species, miRNA.Cluster, Genes, miRNA.Counts, mRNA.VST, Intensity, Venom.Family
  ) %>% 
  # Group by venom gene and summarize expression
  group_by(Genes, Sample.ID) %>% 
  summarise(Total.miRNA.Expression = sum(miRNA.Counts, na.rm = T)) %>% 
  distinct()

# Gene data frame
gene_df <- mi_df %>% 
  # Filter by binding energy and binding score
  filter(
    Total.Energy < -7,
    Total.Score > 155
  ) %>%
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  select(
    Sample.ID, Species, Genes, mRNA.VST, Intensity, Venom.Family
  ) %>% 
  distinct()

# Fuse data frames
total_mirna_gene_df <- left_join(
  gene_df,
  total_mirna_df,
  by = c('Genes', 'Sample.ID')
)


# Create plot that can hopefully show something
total_mirna_plot <- ggplot(
  data = total_mirna_gene_df,
  mapping = aes(x = log(Total.miRNA.Expression + 1), y = mRNA.VST)
) +
  geom_point(
    aes(
      shape = Species, 
      color = log(Intensity + 1),
      size = log(Intensity + 1)
    )
  ) +
  scale_color_gradient2(
    low = 'blue',
    mid = 'purple',
    high = 'red',
    midpoint = (min(log(total_mirna_gene_df$Intensity + 1)) + max(log(total_mirna_gene_df$Intensity + 1))) / 2 
  ) +
  scale_size_continuous(range = c(1, 3)) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miRNA vs mRNA expression',
    x = 'miRNA expression',
    y = 'mRNA expression',
    color = 'Protein expression (log)',
    size = 'Protein expression (log)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12),
    legend.position = "right"
  ) +
  facet_wrap( ~ Genes)
total_mirna_plot



#### Number of miRNAs and mRNA ####

# Create data frame with number of miRNAs per venom gene ####
# Get data and format it
mirna_number_df <- mi_df %>%
  # # Filter by binding energy and binding score
  # filter(
  #   Total.Energy < -7,
  #   Total.Score > 155
  # ) %>% 
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  distinct(
    miRNA.Cluster, Genes
  ) %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Number = n()
  ) %>% 
  distinct()

# Get expression levels for genes
gene_df <- mi_df %>% 
  # # Filter by binding energy and binding score
  # filter(
  #   Total.Energy < -7,
  #   Total.Score > 155
  # ) %>%
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  select(
    Sample.ID, Species, Genes, mRNA.VST, Intensity, Venom.Family
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>%
  summarise(
    Mean.mRNA.Exp = mean(mRNA.VST),
    Var.mRNA.Exp = var(mRNA.VST),
    Mean.Prot.Exp = mean(Intensity),
    Var.Prot.Exp = var(Intensity)
  )

# Join the two data frames
mirna_count_df <- left_join(
  gene_df,
  mirna_number_df,
  by = c('Genes')
) %>% 
mutate(
  Venom.Family = case_when(
    grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
    grepl('VEGF', Genes) ~ 'VEGF',
    grepl('ohanin', Genes) ~ 'Ohanin',
    grepl('vQC', Genes) ~ 'vQC',
    grepl('SVSP', Genes) ~ 'SVSP',
    grepl('PLA2', Genes) ~ 'PLA2',
    grepl('CRISP', Genes) ~ 'CRISP',
    grepl('CTL', Genes) ~ 'CTL',
    grepl('EXO', Genes) ~ 'EXO',
    grepl('LAAO', Genes) ~ 'LAAO',
    grepl('myotoxin', Genes) ~ 'Myotoxin',
    grepl('BPP', Genes) ~ 'BPP',
    TRUE ~ 'others'
  )
)


# Create plot of mean protein and mRNA expression against miRNA numbers
mirna_count_plot <- ggplot(
  data = mirna_count_df,
  mapping = aes(x = miRNA.Number, y = Mean.mRNA.Exp)
) +
  geom_point(
    aes(
      color = log(Mean.Prot.Exp + 1),
      size = log(Mean.Prot.Exp + 1)
    )
  ) +
  scale_color_gradient2(
    low = 'blue',
    mid = 'purple',
    high = 'red',
    midpoint = (min(log(mirna_count_df$Mean.Prot.Exp + 1)) + max(log(mirna_count_df$Mean.Prot.Exp + 1))) / 2 
  ) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'solid',
    formula = y ~ x,
    linewidth = 0.5
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_size_continuous(range = c(2, 5)) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miRNA vs mRNA expression',
    x = 'Number of miRNAs per gene',
    y = 'Mean mRNA expression',
    color = 'Mean protein expression (log)',
    size = 'Mean protein expression (log)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )
mirna_count_plot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/miRNA_Numbers_vs_mean_mRNA_expression_2024.12.4.pdf', plot = mirna_count_plot, height = 15, width = 20, create.dir = T)

# Create plot to show how miRNA numbers relate to mRNA variance and protein variance
mirna_count_plot2 <- ggplot(
  data = mirna_count_df,
  mapping = aes(x = miRNA.Number, y = Var.mRNA.Exp)
) +
  geom_point(
    aes(
      color = log(Var.Prot.Exp + 1),
      size = log(Var.Prot.Exp + 1)
    )
  ) +
  scale_color_gradient2(
    low = 'blue',
    mid = 'purple',
    high = 'red',
    midpoint = (min(log(mirna_count_df$Var.Prot.Exp + 1)) + max(log(mirna_count_df$Var.Prot.Exp + 1))) / 2 
  ) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'solid',
    formula = y ~ x,
    linewidth = 0.5
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_size_continuous(range = c(2, 5)) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miRNA vs mRNA expression',
    x = 'Number of miRNAs per gene',
    y = 'Variance in mRNA expression',
    color = 'Variance protein expression (log)',
    size = 'Variance protein expression (log)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )
mirna_count_plot2
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/miRNA_Numbers_vs_variance_in_mRNA_expression_2024.12.4.pdf', plot = mirna_count_plot2, height = 15, width = 20, create.dir = T)




#### Total miRNAs expression and mRNA ####

# Create data frame with number of miRNAs per venom gene ####
# Get data and format it
mirna_total_expression_df <- mi_df %>%
  # # Filter by binding energy and binding score
  # filter(
  #   Total.Energy < -7,
  #   Total.Score > 155
  # ) %>% 
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  distinct(
    Sample.ID, Genes, miRNA.Counts, 
  ) %>% 
  group_by(Genes) %>% 
  summarise(Total.miRNA.Expression = sum(miRNA.Counts, na.rm = T)) %>% 
  distinct()

# Get expression levels for genes
gene_df <- mi_df %>% 
  # # Filter by binding energy and binding score
  # filter(
  #   Total.Energy < -7,
  #   Total.Score > 155
  # ) %>%
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  select(
    Sample.ID, Species, Genes, mRNA.VST, Intensity, Venom.Family
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>%
  summarise(
    Mean.mRNA.Exp = mean(mRNA.VST),
    Var.mRNA.Exp = var(mRNA.VST),
    Mean.Prot.Exp = mean(Intensity),
    Var.Prot.Exp = var(Intensity)
  )

# Join the two data frames
mirna_gene_expression_df <- left_join(
  gene_df,
  mirna_total_expression_df,
  by = c('Genes')
) %>% 
  mutate(
    Venom.Family = case_when(
      grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('VEGF', Genes) ~ 'VEGF',
      grepl('ohanin', Genes) ~ 'Ohanin',
      grepl('vQC', Genes) ~ 'vQC',
      grepl('SVSP', Genes) ~ 'SVSP',
      grepl('PLA2', Genes) ~ 'PLA2',
      grepl('CRISP', Genes) ~ 'CRISP',
      grepl('CTL', Genes) ~ 'CTL',
      grepl('EXO', Genes) ~ 'EXO',
      grepl('LAAO', Genes) ~ 'LAAO',
      grepl('myotoxin', Genes) ~ 'Myotoxin',
      grepl('BPP', Genes) ~ 'BPP',
      TRUE ~ 'others'
    )
  )


# Create plot of mean mRNA and mean protein expression
mirna_total_plot <- ggplot(
  data = mirna_gene_expression_df,
  mapping = aes(x = log(Total.miRNA.Expression + 1), y = Mean.mRNA.Exp)
) +
  geom_point(
    aes(
      color = log(Mean.Prot.Exp + 1),
      size = log(Mean.Prot.Exp + 1)
    )
  ) +
  scale_color_gradient2(
    low = 'blue',
    mid = 'purple',
    high = 'red',
    midpoint = (min(log(mirna_gene_expression_df$Mean.Prot.Exp + 1)) + max(log(mirna_gene_expression_df$Mean.Prot.Exp + 1))) / 2 
  ) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'solid',
    formula = y ~ x,
    linewidth = 0.5
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_size_continuous(range = c(1, 8)) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miRNA vs mRNA expression (log)',
    x = 'Total miRNA expression per gene',
    y = 'Mean mRNA expression',
    color = 'Mean protein expression (log)',
    size = 'Mean protein expression (log)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )
mirna_total_plot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/miRNA_total_expression_vs_mean_mRNA_expression_2024.12.4.pdf', plot = mirna_total_plot, height = 15, width = 20, create.dir = T)


# Create plot of mean mRNA and mean protein expression
mirna_total_plot2 <- ggplot(
  data = mirna_gene_expression_df,
  mapping = aes(x = log(Total.miRNA.Expression + 1), y = Var.mRNA.Exp)
) +
  geom_point(
    aes(
      color = Var.Prot.Exp,
      size = Var.Prot.Exp
    )
  ) +
  scale_color_gradient2(
    low = 'blue',
    mid = 'purple',
    high = 'red',
    midpoint = (min(mirna_gene_expression_df$Var.Prot.Exp) + max(mirna_gene_expression_df$Var.Prot.Exp)) / 2 
  ) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'solid',
    formula = y ~ x,
    linewidth = 0.5
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_size_continuous(range = c(1, 8)) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miRNA vs mRNA expression (log)',
    x = 'Total miRNA expression per gene',
    y = 'Variance in mRNA expression',
    color = 'Variance in protein expression (log)',
    size = 'Variance in protein expression (log)'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, colour ="black"),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )
mirna_total_plot2
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Dotplots/miRNA_total_expression_vs_variance_in_mRNA_2024.12.4.pdf', plot = mirna_total_plot2, height = 15, width = 20, create.dir = T)




#### Are highly expressed genes less targeted by miRNAs ####

# Let's check if more highly expressed genes have significantly more miRNAs targeting them.
# Get number of miRNAs per gene
mirna_number_df <- mi_df %>%
  # # Filter by binding energy and binding score
  # filter(
  #   Total.Energy < -7,
  #   Total.Score > 155
  # ) %>% 
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  distinct(
    miRNA.Cluster, Genes
  ) %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Number = n()
  ) %>% 
  distinct()

# Get expression levels for genes
gene_df <- mi_df %>% 
  # # Filter by binding energy and binding score
  # filter(
  #   Total.Energy < -7,
  #   Total.Score > 155
  # ) %>%
  # Filter out mRNAs that weren't in the library
  filter(
    In.Library == 'Yes',
    !str_detect(Genes, 'ADAM')
  ) %>% 
  select(
    Sample.ID, Species, Genes, mRNA.VST, Intensity, Venom.Family
  ) %>% 
  distinct()

# Join the two data frames
mirna_count_df <- left_join(
  gene_df,
  mirna_number_df,
  by = c('Genes')
) %>% 
  mutate(Log.Intensity = log(Intensity + 1)) %>% 
  mutate(
    Venom.Family = case_when(
      grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('VEGF', Genes) ~ 'VEGF',
      grepl('ohanin', Genes) ~ 'Ohanin',
      grepl('vQC', Genes) ~ 'vQC',
      grepl('SVSP', Genes) ~ 'SVSP',
      grepl('PLA2', Genes) ~ 'PLA2',
      grepl('CRISP', Genes) ~ 'CRISP',
      grepl('CTL', Genes) ~ 'CTL',
      grepl('EXO', Genes) ~ 'EXO',
      grepl('LAAO', Genes) ~ 'LAAO',
      grepl('myotoxin', Genes) ~ 'Myotoxin',
      grepl('BPP', Genes) ~ 'BPP',
      TRUE ~ 'others'
    )
  )


## mRNA expression
# Get mean expression level for binning
mean_mrna_experssion <- mean(mirna_count_df$mRNA.VST)

# Create bins expression levels
# mirna_count_df$Bin <- cut(mirna_count_df$mRNA.VST, breaks = quantile(mirna_count_df$mRNA.VST, probs = c(0, 0.33, 0.67, 1)), labels = c('Low', 'Medium', 'High'), include.lowest = T)
mirna_count_df$Bin <- ifelse(mirna_count_df$mRNA.VST > mean_mrna_experssion, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
mirna_count_df$Bin <- factor(mirna_count_df$Bin, levels = c('Low', 'High'))

# Create a box plot that show whether highly expressed mRNAs have less overall miRNAs
mirna_count_mrna_comparison_boxplot <- ggplot(
  mirna_count_df,
  aes(
    x = Bin, y = miRNA.Number, fill = Bin
  )
) +
  geom_boxplot() +
  # stat_summary(
  #   fun = mean, 
  #   geom = "point", 
  #   shape = 20, 
  #   size = 3, 
  #   color = "red"
  # ) +  # Add mean as a red point
  # scale_fill_viridis_d(option = 'plasma') + 
  scale_color_manual(values = c('blue', 'red')) + 
  labs(
    y = 'Number of miRNAs per gene',
    x = 'mRNA expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw()
mirna_count_mrna_comparison_boxplot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Boxplots/Number_of_miRNAs_against_mRNA_expression_2024.12.4.pdf', plot = mirna_count_mrna_comparison_boxplot, height = 10, width = 8, create.dir = T)








# Calculate the mean and confidence intervals for each bin
mirna_mean_ci <- mirna_count_df %>%
  group_by(Bin) %>%
  summarize(
    Mean = mean(miRNA.Number, na.rm = TRUE),
    CI_Lower = Mean - qt(0.975, df = n() - 1) * sd(miRNA.Number, na.rm = TRUE) / sqrt(n()),
    CI_Upper = Mean + qt(0.975, df = n() - 1) * sd(miRNA.Number, na.rm = TRUE) / sqrt(n())
  )

# Plot mean and confidence intervals
mirna_mean_ci_plot <- ggplot(mirna_mean_ci, aes(x = Bin, y = Mean, color = Bin)) +
  geom_point(size = 4) +  # Plot the means
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +  # Add error bars for CI
  scale_color_manual(values = c('blue', 'red')) + 
  labs(
    y = 'Mean Number of miRNAs per Gene (± 95% CI)',
    x = 'mRNA Expression Level'
  ) +
  ggpubr::stat_compare_means(
    data = mirna_count_df,
    aes(x = Bin, y = miRNA.Number),
    method = 't.test'
  ) +  # Add p-value from t-test
  coord_cartesian(ylim = c(0, 35)) +
  theme_linedraw()
mirna_mean_ci_plot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Boxplots/Number_of_miRNAs_against_mRNA_expression2_2024.12.4.pdf', plot = mirna_mean_ci_plot, height = 10, width = 8, create.dir = T)



## Protein expression
# Get mean expression level for binning
mean_prot_expession <- mean(mirna_count_df$Log.Intensity)

# Create bins for expression levels
mirna_count_df$Bin <- ifelse(mirna_count_df$Log.Intensity > mean_prot_expession, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
mirna_count_df$Bin <- factor(mirna_count_df$Bin, levels = c('Low', 'High'))

# Create a box plot that show whether highly expressed genes have less overall miRNAs
mirna_count_prot_comparison_boxplot <- ggplot(
  mirna_count_df,
  aes(
    x = Bin, y = miRNA.Number, fill = Bin
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  # scale_fill_viridis_d(option = 'plasma') + 
  labs(
    y = 'Number of miRNAs per gene',
    x = 'Protein expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw()
mirna_count_prot_comparison_boxplot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Boxplots/Number_of_miRNAs_against_Protein_expression_2024.12.4.pdf', plot = mirna_count_prot_comparison_boxplot, height = 10, width = 8, create.dir = T)
# Yes, highly expressed genes are less targeted by miRNAs! Eureka, I cry tears of joy!




#### Are highly expressed miRNAs targeting lowly expressed genes? ####

# Format data 
mirna_mrna_prot_expression_df <- mi_df %>% 
  distinct(
    Sample.ID, Species, Genes, mRNA.VST, Intensity, Venom.Family, miRNA.Cluster, miRNA.VST
  ) %>% 
  mutate(Intensity = log(Intensity + 1))

# Get mean mRNA expression
mean_mrna_experssion <- mean(mirna_mrna_prot_expression_df$mRNA.VST)

# Create bins expression levels
# mirna_count_df$Bin <- cut(mirna_count_df$mRNA.VST, breaks = quantile(mirna_count_df$mRNA.VST, probs = c(0, 0.33, 0.67, 1)), labels = c('Low', 'Medium', 'High'), include.lowest = T)
mirna_mrna_prot_expression_df$Bin <- ifelse(mirna_mrna_prot_expression_df$mRNA.VST > mean_mrna_experssion, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
mirna_mrna_prot_expression_df$Bin <- factor(mirna_mrna_prot_expression_df$Bin, levels = c('Low', 'High'))


# Create box plot to show whether miRNA expression is significantly higher for miRNAs that target lowly expressed genes
mirna_expression_mrna_comparison_boxplot <- ggplot(
  mirna_mrna_prot_expression_df,
  aes(
    x = Bin, y = miRNA.VST, fill = Bin
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  # scale_fill_viridis_d(option = 'plasma') + 
  labs(
    y = 'miRNA expression',
    x = 'mRNA expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw() 
mirna_expression_mrna_comparison_boxplot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Boxplots/miRNA_expression_against_mRNA_expression_2024.12.4.pdf', plot = mirna_expression_mrna_comparison_boxplot, height = 10, width = 8, create.dir = T)





# Get mean mRNA expression
mean_prot_experssion <- mean(mirna_mrna_prot_expression_df$Intensity)

# Create bins expression levels
# mirna_count_df$Bin <- cut(mirna_count_df$mRNA.VST, breaks = quantile(mirna_count_df$mRNA.VST, probs = c(0, 0.33, 0.67, 1)), labels = c('Low', 'Medium', 'High'), include.lowest = T)
mirna_mrna_prot_expression_df$Bin <- ifelse(mirna_mrna_prot_expression_df$Intensity > mean_prot_experssion, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
mirna_mrna_prot_expression_df$Bin <- factor(mirna_mrna_prot_expression_df$Bin, levels = c('Low', 'High'))


# Create box plot to show whether miRNA expression is significantly higher for miRNAs that target lowly expressed genes
mirna_expression_protein_comparison_boxplot <- ggplot(
  mirna_mrna_prot_expression_df,
  aes(
    x = Bin, y = miRNA.VST, fill = Bin
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  # scale_fill_viridis_d(option = 'plasma') + 
  labs(
    y = 'miRNA expression',
    x = 'mRNA expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw() 
mirna_expression_protein_comparison_boxplot
ggsave('Figures/Expression_Plots/Qualitative_Analysis/Boxplots/miRNA_expression_against_protein_expression_2024.12.4.pdf', plot = mirna_expression_protein_comparison_boxplot, height = 10, width = 8, create.dir = T)




