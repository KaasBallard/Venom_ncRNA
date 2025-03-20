# Last Edited: 2024/11/1

#### Set up and Read Data ####

# Load in packages
library(compositions)
library(cowplot)
library(tidyverse)
library(scales)
# library(matrixStats)
library(ggrepel)
# library(ggpubr)
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
# Load libraries that are important for the model
# Load libraries that are important for the tree model
library(brms)
library(visreg)
library(ggplot2)
library(ape)
library(plyr)
library(ggtree)
library(caper)
library(phytools)
library(MuMIn)
library(ggtree)
library(tidytree)
library(tidybayes)
library(MuMIn)
library(broom.mixed)
library(lme4)
library(cmdstanr)
library(rstan)
# Source my functions
# source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')



# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")

# Create variable for the fused dataset.
# If this is in Dropbox the file needs to be on the machine
miRNA_mRNA_protein_data <- 'Data/Merged/miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.08.31.tsv'


# Read both in as data frames
miRNA_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T)
# # I am going to exclude the following genes because they are not well annotated
# excluded_genes = c(
#   'maker-scaffold-mi1-augustus-gene-59.13_crovir-transcript-12940',
#   'maker-scaffold-mi1-augustus-gene-59.20_crovir-transcript-12947',
#   'maker-scaffold-mi2-augustus-gene-22.17_crovir-transcript-739',
#   'maker-scaffold-un11-augustus-gene-5.19',
#   'XP_016876419',
#   'XP_011528471'
# )
# # I removed the above method because it isn't capturing every instance of the problem.


# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  dplyr::rename(
    'miRNA.Cluster.Original' = 'miRNA.Cluster',
    'Genes' = 'Converted.Gene.IDs',
    'miRNA.Cluster' = 'Putative.miRNA.Name'
  ) %>%
  dplyr::select(miRNA.Cluster, everything()) %>% # Move the new miRNA.Clusters to the front
  filter(!str_detect(Genes, 'maker-scaffold|augustus|XP_')) %>% # This should get rid of all of the weirdly annotated genes
  # filter(!(Genes %in% excluded_genes)) %>% 
  # filter(str_detect(Origin, 'three_prime_utr|five_prime_utr')) %>% # Filter out CDS targeting miRNAs
  dplyr::select(-'gtf.gene', -'crovir.transcript', -'Protein.Probability', -'Top.Peptide.Probability', -'Blast.Alignment.Length') %>%  # Remove columns to save memory
  # # Filter out unobserved proteins
  # dplyr::filter(
  #   Protein.Observed == 'Yes'
  # ) %>%
  # Change the name of the venom adam genes so that they are correct
  mutate(
    Genes = ifelse(Genes == "Venom_ADAM28_1", 'Venom_SVMP12', Genes),
    Genes = ifelse(Genes == 'Venom_ADAM28_2', 'ADAM28', Genes)
  ) %>% 
  # Filter out EXO because the weren't expressed at all
  dplyr::filter(!str_detect(Genes, 'EXO'))
rm(miRNA_mRNA_protein_df)


# CLR file
clr_data <- 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_2024.11.1.csv'

# Read file
all_samples_clr_df <- read.csv(clr_data) %>% 
  select(-Scaled.Protein.Expression, -Scaled.mRNA.Expression, -Intensity, -RNA.Raw.Counts)

# Also create a data frame that has only the detected proteins
all_samples_detected_proteins_only_df <- all_samples_clr_df %>% 
  filter(Protein.Detected == 'Yes')


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
# five_prime_color <- '#0072b2'
five_prime_color <- '#1B9E77'
# cds_color <- '#d55e00'
cds_color <- '#4A70B5'



# Save model files
model8.1 <- 'Data/Models/venom_model8.1_with_5e+06_iterations_from_brms_2024.11.1.rds'
model8.2 <- 'Data/Models/venom_model8.2_with_5e+05_iterations_from_brms_2024.11.1.rds'
model8.3 <- 'Data/Models/venom_model8.3_with_250000_iterations_from_brms_2024.11.1.rds'
model8.4 <- 'Data/Models/venom_model8.4_with_1e+05_iterations_from_brms_2024.11.1.rds'
model8.5 <- 'Data/Models/venom_model8.5_with_50000_iterations_from_brms_2024.11.1.rds'
model8.6 <- 'Data/Models/venom_model8.4_with_5000_iterations_from_brms_2024.10.30.rds'
model8.7 <- 'Data/Models/venom_model8.5_with_500_iterations_from_brms_2024.10.30.rds'


# Read models in
venom_model8.1 <- readRDS(model8.1)
venom_model8.2 <- readRDS(model8.2)
venom_model8.3 <- readRDS(model8.3)
venom_model8.4 <- readRDS(model8.4)
venom_model8.5 <- readRDS(model8.5)
venom_model8.6 <- readRDS(model8.6)
venom_model8.7 <- readRDS(model8.7)



#### Data manipulation ####

# Generate predicted draws for each model
draws_5M <- add_epred_draws(venom_model8.1, newdata = all_samples_clr_df, .draws = TRUE) %>% mutate(Iterations = "5,000,000")
draws_500K <- add_epred_draws(venom_model8.2, newdata = all_samples_clr_df, .draws = TRUE) %>% mutate(Iterations = "500,000")
draws_250K <- add_epred_draws(venom_model8.3, newdata = all_samples_clr_df, .draws = TRUE) %>% mutate(Iterations = "250,000")
draws_100K <- add_epred_draws(venom_model8.4, newdata = all_samples_clr_df, .draws = TRUE) %>% mutate(Iterations = "100,000")
draws_50K <- add_epred_draws(venom_model8.5, newdata = all_samples_clr_df, .draws = TRUE) %>% mutate(Iterations = "50,000")
draws_5K <- add_epred_draws(venom_model8.6, newdata = all_samples_clr_df, .draws = TRUE) %>% mutate(Iterations = "5,000")
draws_500 <- add_epred_draws(venom_model8.7, newdata = all_samples_clr_df, .draws = TRUE) %>% mutate(Iterations = "500")


# Combine all draws
all_draws <- bind_rows(draws_5M, draws_500K, draws_250K, draws_100K, draws_50K, draws_5K, draws_500) %>% 
  mutate(Iterations = factor(Iterations, levels = c("5,000,000", "500,000", "250,000", "100,000", "50,000", "5,000", "500")))

# Filter by families
first_half_families <- all_draws %>% 
  filter(
    Venom.Family == c(
      'VEGF', 'vQC', 'ohanin', 'CRISP', 'CTL', 'LAAO'
    )
  )

second_half_families <- all_draws %>% 
  filter(
    Venom.Family == c(
      'BPP', 'SVMP', 'SVSP', 'myotoxin', 'PLA2'
    )
  )

#### Create Density Plot ####

# Create color scheme for the venom genes
venom_colors <- c(
  SVMP = SVMP_color,
  ADAM = ADAM_color,
  SVSP = SVSP_color,
  PLA2 = PLA2_color,
  miRNA = miRNA_color,
  VEGF = VEGF_color,
  ohanin = ohanin_color,
  myotoxin = myotoxin_color,
  vQC = vQC_color,
  CRISP = CRISP_color,
  CTL = CTL_color,
  EXO = EXO_color,
  LAAO = LAAO_color,
  BPP = BPP_color,
  others = other_color
)

# Plot
draws_violin_plot <- ggplot(all_draws, aes(x = Genes, y = .epred, fill = Venom.Family)) +
  geom_violin(
    scale = "width",  # Adjust width to be proportional to the data
    trim = TRUE,      # Trims the tails of the violins
    alpha = 0.6       # Set transparency
  ) +
  labs(
    title = "Distribution of Predicted Draws by Iteration Count",
    x = 'Genes',
    y = 'Predicted Draws',
    fill = 'Venom Family'
  ) +
  scale_fill_manual(values = venom_colors) +   # Apply color scheme for fill
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top',
    axis.text.x = element_text(angle = 70, hjust = 1)
  ) +
  facet_wrap(~ Iterations)
draws_cluster_plot
ggsave('Figures/Model_Comparison/Iterations_violin_plot_2024.11.1.png', plot = draws_violin_plot, height = 15, width = 20, create.dir = T)


# Create a smaller graph with half of the families
draws_violin_plot_first_half <- ggplot(first_half_families, aes(x = Genes, y = .epred, fill = Venom.Family)) +
  geom_violin(
    scale = "width",  # Adjust width to be proportional to the data
    trim = TRUE,      # Trims the tails of the violins
    alpha = 0.6       # Set transparency
  ) +
  labs(
    title = "Distribution of Predicted Draws by Iteration Count",
    x = 'Genes',
    y = 'Predicted Draws',
    fill = 'Venom Family'
  ) +
  scale_fill_manual(values = venom_colors) +   # Apply color scheme for fill
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top',
    axis.text.x = element_text(angle = 70, hjust = 1)
  ) +
  facet_wrap(~ Iterations)
draws_violin_plot_first_half
ggsave('Figures/Model_Comparison/Iterations_violin_plot_first_half_of_families_2024.11.1.png', plot = draws_violin_plot_first_half, height = 15, width = 20, create.dir = T)

# Create a plot with the other half of families
draws_violin_plot_second_half <- ggplot(second_half_families, aes(x = Genes, y = .epred, fill = Venom.Family)) +
  geom_violin(
    scale = "width",  # Adjust width to be proportional to the data
    trim = TRUE,      # Trims the tails of the violins
    alpha = 0.6       # Set transparency
  ) +
  labs(
    title = "Distribution of Predicted Draws by Iteration Count",
    x = 'Genes',
    y = 'Predicted Draws',
    fill = 'Venom Family'
  ) +
  scale_fill_manual(values = venom_colors) +   # Apply color scheme for fill
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top',
    axis.text.x = element_text(angle = 70, hjust = 1)
  ) +
  facet_wrap(~ Iterations)
draws_violin_plot_second_half
ggsave('Figures/Model_Comparison/Iterations_violin_plot_second_half_of_families_2024.11.1.png', plot = draws_violin_plot_second_half, height = 15, width = 20, create.dir = T)
