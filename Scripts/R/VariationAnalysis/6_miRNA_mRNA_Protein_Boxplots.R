# Last Edited: 2025/01/29

# Set up and Read Data ----

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
library(arrow)
# Source my functions
source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')


# Create variables for the fused dataset
miRNA_number_data <- 'Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filtered_miRNA_number_data <- 'Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'

# Read data
miRNA_num_df <- read_parquet(file = miRNA_number_data) %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    in.library == 'Yes'
  )
glimpse(miRNA_num_df)

# Read the filtered data
filt_miRNA_num_df <- read_parquet(file = filtered_miRNA_number_data) %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    in.library == 'Yes'
  )
glimpse(filt_miRNA_num_df)

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
  'PLA2B1', 'PLA2K', 'PLA2C1', 'PLA2A1',
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



# Number of miRNAs vs mRNA and Protein expression ----

## Get mean mRNA expression ----

# Create a variable to hold the data so that I am not overwriting it
mRNA_df <- miRNA_num_df %>% 
  filter(
    str_detect(genes, 'Venom_'),
    # feature.type == 'three_prime_utr'
  )

# Get mean expression level for binning
mean_mrna_experssion <- mean(mRNA_df$mRNA.vst)

# Create bins expression levels
mRNA_df$bin <- ifelse(mRNA_df$mRNA.vst > mean_mrna_experssion, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
mRNA_df$bin <- factor(mRNA_df$bin, levels = c('Low', 'High'))


## Create plot for mRNA ----

# Create a box plot that show whether highly expressed mRNAs have less overall miRNAs
miRNA_numb_mRNA_exp_boxplot <- ggplot(
  mRNA_df,
  aes(
    x = bin, y = number.of.miRNAs, fill = bin
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  labs(
    y = 'Number of miRNAs per gene',
    x = 'mRNA expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw()
miRNA_numb_mRNA_exp_boxplot
ggsave('Figures/Expression_Plots/Boxplots/Number_of_miRNAs_against_mRNA_expression_2025.01.29.pdf', plot = miRNA_numb_mRNA_exp_boxplot, height = 10, width = 8, create.dir = T)
# Yes, highly expressed genes are less targeted by miRNAs! Eureka, I cry tears of joy!

# Check what it looks like if you use a facet for the gene families
venom_facet_miRNA_mRNA_boxplot <- miRNA_numb_mRNA_exp_boxplot + facet_wrap( ~ feature.type)
venom_facet_miRNA_mRNA_boxplot
ggsave('Figures/Expression_Plots/Boxplots/Number_of_miRNAs_against_mRNA_expression_by_type_2025.01.29.pdf', plot = venom_facet_miRNA_mRNA_boxplot, height = 10, width = 8, create.dir = T)

## Get mean protein expression ----

# Create a variable to hold the data so that I am not overwriting it
protein_df <- miRNA_num_df %>% 
  filter(
    str_detect(genes, 'Venom_')
  )

# Get mean expression level for binning
mean_protein_experssion <- mean(protein_df$intensity)

# Create bins expression levels
protein_df$bin <- ifelse(protein_df$intensity > mean_protein_experssion, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
protein_df$bin <- factor(protein_df$bin, levels = c('Low', 'High'))


## Create plot for protein ----

# Create a box plot that show whether highly expressed mRNAs have less overall miRNAs
miRNA_numb_protein_exp_boxplot <- ggplot(
  protein_df,
  aes(
    x = bin, y = number.of.miRNAs, fill = bin
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  labs(
    y = 'Number of miRNAs per gene',
    x = 'Protein expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw()
miRNA_numb_protein_exp_boxplot
ggsave('Figures/Expression_Plots/Boxplots/Number_of_miRNAs_against_protein_expression_2025.01.29.pdf', plot = miRNA_numb_protein_exp_boxplot, height = 10, width = 8, create.dir = T)

# By target type
venom_facet_miRNA_protein_boxplot <- miRNA_numb_protein_exp_boxplot + facet_wrap( ~ feature.type)
venom_facet_miRNA_protein_boxplot
ggsave('Figures/Expression_Plots/Boxplots/Number_of_miRNAs_against_protein_expression_by_type_2025.01.29.pdf', plot = venom_facet_miRNA_protein_boxplot, height = 10, width = 8, create.dir = T)



# Filtered number of miRNAs vs mRNA and Protein expression ----

## Get mean mRNA expression ----

# Create a variable to hold the data so that I am not overwriting it
mRNA_df2 <- filt_miRNA_num_df %>% 
  filter(
    str_detect(genes, 'Venom_'),
    !is.na(feature.type)
  )

# Get mean expression level for binning
mean_mrna_experssion2 <- mean(mRNA_df2$mRNA.vst)

# Create bins expression levels
mRNA_df2$bin <- ifelse(mRNA_df2$mRNA.vst > mean_mrna_experssion2, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
mRNA_df2$bin <- factor(mRNA_df2$bin, levels = c('Low', 'High'))


## Create plot for mRNA ----

# Create a box plot that show whether highly expressed mRNAs have less overall miRNAs
miRNA_numb_mRNA_exp_boxplot2 <- ggplot(
  mRNA_df2,
  aes(
    x = bin, y = number.of.miRNAs, fill = bin
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  labs(
    y = 'Number of miRNAs per gene',
    x = 'mRNA expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw()
miRNA_numb_mRNA_exp_boxplot2
ggsave('Figures/Expression_Plots/Boxplots/Filtered_number_of_miRNAs_against_mRNA_expression_2025.01.29.pdf', plot = miRNA_numb_mRNA_exp_boxplot2, height = 10, width = 8, create.dir = T)
# Yes, highly expressed genes are less targeted by miRNAs! Eureka, I cry tears of joy!

# Check what it looks like if you use a facet for the gene families
venom_facet_miRNA_mRNA_boxplot2 <- miRNA_numb_mRNA_exp_boxplot2 + facet_wrap( ~ feature.type)
venom_facet_miRNA_mRNA_boxplot2
ggsave('Figures/Expression_Plots/Boxplots/Filtered_number_of_miRNAs_against_mRNA_expression_by_type_2025.01.29.pdf', plot = venom_facet_miRNA_mRNA_boxplot2, height = 10, width = 8, create.dir = T)


## Get mean protein expression ----

# Create a variable to hold the data so that I am not overwriting it
protein_df2 <- filt_miRNA_num_df %>% 
  filter(
    str_detect(genes, 'Venom_'),
    !is.na(feature.type)
  )

# Get mean expression level for binning
mean_protein_experssion2 <- mean(protein_df2$intensity)

# Create bins expression levels
protein_df2$bin <- ifelse(protein_df2$intensity > mean_protein_experssion2, 'High', 'Low')

# Convert Bin to a factor and set levels to display Low on the left and High on the right
protein_df2$bin <- factor(protein_df2$bin, levels = c('Low', 'High'))


## Create plot for protein ----

# Create a box plot that show whether highly expressed mRNAs have less overall miRNAs
miRNA_numb_protein_exp_boxplot2 <- ggplot(
  protein_df2,
  aes(
    x = bin, y = number.of.miRNAs, fill = bin
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  labs(
    y = 'Number of miRNAs per gene',
    x = 'Protein expression level'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw()
miRNA_numb_protein_exp_boxplot2
ggsave('Figures/Expression_Plots/Boxplots/Filtered_number_of_miRNAs_against_protein_expression_2025.01.29.pdf', plot = miRNA_numb_protein_exp_boxplot2, height = 10, width = 8, create.dir = T)

# By target type
venom_facet_miRNA_protein_boxplot2 <- miRNA_numb_protein_exp_boxplot2 + facet_wrap( ~ feature.type)
venom_facet_miRNA_protein_boxplot2
ggsave('Figures/Expression_Plots/Boxplots/Filtered_number_of_miRNAs_against_protein_expression_by_type_2025.01.29.pdf', plot = venom_facet_miRNA_protein_boxplot2, height = 10, width = 8, create.dir = T)



# # Calculate the mean and confidence intervals for each bin
# mirna_mean_ci <- mirna_count_df %>%
#   group_by(Bin) %>%
#   summarize(
#     Mean = mean(miRNA.Number, na.rm = TRUE),
#     CI_Lower = Mean - qt(0.975, df = n() - 1) * sd(miRNA.Number, na.rm = TRUE) / sqrt(n()),
#     CI_Upper = Mean + qt(0.975, df = n() - 1) * sd(miRNA.Number, na.rm = TRUE) / sqrt(n())
#   )
# 
# # Plot mean and confidence intervals
# mirna_mean_ci_plot <- ggplot(mirna_mean_ci, aes(x = Bin, y = Mean, color = Bin)) +
#   geom_point(size = 4) +  # Plot the means
#   geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +  # Add error bars for CI
#   scale_color_manual(values = c('blue', 'red')) + 
#   labs(
#     y = 'Mean Number of miRNAs per Gene (± 95% CI)',
#     x = 'mRNA Expression Level'
#   ) +
#   ggpubr::stat_compare_means(
#     data = mirna_count_df,
#     aes(x = Bin, y = miRNA.Number),
#     method = 't.test'
#   ) +  # Add p-value from t-test
#   coord_cartesian(ylim = c(0, 35)) +
#   theme_linedraw()
# mirna_mean_ci_plot
# ggsave('Figures/Expression_Plots/Qualitative_Analysis/Boxplots/Number_of_miRNAs_against_mRNA_expression2_2025.01.29.pdf', plot = mirna_mean_ci_plot, height = 10, width = 8, create.dir = T)


