# Last Edited: 2025/03/21

# Set up and Read Data ----

## Load in packages ----
library(arrow)
library(compositions)
library(cowplot)
library(tidyverse)
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
library(emmeans)
# Also load httpgd since I'm running this on server
library(httpgd)

## Read data ----
# Set working directory
setwd("/home/administrator/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA")

# Expression data for venom genes and their miRNAs
mRNA_protein_data <- 'Data/Merged/mRNA_Protein_Combined_Data_2025.01.22.parquet'
file.exists(mRNA_protein_data)

# Read in the mRNA and protein data
mRNA_protein_df <- read_parquet(file = mRNA_protein_data) %>%
  dplyr::filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    # I am not going to filter out the genes that were detected but didn't have any intensity values
    # I will filter out EXOs because they have no observed protein expression at all in even a single sample
    venom.family != 'EXO'
  ) %>%
  dplyr::mutate(
    log.protein = log(intensity + 1)
  )
glimpse(mRNA_protein_df)


# Expression data for venom genes and their miRNAs
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_Venom_2025.01.22.parquet'
file.exists(miRNA_mRNA_protein_data)

# Read reference data in as a DataFrame
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>%
  dplyr::filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    # Filter out EXOs because they have no observed protein expression
    venom.family != 'EXO'
  ) %>%
  dplyr::select(
    sample.id, genes, venom.family, mRNA.counts, mRNA.vst, intensity, protein.observed, in.library, total.energy, total.score, feature.type, miRNA.cluster, miRNA.counts, miRNA.rpm, miRNA.vst
  )
glimpse(mi_df)

## Create color scheme for the venom genes ----
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



## Read in the model ----


### Load the saved models ----
# mRNA vs Protein model with venom family as an interaction term
model2 <- readRDS('Data/Models/mRNA_vs_Protein/Gaussian/Model2_with_500000_iterations_2025.03.13.rds')

# Summarize model 2
model2_summary <- summary(model2)
model2_summary
# > model2_summary
#  Family: student 
#   Links: mu = identity; sigma = identity; nu = identity 
# Formula: log.protein ~ mRNA.vst * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)) 
#    Data: mRNA_protein_df (Number of observations: 178) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.19      0.20     0.01     0.73 1.00     4437     4365

# Regression Coefficients:
#                               Estimate Est.Error l-95% CI u-95% CI Rhat
# Intercept                        12.20      3.54     5.34    18.94 1.00
# mRNA.vst                          0.66      0.17     0.32     0.99 1.00
# venom.familyCRISP                 4.65      7.03    -9.45    17.81 1.00
# venom.familyCTL                   1.26      4.01    -6.45     9.01 1.00
# venom.familyLAAO                  3.53      7.76   -11.76    18.61 1.00
# venom.familymyotoxin             -6.65      9.16   -23.06    12.51 1.00
# venom.familyohanin               -4.72      5.64   -15.74     6.43 1.00
# venom.familyPLA2                 -5.17      3.85   -12.59     2.51 1.00
# venom.familySVMP                  2.11      3.61    -4.91     9.26 1.00
# venom.familySVSP                  4.13      3.88    -3.48    11.75 1.00
# venom.familyVEGF                 -2.86      6.62   -15.76    10.04 1.00
# venom.familyvQC                   1.16      7.06   -13.22    15.04 1.00
# mRNA.vst:venom.familyCRISP       -0.18      0.50    -1.13     0.86 1.00
# mRNA.vst:venom.familyCTL         -0.13      0.24    -0.61     0.33 1.00
# mRNA.vst:venom.familyLAAO        -0.10      0.47    -1.00     0.81 1.00
# mRNA.vst:venom.familymyotoxin     0.35      0.49    -0.67     1.22 1.00
# mRNA.vst:venom.familyohanin       0.17      0.38    -0.59     0.91 1.00
# mRNA.vst:venom.familyPLA2         0.33      0.20    -0.07     0.72 1.00
# mRNA.vst:venom.familySVMP        -0.09      0.18    -0.44     0.27 1.00
# mRNA.vst:venom.familySVSP        -0.14      0.20    -0.54     0.26 1.00
# mRNA.vst:venom.familyVEGF        -0.07      0.50    -1.08     0.89 1.00
# mRNA.vst:venom.familyvQC         -0.06      0.52    -1.09     1.01 1.00
#                               Bulk_ESS Tail_ESS
# Intercept                         4421     4077
# mRNA.vst                          4412     4182
# venom.familyCRISP                 4549     4517
# venom.familyCTL                   4150     3701
# venom.familyLAAO                  4226     4089
# venom.familymyotoxin              4378     4279
# venom.familyohanin                4717     4441
# venom.familyPLA2                  4482     4367
# venom.familySVMP                  4555     4247
# venom.familySVSP                  4356     4318
# venom.familyVEGF                  4590     4375
# venom.familyvQC                   4449     4117
# mRNA.vst:venom.familyCRISP        4556     4496
# mRNA.vst:venom.familyCTL          4321     4231
# mRNA.vst:venom.familyLAAO         4246     4127
# mRNA.vst:venom.familymyotoxin     4355     4289
# mRNA.vst:venom.familyohanin       4700     4528
# mRNA.vst:venom.familyPLA2         4502     4482
# mRNA.vst:venom.familySVMP         4629     4288
# mRNA.vst:venom.familySVSP         4432     4212
# mRNA.vst:venom.familyVEGF         4620     4253
# mRNA.vst:venom.familyvQC          4392     4161

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.93      0.12     0.70     1.19 1.00     4700     4444
# nu        2.29      0.60     1.43     3.71 1.00     4666     4398

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# Check posterior distribution
pp_check(model2)

## Extract important information form the brm modle using tidybayes (model2) ----

# Extract the data
tidy_2 <- tidy(model2)

# Create a data frame to contain epred (expected values of the response variable that don't include observation noise) draws
model2_df <- mRNA_protein_df %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(
    model2,
    .draw = TRUE,
    .chain = TRUE,
    .iteration = TRUE,
    # ndraws = 50
  )
glimpse(model2_df)

# Create a new data frame that contains the means of the epreds
model2_mean_epred_df <- model2_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = log.protein - mean.epred)
glimpse(model2_mean_epred_df)

# Create a data frame to hold the R^2 values
r2_model2 <- mRNA_protein_df %>%
  dplyr::filter(
    # Filter out any of the non-observed proteins so that they add_epred_draws can be used
    protein.observed == "Yes",
    intensity > 0
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(model2, .draw = TRUE) %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = log.protein - mean.epred) %>%
  ungroup() %>%
  group_by(venom.family) %>%
  # Calculate the non-bayesian R^2
  summarise(
    R2 = 1 - sum((log.protein - mean.epred)^2) / sum((log.protein - mean(log.protein))^2)
  )
glimpse(r2_model2)

# Create a data frame to hold the R^2 values
r2_model2_bayes <- mRNA_protein_df %>%
  dplyr::filter(
    # Filter out any of the non-observed proteins so that they add_epred_draws can be used
    protein.observed == "Yes",
    intensity > 0
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(model2, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = log.protein - .epred
  ) %>%
  ungroup() %>%
  group_by(venom.family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
r2_model2_bayes


## Create Plots ----

# Create figure for mRNA and Protein expression for the simple brms model
mRNA_vs_protein_plot <- ggplot(data = mRNA_protein_df, aes(x = mRNA.vst, y = log.protein, color = venom.family)) +
  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(mRNA_protein_df, protein.observed == "No"),
    aes(fill = venom.family),
    shape = 21, size = 2.5, show.legend = FALSE
  ) +
  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(mRNA_protein_df, protein.observed == "Yes"),
    aes(fill = venom.family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line
  geom_line(
    data = model2_df,
    aes(y = .epred, group = .draw),
    linewidth = 1,
    linetype = 'solid',
    show.legend = FALSE,
    alpha = 1/4
  ) +  # Add epreds lines

  # Add R² text
  geom_text_repel(
    data = r2_model2_bayes,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = TRUE, show.legend = FALSE
  ) +
  # # Make the coordinates fixed
  # coord_fixed() +
  # # Add a 1:1 line
  # geom_abline(slope = 1, intercept = 0, linetype = 'solid', color = 'black', linewidth = 0.25) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'peak intensity (log scaled)',
    x = 'mRNA expression (VST)',
    fill = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  ) +
  facet_wrap(~ venom.family)
mRNA_vs_protein_plot
ggsave('Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Main/Venom_model2_mRNA_vs_protein_line.2025.03.20.png', plot = mRNA_vs_protein_plot, height = 15, width = 20, create.dir = TRUE)

# # Animate the plot because cool
# all_gene_expression_plot_brm_epred3 <- all_gene_expression_plot_brm_epred2 +
#   transition_states(.draw, 0, 1) +
#   shadow_mark(past = T, future = T, alpha = 1/20, color = 'gray50')
#
# # Set the number of draws
# ndraws1 <- 50
#
# # Animate figure
# all_gene_expression_plot_brm_epred4 <- animate(
#   all_gene_expression_plot_brm_epred3,
#   nframes = ndraws1,
#   fps = 2.5,
#   width = 10000,
#   height = 10000,
#   res = 600,
#   renderer = gifski_renderer("all_gene_expression_animation.gif")
# )
# # all_gene_expression_plot_brm_epred4
# anim_save('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model8_brm_epred_animated_2024.11.13.gif', animation = all_gene_expression_plot_brm_epred4, create.dir = T)


# Create figure for mRNA and Protein expression for the simple brms model (venom_model8) using the stat_distribution method and epreds
mRNA_vs_protein_plot_ribbon <- ggplot(data = mRNA_protein_df, aes(x = mRNA.vst, y = log.protein, color = venom.family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(mRNA_protein_df, protein.observed == "No"),
    aes(fill = venom.family),
    shape = 21, size = 2.5, show.legend = FALSE
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(mRNA_protein_df, protein.observed == "Yes"),
    aes(fill = venom.family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8 with proper grouping and fill for venom.family
  ggdist::stat_lineribbon(
    data = model2_df,
    aes(y = .epred, fill = venom.family),  # Color ribbons by venom.family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = FALSE
  ) +

  # Add R² text
  geom_text_repel(
    data = r2_model2_bayes,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = TRUE, show.legend = FALSE
  ) +

  # Make the coordinates fixed
  coord_fixed() +
  # Add a 1:1 line
  geom_abline(slope = 1, intercept = 0, linetype = 'solid', color = '#AAAAAA', linewidth = 0.25) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'peak intensity (log scaled)',
    x = 'mRNA expression (VST)',
    fill = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  ) +
  facet_wrap(~ venom.family)
mRNA_vs_protein_plot_ribbon
ggsave('Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Main/Venom_model2_mRNA_vs_protein_ribbon.2025.03.20.png', plot = mRNA_vs_protein_plot_ribbon, height = 15, width = 20, create.dir = TRUE)
ggsave('Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Main/Venom_model2_mRNA_vs_protein_ribbon.2025.03.20.pdf', plot = mRNA_vs_protein_plot_ribbon, height = 15, width = 20, create.dir = TRUE)

# Create a mean epred plot
mRNA_vs_protein_plot_mean <- ggplot(data = mRNA_protein_df, aes(x = mRNA.vst, y = log.protein, color = venom.family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(mRNA_protein_df, protein.observed == "No"),
    aes(fill = venom.family),
    shape = 21, size = 2.5, show.legend = FALSE
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(mRNA_protein_df, protein.observed == "Yes"),
    aes(fill = venom.family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add mean predictions
  geom_line(data = model2_mean_epred_df, aes(y = mean.epred, color = venom.family), linewidth = 1, linetype = 'solid', show.legend = FALSE) +  # Add mean predicted lines

  # Add R² text
  geom_text_repel(
    data = r2_model2_bayes,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = TRUE, show.legend = FALSE
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'peak intensity (log scaled)',
    x = 'mRNA expression (VST)',
    color = 'Venom Family',
    fill = 'Venom Family'  # Add label for the ribbon legend
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme to points/lines
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  ) +
  facet_wrap( ~ venom.family)
mRNA_vs_protein_plot_mean
ggsave('Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Main/Venom_model2_mRNA_vs_protein_mean_line.2025.03.20.png', plot = mRNA_vs_protein_plot_mean, height = 15, width = 20, create.dir = TRUE)



## Format residuals data ----

# Create name intersection for data
names <- intersect(names(model2_mean_epred_df), names(mRNA_protein_df))
names2 <- intersect(names(mi_df), names(mRNA_protein_df))

# Combine the data frames
model2_residuals_df <- left_join(
  mRNA_protein_df,
  model2_mean_epred_df,
  by = names
) %>%
  left_join(
    mi_df,
    by = names2
  ) %>%
  select(
    -feature.type, -miRNA.counts, -total.score, -total.energy
  ) %>%
  distinct()
glimpse(model2_residuals_df)


## miRNA vs residual protein expression Graphs (based on VST) ----

# Get the number of available cores
num_cores <- detectCores()

# Source plotting functions
source('Scripts/R/Functions/Residuals_Correlation_Function.R')

# Create a list of venom genes for a function to iterate through
venom_genes <- model2_residuals_df %>% distinct(genes)
venom_genes <- venom_genes$genes

# Create path for the plots to go in
path <- 'Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Residual_Plots/All'

# Add only the plots with R squared higher than
# Create a second path
path2 <- 'Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Residual_Plots/Good_R2'

# Define your batch size
batch_size <- 5  # Adjust this depending on your system's memory limits

# Split genes into batches
gene_batches <- split(venom_genes, ceiling(seq_along(venom_genes)/batch_size))

# Parallelize the execution of batches and do all plots
mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot3(model2_residuals_df, gene, path, '2025.03.20', filter_r_squared = FALSE)
    gc()
  }
}, mc.cores = num_cores)

# For good R squared values
mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot3(model2_residuals_df, gene, path2, '2025.03.20', filter_r_squared = TRUE)
    gc()
  }
}, mc.cores = num_cores)



## Bubble plots ----

### Format data for the bubble plots ----
# Calculate the covariance
covariance_df <- model2_residuals_df %>%
  ungroup() %>%
  select(sample.id, genes, venom.family, residuals, miRNA.cluster, miRNA.vst) %>%
  group_by(venom.family, genes, miRNA.cluster) %>%
  summarise(
    pearson.cov = cov(miRNA.vst, residuals, method = 'pearson'),
    kendall.cov = cov(miRNA.vst, residuals, method = 'kendall'),
    spearman.cov = cov(miRNA.vst, residuals, method = 'spearman'),
    pearson.cor = cor(miRNA.vst, residuals, method = 'pearson'),
    # Calculate the variance in miRNA.vst and residuals to filter out genes that didn't vary enough for correlation to be meaningful
    miRNA.variance = var(miRNA.vst),
    residual.variance = var(residuals)
  )
glimpse(covariance_df)

# Read reference data in as a DataFrame
mi_df2 <- read_parquet(file = miRNA_mRNA_protein_data) %>%
  dplyr::filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_')
  ) %>%
  select(
    -contains('target'), -contains('start'), -contains('end'), -contains('strand'), -contains('chrom'), -contains('max'), -contains('length'), -contains('sequence'),
    -best.miRNA.ortholog, -miRNA.cluster.original, -E.value, -bit.score, -blast.percent.identity, -miRNA.name.probability
  )
glimpse(mi_df2)

# Create a shared set of names
shared_column_names <- intersect(names(mi_df2), names(covariance_df))

# Join the data to recover the feature type and score informations
miRNA_covariance_df <- left_join(
  covariance_df,
  mi_df2,
  by = shared_column_names,
  relationship = 'many-to-many'
) %>%
  dplyr::select(
    genes, venom.family, contains('cov'), pearson.cor, miRNA.cluster, total.score, total.energy, feature.type, contains('variance')
  )  %>%
  mutate(genes = str_remove(genes, '^Venom_')) %>%
  distinct()
glimpse(miRNA_covariance_df)


### Covariance against binding energy ----
# First lets try and see if binding score and binding energy correlate with the covariance values at all
# Binding energy
be_vs_covariance_plot <- ggplot(miRNA_covariance_df, aes(x = abs(total.energy), y = pearson.cov)) +
  geom_point(aes(color = venom.family, shape = feature.type)) +
  geom_smooth(
    method = 'lm',
    se = TRUE,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_color_manual(values = venom_colors) +  # Use custom colors
  scale_y_continuous(breaks = pretty(miRNA_covariance_df$pearson.cov, n = 20)) +  # Increase y-axis labels
  theme_linedraw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Energy',
    y = 'Covariance',
    title = 'Binding Energy vs R Squared'
  )
be_vs_covariance_plot
# They do not correlate at all
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Covariance/miRNA_expression_covariance_against_protein_residuals_vs_binding_energy.2025.03.20.pdf", plot = be_vs_covariance_plot, width = 10, height = 10, dpi = 900, create.dir = TRUE)

# Same as the above, but labeled with the miRNA cluster and gene the dot represents
be_vs_covariance_plot_labeled <- be_vs_covariance_plot +
  geom_text(
    aes(label = paste(genes, miRNA.cluster, sep = ": ")), # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  )
be_vs_covariance_plot_labeled
# They do not correlate at all
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Covariance/Labeled_miRNA_expression_covariance_against_protein_residuals_vs_binding_energy.2025.03.20.pdf", plot = be_vs_covariance_plot_labeled, width = 10, height = 10, dpi = 900, create.dir = TRUE)


# Binding Score
bs_vs_covariance_plot <- ggplot(miRNA_covariance_df, aes(x = abs(total.score), y = pearson.cov)) +
  geom_point(aes(color = venom.family, shape = feature.type)) +
  geom_smooth(
    method = 'lm',
    se = TRUE,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_color_manual(values = venom_colors) +  # Use custom colors
  scale_y_continuous(breaks = pretty(miRNA_covariance_df$pearson.cov, n = 20)) +  # Increase y-axis labels
  theme_linedraw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Score',
    y = 'Covariance',
    title = 'Binding Energy vs R Squared'
  )
bs_vs_covariance_plot
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Covariance/miRNA_expression_covariance_against_protein_residuals_vs_binding_score.2025.03.20.pdf", plot = bs_vs_covariance_plot, width = 10, height = 10, dpi = 900, create.dir = TRUE)


# Same as the above, but labeled with the miRNA cluster and gene the dot represents
bs_vs_covariance_plot_labeled <- bs_vs_covariance_plot +
  geom_text(
    aes(label = paste(genes, miRNA.cluster, sep = ": ")), # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  )
bs_vs_covariance_plot_labeled
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Covariance/Labeled_miRNA_expression_covariance_against_protein_residuals_vs_binding_score.2025.03.20.pdf", plot = bs_vs_covariance_plot_labeled, width = 10, height = 10, dpi = 900, create.dir = TRUE)


### Bubble plot of correlation coefficients ----
# Create bubble plot that encodes R as the color and Binding Score as the size
r_and_binding_score_bubble_plot <- ggplot(miRNA_covariance_df, aes(x = genes, y = miRNA.cluster)) +
  geom_point(aes(size = total.score, fill = pearson.cor), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = seq(-1, 1, by = 0.2)  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$feature.type.Color, r_squared_df2$feature.type)) +  # Use custom colors
  scale_size_continuous(range = c(1, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'Pearson correlation coefficient',
    title = 'Correlation and Binding Score'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face = "bold", colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    feature.type ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      feature.type = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_and_binding_score_bubble_plot
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Correlation/miRNA_residauls_correlation_and_bs_bubble_plot_2025.03.20.pdf", plot = r_and_binding_score_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)


# Create bubble plot that encodes R as the color and Binding Score as the size
r_and_binding_energy_bubble_plot <- ggplot(miRNA_covariance_df, aes(x = genes, y = miRNA.cluster)) +
  geom_point(aes(size = total.energy, fill = pearson.cor), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = seq(-1, 1, by = 0.2)  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$feature.type.Color, r_squared_df2$feature.type)) +  # Use custom colors
  scale_size_continuous(range = c(15, 1)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'Pearson correlation coefficient',
    title = 'Correlation and Binding Energy'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face = "bold", colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    feature.type ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      feature.type = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_and_binding_energy_bubble_plot
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Correlation/miRNA_residauls_correlation_and_be_bubble_plot_2025.03.20.pdf", plot = r_and_binding_energy_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)

### Format data for the filtered bubble plots ----

# Set limits for grey and red color scale
min_correlation <- 0.5
max_total_energy <- -7
min_total_score <- 155

miRNA_covariance_df2 <- miRNA_covariance_df %>%
  filter(
    # # Filter out any of the miRNAs that didn't vary enough for correlation to be meaningful
    # residual.variance > 0.1,
    # Filter out any pearson's correlation below the minimum to look at only interesting correlations
    abs(pearson.cor) > min_correlation,
    # Filter out any low binding miRNAs
    total.energy <= max_total_energy,
    total.score >= min_total_score
  )

### Filtered bubble plots ----

# Create bubble plot that encodes R as the color and Binding Score as the size
filt_r_and_binding_score_bubble_plot <- ggplot(miRNA_covariance_df2, aes(x = genes, y = miRNA.cluster)) +
  geom_point(aes(size = total.score, fill = pearson.cor), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = seq(-1, 1, by = 0.2)  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$feature.type.Color, r_squared_df2$feature.type)) +  # Use custom colors
  scale_size_continuous(range = c(1, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'Pearson correlation coefficient',
    title = 'Correlation and Binding Score'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face = "bold", colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    feature.type ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      feature.type = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
filt_r_and_binding_score_bubble_plot
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Correlation/filtered_miRNA_residauls_correlation_and_bs_bubble_plot_2025.03.20.pdf", plot = filt_r_and_binding_score_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)


# Create bubble plot that encodes R as the color and Binding Score as the size
filt_r_and_binding_energy_bubble_plot <- ggplot(miRNA_covariance_df2, aes(x = genes, y = miRNA.cluster)) +
  geom_point(aes(size = total.energy, fill = pearson.cor), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = seq(-1, 1, by = 0.2)  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$feature.type.Color, r_squared_df2$feature.type)) +  # Use custom colors
  scale_size_continuous(range = c(15, 1)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'Pearson correlation coefficient',
    title = 'Correlation and Binding Energy'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face = "bold", colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    feature.type ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      feature.type = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
filt_r_and_binding_energy_bubble_plot
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Correlation/filtered_miRNA_residauls_correlation_and_be_bubble_plot_2025.03.20.pdf", plot = filt_r_and_binding_energy_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)

#### Create a version with only 3' UTR targeting ----

# Remove the 5' UTR and CDS targeting
miRNA_covariance_3UTR_df <- miRNA_covariance_df2 %>%
  filter(feature.type == 'three_prime_utr') %>%
  distinct()

# Create bubble plot
# Create bubble plot that encodes R as the color and Binding Score as the size
r_and_binding_energy_bubble_plot_3UTR <- ggplot(miRNA_covariance_3UTR_df, aes(x = genes, y = miRNA.cluster)) +
  geom_point(aes(size = total.energy, fill = pearson.cor), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = seq(-1, 1, by = 0.2)  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$feature.type.Color, r_squared_df2$feature.type)) +  # Use custom colors
  scale_size_continuous(range = c(15, 1)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'Pearson correlation coefficient',
    title = 'Correlation and Binding Energy'
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face = "bold", colour = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title.position = "top"
  )
r_and_binding_energy_bubble_plot_3UTR
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Model2/Correlation/2filtered_3UTR_miRNA_residauls_correlation_and_be_bubble_plot_2025.04.01.pdf", plot = r_and_binding_energy_bubble_plot_3UTR, width = 4, height = 7, dpi = 900, create.dir = TRUE)
