# Last Edited: 2024/11/8

#### Set up and Read Data ####

# Load in packages
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
library(emmeans)
# Source my functions
# source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')

# Also load httpgd since I'm running this on server
library(httpgd)


# Set working directory
setwd('~/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/')
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
# setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")

# Create variable for the fused dataset.
# If this is in Dropbox the file needs to be on the machine
miRNA_mRNA_protein_data <- 'Data/Merged/miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.08.31.tsv'

# CLR file
clr_data <- 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_2024.11.13.csv'

# miRNA, mRNA, protein file
mirna_mrna_protein_file <- 'Data/miRNA/miRNA_mRNA_protein_data_for_regression_analysis_2024.11.13.csv'

# Read both in as data frames
miRNAaf_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T)
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
mi_df <- miRNAaf_mRNA_protein_df %>% 
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
rm(miRNAaf_mRNA_protein_df)

# Tree file for the 12 snakes samples
tree_file <- 'Data/Tree/11Snake_Tree_ViridisPolytomy.phy'

# Read tree file in
snake_tree <- ape::read.tree(tree_file)


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



#### Format the phylogeny and get the matrix ####

# View tip labels
snake_tree$tip.label

# List all the samples not in the analysis
missing_samples <- c(
  "'CV1084_viridis_Mid_M'", "'CV1082_viridis_South_M'", "'CV1095_viridis_North_M'",
  "'CV1089_viridis_South_M'", "'CV1090_cerberus_Other_M'"
)

# Drop tips that are not in the study
snake_tree <- drop.tip(snake_tree, missing_samples)
snake_tree$tip.label

# After dropping tips, check for NaN values
any(is.nan(snake_tree$edge.length) | is.na(snake_tree$edge.length))

# Recompute branch lengths after pruning
snake_tree <- ape::compute.brlen(snake_tree)

# Clean the tip labels by removing single quotes
snake_tree$tip.label <- gsub("'", "", snake_tree$tip.label) 

# Check the cleaned tip labels
snake_tree$tip.label

# Create a vector that contains names to rename the old ones
tip_names <- c(
  "CV0857_viridis_North_M" = 'LVG.4.CV0857.viridis.North.M', 
  "CV1081_viridis_Mid_M" = 'LVG.2.CV1081.viridis.Mid.M', 
  "CV1087_viridis_North_F" = 'RVG.5S.CV1087.viridis.North.F', 
  "CV1086_viridis_South_M" = 'LVG.9.CV1086.viridis.South.M', 
  "CV0985_concolor_Other_F" = 'RVG.7S.CV0985.concolor.Other.F', 
  "CV0987_lutosus_Other_F" = 'RVG.6S.CV0987.lutosus.Other.M'
)


# Rename tips correctly using a match
for (i in seq_along(tip_names)) {
  old_name <- names(tip_names)[i]
  new_name <- tip_names[i]
  snake_tree$tip.label[snake_tree$tip.label == old_name] <- new_name
}

# Check tree
snake_tree$tip.label
tree_plot <- ggtree(snake_tree) + geom_tiplab()
tree_plot

# Get tip labels from the phylogenetic tree
tree_tips <- snake_tree$tip.label

# Generate the phylogenetic covariance matrix from the tree
phylo_cov_matrix <- ape::vcv.phylo(snake_tree)




#### Read in the CLR and miRNA data ####

# Read file
all_samples_clr_df <- read.csv(clr_data)

# Also create a data frame that has only the detected proteins
all_samples_detected_proteins_only_df <- all_samples_clr_df %>% 
  filter(Protein.Detected == 'Yes')

# Read the file in for later regressions
all_samples_plot_miRNA_RNA_protein_df <- read.csv(mirna_mrna_protein_file)

# Let's remove the intensity column for later fusion
miRNA_RNA_protein_plot_df <- all_samples_plot_miRNA_RNA_protein_df %>% 
  select(-Intensity)




#### Linear Venom modeling ####


# Create different models for venom and then compare
venom_model1 <- lme4::lmer(
  CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (1 | Sample.ID),
  data = all_samples_detected_proteins_only_df
)
# boundary (singular) fit: see help('isSingular')
venom_model2 <- lme4::lmer(
  CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (1 | Genes),
  data = all_samples_detected_proteins_only_df
)
venom_model3 <- lme4::lmer(
  CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (1 | Genes) + (1 | Sample.ID),
  data = all_samples_detected_proteins_only_df
)
venom_model4 <- lme4::lmer(
  CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (CLR.Scaled.mRNA.Expression | Genes) + (CLR.Scaled.mRNA.Expression | Sample.ID),
  data = all_samples_detected_proteins_only_df
)
venom_model5 <- lm(
  CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + Sample.ID,
  data = all_samples_detected_proteins_only_df
)
venom_model6 <- lm(
  CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + Genes,
  data = all_samples_detected_proteins_only_df
)
# boundary (singular) fit: see help('isSingular')
venom_model7 <- lm(
  CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family,
  data = all_samples_detected_proteins_only_df
)


# Load the saved models
venom_model8 <- readRDS("Data/Models/mRNA_vs_Protein/Gaussian/venom_model8_with_500000_iterations_from_brms_2024.11.13.rds")

# Load models that includes miRNAs in the modeling
venom_model9 <- readRDS('Data/Models/ProtVSmRNA_with_miRNA/Gaussian/venom_model9_with_25000_iterations_from_brms_2024.11.13.rds')
venom_model10 <- readRDS('Data/Models/ProtVSmRNA_with_miRNA/Gaussian/venom_model10_with_25000_iterations_from_brms_2024.11.13.rds')
venom_model11 <- readRDS('Data/Models/ProtVSmRNA_with_miRNA/Gaussian/venom_model11_with_25000_iterations_from_brms_2024.11.13.rds')
venom_model12 <- readRDS('Data/Models/ProtVSmRNA_with_miRNA/Gaussian/venom_model12_with_25000_iterations_from_brms_2024.11.13.rds')


# Compare models
anova(venom_model1, venom_model2, venom_model3, venom_model4, venom_model5, venom_model6, venom_model7)
# > anova(venom_model1, venom_model2, venom_model3, venom_model4, venom_model5, venom_model6, venom_model7)
# refitting model(s) with ML (instead of REML)
# Data: all_samples_detected_proteins_only_df
# Models:
# ..7: CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family
# ..1: CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (1 | Sample.ID)
# ..2: CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (1 | Genes)
# ..3: CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (1 | Genes) + (1 | Sample.ID)
# ..5: CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + Sample.ID
# ..4: CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (CLR.Scaled.mRNA.Expression | Genes) + (CLR.Scaled.mRNA.Expression | Sample.ID)
# ..6: CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + Genes
#    npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# ..7   23 670.16 740.89 -312.08   624.16
# ..1   24 672.16 745.96 -312.08   624.16  0.000  1          1
# ..2   24 627.62 701.42 -289.81   579.62 44.542  0
# ..3   25 629.62 706.50 -289.81   579.62  0.000  1          1
# ..5   28 674.16 760.26 -309.08   618.16  0.000  3          1
# ..4   29 625.95 715.13 -283.97   567.95 50.210  1  1.381e-12 ***
# ..6   40 564.78 687.79 -242.39   484.78 83.166 11  3.592e-13 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

AIC(venom_model1, venom_model2, venom_model3, venom_model4, venom_model5, venom_model6, venom_model7)
# > AIC(venom_model1, venom_model2, venom_model3, venom_model4, venom_model5, venom_model6, venom_model7)
# df      AIC
# venom_model1 24 599.7457
# venom_model2 24 542.7037
# venom_model3 25 544.7037
# venom_model4 29 536.8653
# venom_model5 28 674.1581
# venom_model6 40 564.7814
# venom_model7 23 670.1605

# Run waic to see model complexity
waic(venom_model8)


# Run waic for the venom_model9
waic(venom_model9, venom_model10, venom_model11, venom_model12)

# # Add loo
# add_loo(venom_model8)
# add_loo(venom_model9)
# add_loo(venom_model10)
# add_loo(venom_model11)
# add_loo(venom_model12)


# # Run loo
# loo(venom_model9, venom_model10, venom_model11, venom_model12, moment_match = T)

# Model summaries
venom_model8_summary <- summary(venom_model8)
venom_model8_summary

venom_model9_summary <- summary(venom_model9)
venom_model9_summary

venom_model10_summary <- summary(venom_model10)
venom_model10_summary

venom_model11_summary <- summary(venom_model11)
venom_model11_summary

venom_model12_summary <- summary(venom_model12)
venom_model12_summary






#### mRNA vs protein expression for all samples ####

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


## Extract important information form the brm modle using tidybayes (venom_model8)

tidy_8 <- tidy(venom_model8)


# Create a data frame to contain epred (expected values of the response variable that don't include observation noise) draws
tidy_clr_epred_df <- all_samples_clr_df %>%
  select(
    -RNA.Raw.Counts, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model8, .draw = TRUE)

# Create a new data frame that contains the means of the epreds
tidy_clr_mean_epred_df <- tidy_clr_epred_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = CLR.Scaled.Protein.Expression - mean.epred)

# Create a data frame to hold the R^2 values
tidy_r_squared_df <- all_samples_detected_proteins_only_df %>%
  select(
    -RNA.Raw.Counts, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model8, .draw = TRUE) %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = CLR.Scaled.Protein.Expression - mean.epred) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    R2 = 1 - sum((CLR.Scaled.Protein.Expression - mean.epred)^2) / sum((CLR.Scaled.Protein.Expression - mean(CLR.Scaled.Protein.Expression))^2)
  )
tidy_r_squared_df


# Create a data frame to hold the R^2 values
tidy_r_squared_df2 <- all_samples_detected_proteins_only_df %>%
  select(
    -RNA.Raw.Counts, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model8, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = CLR.Scaled.Protein.Expression - .epred
  ) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
tidy_r_squared_df2





# Create figure for mRNA and Protein expression for the simple brms model (venom_model8)
all_gene_expression_plot_brm_epred <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +
  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8
  geom_line(
    data = tidy_clr_epred_df,
    aes(y = .epred, group = .draw),
    linewidth = 1,
    linetype = 'solid',
    show.legend = F,
    alpha = 1/4
  ) +  # Add epreds lines

  # Add R² text
  geom_text_repel(
    data = tidy_r_squared_df2,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
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
  )
all_gene_expression_plot_brm_epred
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_mRNAvsProtein_Expression_Plot_venom_model8_brm_epred_2024.11.13.png', plot = all_gene_expression_plot_brm_epred, height = 15, width = 20, create.dir = T)

# Add facet wrap
all_gene_expression_plot_brm_epred2  <- all_gene_expression_plot_brm_epred + facet_wrap(~ Venom.Family)
all_gene_expression_plot_brm_epred2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model8_brm_epred_2024.11.13.png', plot = all_gene_expression_plot_brm_epred2, height = 15, width = 20, create.dir = T)


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
all_gene_expression_epred_plot <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = tidy_clr_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +

  # Add R² text
  geom_text_repel(
    data = tidy_r_squared_df2,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
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
  )
all_gene_expression_epred_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_mRNAvsProtein_Expression_Plot_With_BRM_Epred_Method_venom_model8_2024.11.13.png', plot = all_gene_expression_epred_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
all_gene_expression_epred_plot2 <- all_gene_expression_epred_plot + facet_wrap(~ Venom.Family)
all_gene_expression_epred_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_mRNAvsProtein_Expression_Plot_With_BRM_Epred_Method_Separated_by_Venom_Family_venom_model8_2024.11.13.png', plot = all_gene_expression_epred_plot2, height = 15, width = 20, create.dir = T)



# Create a mean epred plot
all_gene_expression_mean_epred_plot <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add mean predictions
  geom_line(data = tidy_clr_mean_epred_df, aes(y = mean.epred), linewidth = 1, linetype = 'solid', show.legend = F) +  # Add mean predicted lines

  # Add R² text
  geom_text_repel(
    data = tidy_r_squared_df2,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
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
  )
all_gene_expression_mean_epred_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_mRNAvsProtein_Expression_Plot_With_BRM_Mean_Epred_Method_venom_model8_2024.11.13.pdf', plot = all_gene_expression_mean_epred_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
all_gene_expression_mean_epred_plot2 <- all_gene_expression_mean_epred_plot + facet_wrap(~ Venom.Family)
all_gene_expression_mean_epred_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_mRNAvsProtein_Expression_Plot_With_BRM_Mean_Epred_Method_Separated_by_Venom_Family_venom_model8_2024.11.13.png', plot = all_gene_expression_mean_epred_plot2, height = 15, width = 20, create.dir = T)










#### miRNA vs Protein for all samples ####

# Since zero values are not treated correctly by the CLR, I need to find the lowest non-zero expression level
min_nonzero_mirna <- all_samples_plot_miRNA_RNA_protein_df %>%
  filter(miRNA.RPM > 0) %>%
  summarise(min_value = min(miRNA.RPM)) %>%
  pull(min_value)
print(min_nonzero_mirna)

# Initialize miRNA pseudo count
miRNA_pseudo_count <- 1e-7
protein_pseduo_count <- 1e-7

# Calculate CLR
all_samples_mirna_clr_df <- all_samples_plot_miRNA_RNA_protein_df %>%
  group_by(Sample.ID) %>%
  mutate(Scaled.Protein.Expression = Intensity / sum(Intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.miRNA.Expression = miRNA.RPM / sum(miRNA.RPM)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  # Replace zeros with the small constant before CLR
  mutate(Scaled.Protein.Expression = ifelse(Scaled.Protein.Expression == 0, protein_pseduo_count, Scaled.Protein.Expression)) %>%
  mutate(Scaled.miRNA.Expression = ifelse(Scaled.miRNA.Expression == 0, miRNA_pseudo_count, Scaled.miRNA.Expression)) %>%
  # Perform CLR transformation
  mutate(CLR.Scaled.Protein.Expression = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(CLR.Scaled.miRNA.Expression = as.numeric(compositions::clr(Scaled.miRNA.Expression)))

# Remove any undetected proteins from the data for fitting purposes
all_samples_detected_proteins_only_mirna_clr_df <- all_samples_mirna_clr_df %>%
  filter(Protein.Detected == 'Yes')



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

# Create figure for mRNA and Protein expression
venom_family_specific_expression_plot <- ggplot(
  all_samples_detected_proteins_only_mirna_clr_df,
  aes(x = CLR.Scaled.miRNA.Expression, y = CLR.Scaled.Protein.Expression)
) +
  # Add points used for fitting
  geom_point(
    data = subset(all_samples_mirna_clr_df, Protein.Detected == 'Yes'),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +
  # Add points not used for fitting
  geom_point(
    data = subset(all_samples_mirna_clr_df, Protein.Detected == 'No'),
    aes(color = Venom.Family, fill = Venom.Family),
    size = 2.5, shape = 21, show.legend = F
  ) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(miRNA expression)',
    color = 'Venom Family'
  ) +
  scale_fill_manual(values = venom_colors) +
  scale_color_manual(values = venom_colors) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
venom_family_specific_expression_plot
# As expected, there is no correlation if you do all miRNAs against the protein expression for the genes they target
# Save
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/mRNAvsProteinModeling/All_Individuals_microRNAvsProtein_Expression_Plot_2024.11.13.pdf', plot = venom_family_specific_expression_plot, create.dir = T)





#### Data frame for residuals and miRNA (rpm) ####

# Remove the following characters from the Sample.ID. column
characters_to_remove <- c(
  'LVG.2.',
  'LVG.9.',
  'LVG.4.',
  'RVG.5S.',
  'RVG.6S.',
  'RVG.7S.'
)
characters_to_remove2 <- c(
  '.viridis.Mid.M',
  '.viridis.North.M',
  '.viridis.South.M',
  '.viridis.North.F',
  '.concolor.Other.F',
  '.lutosus.Other.M'
)


# Create name intersection for data
names <- intersect(names(tidy_clr_mean_epred_df), names(all_samples_clr_df))
names2 <- intersect(names(miRNA_RNA_protein_plot_df %>% select(-RNA.Raw.Counts, -contains('Total.'))), names(all_samples_clr_df))

# Combine the data frames
all_samples_residuals_miRNA_df <- left_join(
  all_samples_clr_df,
  tidy_clr_mean_epred_df,
  by = names
) %>%
  left_join(
    miRNA_RNA_protein_plot_df %>% select(-RNA.Raw.Counts, -contains('Total.')),
    by = names2
  ) %>%
  # Shrink size to save RAM
  select(
    -contains('Scaled')
  ) %>%
  distinct() %>%
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove, collapse = "|"), '')) %>%
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove2, collapse = "|"), ''))
glimpse(all_samples_residuals_miRNA_df)

wider_all_samples_residuals_miRNA_df <- all_samples_residuals_miRNA_df %>%
  pivot_wider(
    names_from = 'miRNA.Cluster',
    values_from = 'miRNA.RPM'
  )
glimpse(wider_all_samples_residuals_miRNA_df)



# Let's if residuals can be generally correlated with miRNA expression levels
residuals_mirna_general_plot <- ggplot(
  all_samples_residuals_miRNA_df,
  aes(x = miRNA.VST, y = residuals, fill = Venom.Family, color = Venom.Family)
) +
  # Add points
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(
    method = 'lm',
    se = T,
    linetype = 'dashed',
    formula = y ~ x,
    show.legend = F
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'Protein residuals vs miRNA expression',
    y = 'Residuals of protein',
    x = 'miRNA VST',
    fill = 'Venom Family'
  ) +
  scale_fill_manual(values = venom_colors) +
  scale_color_manual(values = venom_colors) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
residuals_mirna_general_plot

# Create a second graph with a facet
residuals_mirna_general_plot2 <- residuals_mirna_general_plot + facet_wrap( ~ Venom.Family)
residuals_mirna_general_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/General/Residual_protein_expression_vs_miRNA_expression_2024.11.13.pdf', plot = residuals_mirna_general_plot2, create.dir = T)



#### miRNA vs residual protein expression Graphs (based on VST) ####

# Get the number of available cores
num_cores <- detectCores()


source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/Residuals_Correlation_Function.R')

# Create a list of venom genes for a function to iterate through
venom_genes <- wider_all_samples_residuals_miRNA_df %>% distinct(Genes)
venom_genes <- venom_genes$Genes
# Create path for the plots to go in
path <- 'Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Residuals_Plots'

# # Loop the function
# for (gene in venom_genes) {
#  # Run the function in a loop
#  plots <- residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path, '2024.11.13', filter_r_squared = F)
# }


# Add only the plots with R squared higher than
# Create a second path
path2 <- 'Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Residuals_Plots/Good_R_squared'

# # Loop through the function
# for (gene in venom_genes) {
#  # Run the function in a loop
#  plots <- residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path2, '2024.11.13', filter_r_squared = T)
# }



# # Define a helper function to call residuals_vs_miRNA_plot2 with different arguments
# run_plot <- function(gene, filter_r_squared, path) {
#   residuals_vs_miRNA_plot2(
#     data = wider_all_samples_residuals_miRNA_df,
#     Gene = gene,
#     parent_directory = path,
#     Date = '2024.11.13',
#     filter_r_squared = filter_r_squared
#   )
# }

# # Run the function in parallel for the first case (no R-squared filtering)
# mclapply(venom_genes, run_plot, filter_r_squared = FALSE, path = path, mc.cores = num_cores)

# Define your batch size
batch_size <- 5  # Adjust this depending on your system's memory limits

# Split genes into batches
gene_batches <- split(venom_genes, ceiling(seq_along(venom_genes)/batch_size))

# Parallelize the execution of batches
mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path, '2024.11.13', filter_r_squared = FALSE)
    gc()
  }
}, mc.cores = num_cores)

# # Run the function in parallel for the second case (filtered by R-squared)
# mclapply(venom_genes, run_plot, filter_r_squared = TRUE, path = path2, mc.cores = num_cores)

mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path2, '2024.11.13', filter_r_squared = TRUE)
    gc()
  }
}, mc.cores = num_cores)





# Create a dataframe just for SVSP7 and filter at any missing data
vegf1_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_VEGF1') %>%
  dplyr::select(where(~any(!is.na(.)))) %>%
  mutate(Color = case_when(
    grepl('SVMP', Genes) ~ SVMP_color,
    grepl('ADAM', Genes) ~ ADAM_color,
    grepl('VEGF', Genes) ~ VEGF_color,
    grepl('ohanin', Genes) ~ ohanin_color,
    grepl('vQC', Genes) ~ vQC_color,
    grepl('SVSP', Genes) ~ SVSP_color,
    grepl('PLA2', Genes) ~ PLA2_color,
    grepl('CRISP', Genes) ~ CRISP_color,
    grepl('CTL', Genes) ~ CTL_color,
    grepl('EXO', Genes) ~ EXO_color,
    grepl('LAAO', Genes) ~ LAAO_color,
    grepl('myotoxin', Genes) ~ myotoxin_color,
    grepl('BPP', Genes) ~ BPP_color,
    TRUE ~ other_color
  ))

# Create a color vector
colors <- setNames(vegf1_df$Color, vegf1_df$Venom.Family)

# Create graph of values
miR_181c_5p_vegf1_plot <- ggplot(data = vegf1_df, aes(x = `cvi-miR-181c-5p`, y = residuals)) +
  geom_point(aes(color = Venom.Family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miR-181c-5p miRNA rpm vs residual protein expression',
    y = 'Residual Protein Expression',
    x = 'miRNA RPM',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
miR_181c_5p_vegf1_plot



# Create a dataframe just for SVSP7 and filter at any missing data
svmp6_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_SVMP6') %>%
  dplyr::select(where(~any(!is.na(.)))) %>%
  mutate(Color = case_when(
    grepl('SVMP', Genes) ~ SVMP_color,
    grepl('ADAM', Genes) ~ ADAM_color,
    grepl('VEGF', Genes) ~ VEGF_color,
    grepl('ohanin', Genes) ~ ohanin_color,
    grepl('vQC', Genes) ~ vQC_color,
    grepl('SVSP', Genes) ~ SVSP_color,
    grepl('PLA2', Genes) ~ PLA2_color,
    grepl('CRISP', Genes) ~ CRISP_color,
    grepl('CTL', Genes) ~ CTL_color,
    grepl('EXO', Genes) ~ EXO_color,
    grepl('LAAO', Genes) ~ LAAO_color,
    grepl('myotoxin', Genes) ~ myotoxin_color,
    grepl('BPP', Genes) ~ BPP_color,
    TRUE ~ other_color
  ))

# Create a color vector
colors <- setNames(svmp6_df$Color, svmp6_df$Venom.Family)


# Create graph of values
miR_32_3p_svmp6_plot <- ggplot(data = svmp6_df,aes(x = `cvi-miR-32-3p`, y = residuals)) +
  geom_point(aes(color = Venom.Family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'miR-32-3p miRNA RPM vs SVMP6 protein expression',
    y = 'Residual Protein Expression',
    x = 'miRNA RPM',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
miR_32_3p_svmp6_plot


ctl4_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_CTL4') %>%
  dplyr::select(where(~any(!is.na(.)))) %>%
  mutate(Color = case_when(
    grepl('SVMP', Genes) ~ SVMP_color,
    grepl('ADAM', Genes) ~ ADAM_color,
    grepl('VEGF', Genes) ~ VEGF_color,
    grepl('ohanin', Genes) ~ ohanin_color,
    grepl('vQC', Genes) ~ vQC_color,
    grepl('SVSP', Genes) ~ SVSP_color,
    grepl('PLA2', Genes) ~ PLA2_color,
    grepl('CRISP', Genes) ~ CRISP_color,
    grepl('CTL', Genes) ~ CTL_color,
    grepl('EXO', Genes) ~ EXO_color,
    grepl('LAAO', Genes) ~ LAAO_color,
    grepl('myotoxin', Genes) ~ myotoxin_color,
    grepl('BPP', Genes) ~ BPP_color,
    TRUE ~ other_color
  ))
# Create a color vector
colors <- setNames(ctl4_df$Color, ctl4_df$Venom.Family)

# Create graph of values
cluster_341_CTL4_plot <- ctl4_df %>%
  ggplot(aes(x = `Cluster_341`, y = residuals)) +
  geom_point(aes(color = Venom.Family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'Cluster 341 miRNA RPM vs residual protein expression',
    y = 'Residual Protein Expression',
    x = 'miRNA RPM',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
cluster_341_CTL4_plot


# Create a dataframe just for SVSP7 and filter at any missing data
myotoxin_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_myotoxin') %>%
  dplyr::select(where(~any(!is.na(.)))) %>%
  mutate(Color = case_when(
    grepl('SVMP', Genes) ~ SVMP_color,
    grepl('ADAM', Genes) ~ ADAM_color,
    grepl('VEGF', Genes) ~ VEGF_color,
    grepl('ohanin', Genes) ~ ohanin_color,
    grepl('vQC', Genes) ~ vQC_color,
    grepl('SVSP', Genes) ~ SVSP_color,
    grepl('PLA2', Genes) ~ PLA2_color,
    grepl('CRISP', Genes) ~ CRISP_color,
    grepl('CTL', Genes) ~ CTL_color,
    grepl('EXO', Genes) ~ EXO_color,
    grepl('LAAO', Genes) ~ LAAO_color,
    grepl('myotoxin', Genes) ~ myotoxin_color,
    grepl('BPP', Genes) ~ BPP_color,
    TRUE ~ other_color
  ))

# Create a color vector
colors <- setNames(myotoxin_df$Color, myotoxin_df$Venom.Family)


# Create graph of values
miR_34a_5p_myotoxin_plot <- myotoxin_df %>%
  ggplot(aes(x = `cvi-miR-34a-5p`, y = residuals)) +
  geom_point(aes(color = Venom.Family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'cvi-miR-34a-5p miRNA VST vs residual protein expression',
    y = 'Residual Protein Expression',
    x = 'miRNA RPM',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
miR_34a_5p_myotoxin_plot





#### Bubble Plot of R^2 values and Binding Score ####

# Create a data frame that contains the origin information
origin_df <- mi_df %>%
  dplyr::select(
    miRNA.Cluster, miRNA.Start, miRNA.End, miRNA.Strandedness, Genes, miRNA.Target.Start, miRNA.Target.End, Total.Score, Total.Energy, Max.Score, Max.Energy, Origin
  ) %>%
  distinct()

# Create a shared set of names
shared_column_names_for_origin_fusion = intersect(names(origin_df), names(all_samples_residuals_miRNA_df))


# Create a data frame that contains residauls and other important data
miRNA_residuals_target_type_df <- left_join(
  all_samples_residuals_miRNA_df,
  origin_df,
  by = shared_column_names_for_origin_fusion,
  relationship = 'many-to-many'
) %>%
  dplyr::select(
    Sample.ID, Genes, Venom.Family, residuals, miRNA.Cluster, Total.Score, Total.Energy, miRNA.RPM, Origin
  ) %>%
  mutate(
    Venom.Color = case_when(
      grepl('SVMP', Genes) ~ SVMP_color,
      grepl('VEGF', Genes) ~ VEGF_color,
      grepl('ohanin', Genes) ~ ohanin_color,
      grepl('vQC', Genes) ~ vQC_color,
      grepl('SVSP', Genes) ~ SVSP_color,
      grepl('PLA2', Genes) ~ PLA2_color,
      grepl('CRISP', Genes) ~ CRISP_color,
      grepl('CTL', Genes) ~ CTL_color,
      grepl('EXO', Genes) ~ EXO_color,
      grepl('LAAO', Genes) ~ LAAO_color,
      grepl('myotoxin', Genes) ~ myotoxin_color,
      grepl('BPP', Genes) ~ BPP_color,
      TRUE ~ other_color
      )
  ) %>%
  mutate(Genes = str_remove(Genes, '^Venom_')) %>%
  distinct()

# Define a function to calculate the R squared value for each miRNA-Gene pair
calc_r_squared_and_slope <- function(data) {
  model <- lm(residuals ~ miRNA.RPM, data = data)
  r_squared <- summary(model)$r.squared
  slope <- coef(model)['miRNA.RPM'] # Extract the slope (coefficient for miRNA.RPM)

  return(list(r_squared = r_squared, slope = slope))
}

# Calculate the R^2 values grouped by Genes and miRNA.Cluster
r_squared_df <- miRNA_residuals_target_type_df %>%
  group_by(Genes, miRNA.Cluster) %>%
  summarize(result = list(calc_r_squared_and_slope(pick(everything()))), .groups = 'drop') %>%
  unnest_wider(result, names_sep = "_") %>%  # Unpack the list into separate columns
  rename(
    'R.Squared' = 'result_r_squared',
    'Slope' = 'result_slope'
  )

# Merge the R^2 values back into the original df
r_squared_df2 <- left_join(
  miRNA_residuals_target_type_df,
  r_squared_df,
  by = c('Genes', 'miRNA.Cluster'),
) %>%
  # dplyr::filter(
  #   miRNA.Cluster == c('cvi-miR-148a-3p', 'cvi-miR-737-5p'),
  #   Genes == 'Venom_CRISP3'
  # ) %>%  # This just checks that the R squared values are the same for the same miRNA-Gene relationship
  dplyr::select(
    -Sample.ID, -residuals, -miRNA.RPM
  ) %>% # This just removes unneccessary columns that could create repeats
  distinct() %>%
  mutate(
    Origin.Color = case_when(
      grepl('three_prime_utr', Origin) ~ three_prime_color,
      grepl('five_prime_utr', Origin) ~ five_prime_color,
      grepl('CDS', Origin) ~ cds_color,
      TRUE ~ other_color
    )
  )
# View the updated data frame
print(head(r_squared_df2, n = 10))



# First lets try and see if binding score and binding energy correlate with the R^2 values at all
# Binding energy
be_vs_r_squared_plot <- ggplot(r_squared_df2, aes(x = abs(Total.Energy), y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_linedraw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Energy',
    y = 'R Squared',
    title = 'Binding Energy vs R Squared'
  )
be_vs_r_squared_plot
# They do not correlate at all
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/miRNA_Binding_Energy_vs_RSquared_Regressions_2024.11.13.pdf", plot = be_vs_r_squared_plot, width = 10, height = 10, dpi = 900, create.dir = T)

# Create the same as above but include lables
be_vs_r_squared_plot2 <- ggplot(r_squared_df2, aes(x = abs(Total.Energy), y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  geom_text(
    aes(label = paste(Genes, miRNA.Cluster, sep = ": ")),  # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_linedraw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Energy',
    y = 'R Squared',
    title = 'Binding Energy vs R Squared'
  )
be_vs_r_squared_plot2
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/Labeled/Labeled_miRNA_Binding_Energy_vs_RSquared_Regressions_2024.11.13.pdf", plot = be_vs_r_squared_plot2, width = 25, height = 25, dpi = 900, limitsize = F, create.dir = T)



# Binding Score
bs_vs_r_squared_plot <- ggplot(r_squared_df2, aes(x = Total.Score, y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_linedraw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Score',
    y = 'R Squared',
    title = 'Binding Score vs R Squared'
  )
bs_vs_r_squared_plot
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/miRNA_Binding_Score_vs_RSquared_Regressions_2024.11.13.pdf", plot = bs_vs_r_squared_plot, width = 10, height = 10, dpi = 900, create.dir = T)

# Same as above, but labeled to hell
bs_vs_r_squared_plot2 <- ggplot(r_squared_df2, aes(x = Total.Score, y = R.Squared)) +
  geom_point(aes(color = Venom.Family, shape = Origin)) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
  ) +
  geom_text(
    aes(label = paste(Genes, miRNA.Cluster, sep = ": ")),  # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  ) +
  scale_color_manual(values = setNames(r_squared_df2$Venom.Color, r_squared_df2$Venom.Family)) +  # Use custom colors
  scale_y_continuous(breaks = pretty(r_squared_df2$R.Squared, n = 20)) +  # Increase y-axis labels
  theme_linedraw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Score',
    y = 'R Squared',
    title = 'Binding Score vs R Squared'
  )
bs_vs_r_squared_plot2
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/Labeled/Labeled_miRNA_Binding_Score_vs_RSquared_Regressions_2024.11.13.pdf", plot = bs_vs_r_squared_plot2, width = 25, height = 25, dpi = 900, create.dir = T)


# Rescale the R squared values so that negatively sloped miRNA gene relationships are negative
r_squared_df3 <- r_squared_df2 %>%
  mutate(
    R.Squared = ifelse(Slope < 0, -R.Squared, R.Squared) # This adjust the sign of the R.Squared based on slope
  )


# Create bubble plot that encodes R^2 as the color and Binding Score as the size
r_squared_binding_score_bubble_plot <- ggplot(r_squared_df3, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = Total.Score, fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  # geom_point(aes(size = Total.Score, fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +  # Use custom colors
  scale_size_continuous(range = c(1, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'R Squared',
    title = 'R Squared and Binding Score'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_score_bubble_plot
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/RSquared_and_Binding_Score_2024.11.13.pdf", plot = r_squared_binding_score_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = T)


# Create bubble plot that encodes R^2 as the color and Binding Score as the size
r_squared_binding_energy_bubble_plot <- ggplot(r_squared_df3, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = abs(Total.Energy), fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  # geom_point(aes(size = abs(Total.Energy), fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +  # Use custom colors
  scale_size_continuous(range = c(1, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'R Squared',
    title = 'R Squared and Binding Energy'
  ) +
  # theme_linedraw2() +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_energy_bubble_plot
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/RSquared_and_Binding_Energy_2024.11.13.pdf", plot = r_squared_binding_energy_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = T)


# Set limits for grey and red color scale
min_r_squared = 0.5
max_total_energy = -7
min_total_score = 155

# Lets try the above again, but this time filtering by R^2 >= 0.5, and Binding Energy <= -7
filtered_r_squared_binding_energy_bubble_plot <- r_squared_df3 %>%
  dplyr::filter(
    abs(R.Squared) >= min_r_squared,
    Total.Energy <= max_total_energy
  ) %>%
  ggplot(aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = abs(Total.Energy), fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +  # 'stroke' controls the width of the outline
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_gradient(low = 'grey', high = 'red', limits = c(min_r_squared, 1)) +
  # scale_fill_viridis_c(option = 'magma') +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +  # Use custom colors
  scale_size_continuous(range = c(5, 15)) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'R Squared',
    title = 'b. R Squared and Binding Energy'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
filtered_r_squared_binding_energy_bubble_plot
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Filtered_RSquared_and_Binding_Energy_2024.11.13.pdf", plot = filtered_r_squared_binding_energy_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)


# Let's do the same but filter out low binding energies instead
# Check the range of sizes in the filtered data
filtered_r_squared_df <- r_squared_df3 %>%
  dplyr::filter(
    abs(R.Squared) >= min_r_squared,
    Total.Score >= min_total_score
  )

range_size <- range(filtered_r_squared_df$Total.Score, na.rm = TRUE)

filtered_r_squared_binding_score_bubble_plot <- ggplot(filtered_r_squared_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = Total.Score, fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +
  # geom_point(aes(size = Total.Score, fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +
  # scale_fill_viridis_c(option = 'magma') +
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_fill_gradient(low = 'grey', high = 'red', limits = c(min_r_squared, 1)) +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +
  scale_size_continuous(
    range = c(1, 15),  # Adjust these values as needed
    breaks = seq(from = min(range_size), to = max(range_size), length.out = 5)  # Adjust breaks to include the size range
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'R Squared',
    title = 'R Squared and Binding Score'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
filtered_r_squared_binding_score_bubble_plot
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Filtered_RSquared_and_Binding_Score_2024.11.13.pdf", plot = filtered_r_squared_binding_score_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)




# Lets do something else, lets use a "contact" heat map to show well Gene-miRNA pairs with strong evidence of an effect on protein expression
# Set binding score and energy
max_total_energy2 <- -20
min_total_score2 <- 160
min_r_squared2 = 0.5


# Create plot
R_squared_heatmap <- r_squared_df3 %>%
  dplyr::filter(
    Total.Score >= min_total_score2,
    Total.Energy <= max_total_energy2,
    abs(R.Squared) >= min_r_squared2
  ) %>% # Filter by binding energy and score
  ggplot(aes(x = Genes, y = miRNA.Cluster)) +
  geom_tile(aes(fill = R.Squared)) +
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    fill = 'R Squared',
    title = 'R Squared Heatmap'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
R_squared_heatmap
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Heatmaps/Filtered_RSquared_Heatmap_(BS=160_BE=-20)_2024.11.13.pdf", plot = R_squared_heatmap, width = 10, height = 15, dpi = 900, create.dir = T)


# Let's try maximizing the scores for this and redoing the figure
max_scores_df <- r_squared_df3 %>%
  mutate(
    Total.Energy.Scaled = scales::rescale(-Total.Energy), # Negate if lower energy is better
    Total.Score.Scaled = scales::rescale(Total.Score)     # Normalize Total.Score between 0 and 1
  ) %>%
  # Combine the scaled scores into a ranking metric (here we use the sum)
  mutate(
    Total.Support = abs(R.Squared) + Total.Energy.Scaled + Total.Score.Scaled
  ) %>%
  # Sort by the combined score and select the top 30
  arrange(desc(Total.Support)) %>%
  slice_head(n = 30)
# View the top 10 genes
max_scores_df

# Remove unneccessary columns and save the data
max_scores_df2 <- max_scores_df %>%
  dplyr::select(
    -contains('Color'), -Total.Energy.Scaled, -Total.Score.Scaled
  ) %>%
  dplyr::rename(
    'Target.Type' = 'Origin'
  ) %>%
  mutate(
    R.Squared = abs(R.Squared)
  )
write_csv(max_scores_df2, file = 'Tables/Expression_Tables/Venom_Family_Specific_BRM/Top30_Best_miRNA-Gene_Pairs_2024.11.13.csv')

# Create plot
R_squared_heatmap_best_30 <- ggplot(max_scores_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_tile(aes(fill = R.Squared)) +
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    fill = 'R Squared',
    title = 'R Squared Heatmap'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key=element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
R_squared_heatmap_best_30
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Heatmaps/Filtered_RSquared_Heatmap_Best_30_2024.11.13.pdf", plot = R_squared_heatmap_best_30, width = 10, height = 15, dpi = 900, create.dir = T)



# Let's try bubble plots with the same data
# Binding score
r_squared_binding_score_bubble_plot_top_30 <- ggplot(max_scores_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = Total.Score, fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +
  # geom_point(aes(size = Total.Score, fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +
  # scale_fill_viridis_c(option = 'magma') +
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +
  scale_size_continuous(
    range = c(1, 15)  # Adjust these values as needed
    ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Score',
    fill = 'R Squared',
    title = 'R Squared and Binding Score'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_score_bubble_plot_top_30
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Top30_RSquared_and_Binding_Score_2024.11.13.pdf", plot = r_squared_binding_score_bubble_plot_top_30, width = 10, height = 15, dpi = 900, create.dir = T)


# Binding energy
r_squared_binding_energy_bubble_plot_top_30 <- ggplot(max_scores_df, aes(x = Genes, y = miRNA.Cluster)) +
  geom_point(aes(size = abs(Total.Energy), fill = R.Squared), alpha = 0.75, shape = 21, stroke = 1) +
  # geom_point(aes(size = abs(Total.Energy), fill = R.Squared, color = Origin), alpha = 0.75, shape = 21, stroke = 1) +
  # scale_fill_viridis_c(option = 'magma') +
  scale_fill_gradient2(
    low = 'red',
    mid = 'white',
    high = 'blue',
    midpoint = 0,
    breaks = seq(-1, 1, by = 0.2), # Increase number of breaks
    labels = abs(seq(-1, 1, by = 0.2))  # Show absolute values in the legend
  ) +
  # scale_color_manual(name = 'Binding Target', values = setNames(r_squared_df2$Origin.Color, r_squared_df2$Origin)) +
  scale_size_continuous(
    range = c(1, 15)  # Adjust these values as needed
  ) +
  labs(
    x = 'Genes',
    y = 'miRNAs',
    size = 'Binding Energy',
    fill = 'R Squared',
    title = 'R Squared and Binding Energy'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(colour = 'black', face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 10, angle = 70, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8),
    legend.text = element_text(size = 10, face ="bold", colour ="black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  facet_grid(
    Origin ~ ., scales = 'free_y', space = 'free',
    labeller = labeller(
      Origin = c(
        'three_prime_utr' = "3' UTR Targeting",
        'five_prime_utr' = "5' UTR Targeting",
        'CDS' = "Coding Sequence Targeting"
      )
    )
  )
r_squared_binding_energy_bubble_plot_top_30
ggsave("Figures/Expression_Plots/Venom_Family_Specific_BRM/Residuals/Bubble_Plot/Top30_RSquared_and_Binding_Energy_2024.11.13.pdf", plot = r_squared_binding_energy_bubble_plot_top_30, width = 10, height = 15, dpi = 900, create.dir = T)






#### Modeling with miRNAs included ####

## Format data for use

# Read in miRNA-CLR data
miRNA_clr_df <- read.csv(clr_data)
glimpse(miRNA_clr_df)

# Remove rows with undetected proteins
miRNA_clr_detected_proteins_only_df <- miRNA_clr_df %>%
  filter(Protein.Detected == 'Yes')


# See the tidy data
tidy_9 <- tidy(venom_model9)

# Create a data frame to contain epreds
tidy_9_epred_df <- miRNA_clr_df %>%
  # Get epreds
  tidybayes::add_epred_draws(venom_model9, .draw = T)
glimpse(tidy_9_epred_df)

# Create a new data frame that contains the means of the epreds
tidy_9_mean_epred_df <- tidy_9_epred_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = CLR.Scaled.Protein.Expression - mean.epred)


# Create a data frame to hold the R^2 values
tidy_9_r_squared_df <- miRNA_clr_detected_proteins_only_df %>%
  select(
    -RNA.Raw.Counts, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model9, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = CLR.Scaled.Protein.Expression - .epred
  ) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
tidy_9_r_squared_df



# Create figure for mRNA and Protein expression for the simple brms model (venom_model9) using the stat_distribution method and epreds
venom_model9_plot <- ggplot(data = miRNA_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = tidy_9_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +

  # Add R² text
  geom_text_repel(
    data = tidy_9_r_squared_df,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
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
  )
venom_model9_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model9_plot_2024.11.13.png', plot = venom_model9_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
venom_model9_plot2 <- venom_model9_plot + facet_wrap(~ Venom.Family)
venom_model9_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model9_plot_individual_families_2024.11.13.png', plot = venom_model9_plot2, height = 15, width = 20, create.dir = T)





# See the tidy data
tidy_10 <- tidy(venom_model10)

# Create a data frame to contain epreds
tidy_10_epred_df <- miRNA_clr_df %>%
  # Get epreds
  tidybayes::add_epred_draws(venom_model10, .draw = T)
glimpse(tidy_10_epred_df)

# Create a new data frame that contains the means of the epreds
tidy_10_mean_epred_df <- tidy_10_epred_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = CLR.Scaled.Protein.Expression - mean.epred)


# Create a data frame to hold the R^2 values
tidy_10_r_squared_df <- miRNA_clr_detected_proteins_only_df %>%
  select(
    -RNA.Raw.Counts, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model10, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = CLR.Scaled.Protein.Expression - .epred
  ) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
tidy_10_r_squared_df



# Create figure for mRNA and Protein expression for the simple brms model (venom_model10) using the stat_distribution method and epreds
venom_model10_plot <- ggplot(data = miRNA_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = tidy_10_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +

  # Add R² text
  geom_text_repel(
    data = tidy_10_r_squared_df,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
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
  )
venom_model10_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model10_plot_2024.11.13.png', plot = venom_model10_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
venom_model10_plot2 <- venom_model10_plot + facet_wrap(~ Venom.Family)
venom_model10_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model10_plot_individual_families_2024.11.13.png', plot = venom_model10_plot2, height = 15, width = 20, create.dir = T)




# See the tidy data
tidy_11 <- tidy(venom_model11)

# Create a data frame to contain epreds
tidy_11_epred_df <- miRNA_clr_df %>%
  # Get epreds
  tidybayes::add_epred_draws(venom_model11, .draw = T)
glimpse(tidy_11_epred_df)

# Create a new data frame that contains the means of the epreds
tidy_11_mean_epred_df <- tidy_11_epred_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = CLR.Scaled.Protein.Expression - mean.epred)


# Create a data frame to hold the R^2 values
tidy_11_r_squared_df <- miRNA_clr_detected_proteins_only_df %>%
  select(
    -RNA.Raw.Counts, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model11, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = CLR.Scaled.Protein.Expression - .epred
  ) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
tidy_11_r_squared_df



# Create figure for mRNA and Protein expression for the simple brms model (venom_model11) using the stat_distribution method and epreds
venom_model11_plot <- ggplot(data = miRNA_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = tidy_11_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +

  # Add R² text
  geom_text_repel(
    data = tidy_11_r_squared_df,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
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
  )
venom_model11_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model11_plot_2024.11.13.png', plot = venom_model11_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
venom_model11_plot2 <- venom_model11_plot + facet_wrap(~ Venom.Family)
venom_model11_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model11_plot_individual_families_2024.11.13.png', plot = venom_model11_plot2, height = 15, width = 20, create.dir = T)






# See the tidy data
tidy_12 <- tidy(venom_model12)

# Create a data frame to contain epreds
tidy_12_epred_df <- miRNA_clr_df %>%
  # Get epreds
  tidybayes::add_epred_draws(venom_model12, .draw = T)
glimpse(tidy_12_epred_df)

# Create a new data frame that contains the means of the epreds
tidy_12_mean_epred_df <- tidy_12_epred_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  dplyr::select(-.row) %>%
  mutate(residuals = CLR.Scaled.Protein.Expression - mean.epred)


# Create a data frame to hold the R^2 values
tidy_12_r_squared_df <- miRNA_clr_detected_proteins_only_df %>%
  select(
    -RNA.Raw.Counts, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model12, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = CLR.Scaled.Protein.Expression - .epred
  ) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
tidy_12_r_squared_df



# Create figure for mRNA and Protein expression for the simple brms model (venom_model12) using the stat_distribution method and epreds
venom_model12_plot <- ggplot(data = miRNA_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(miRNA_clr_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = tidy_12_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +

  # Add R² text
  geom_text_repel(
    data = tidy_12_r_squared_df,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
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
  )
venom_model12_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model12_plot_2024.11.13.png', plot = venom_model12_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
venom_model12_plot2 <- venom_model12_plot + facet_wrap(~ Venom.Family)
venom_model12_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model12_plot_individual_families_2024.11.13.png', plot = venom_model12_plot2, height = 15, width = 20, create.dir = T)





#### Modeling untransformed data #### 

# Get un-transformed data
mirna_mrna_protein_df <- read.csv(mirna_mrna_protein_file)
glimpse(mirna_mrna_protein_df)

# Remove miRNA info
mRNA_protein_df <- mirna_mrna_protein_df %>% 
  select(
    -contains('miRNA'), -contains('Max')
  ) %>% 
  distinct()
glimpse(mRNA_protein_df)  

# Read venom model in
venom_model13 <- readRDS('Data/Models/mRNA_vs_Protein/Gaussian/venom_model13_with_500000_iterations_from_brms_2024.11.13.rds')

# Model summary
venom_model13_summary <- summary(venom_model13)
venom_model13_summary

# See the tidy data
tidy_13 <- tidy(venom_model13)
glimpse(tidy_13)

# Create a data frame to contain epreds
tidy_13_epred_df <- mRNA_protein_df %>%
  # Get epreds
  tidybayes::add_epred_draws(venom_model13, .draw = T)
glimpse(tidy_13_epred_df)

# Create a new data frame that contains the means of the epreds
tidy_13_mean_epred_df <- tidy_13_epred_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  dplyr::select(-.row) %>%
  mutate(residuals = Intensity - mean.epred)
glimpse(tidy_13_mean_epred_df)

# Create a data frame to hold the R^2 values
tidy_13_r_squared_df <- mRNA_protein_df %>%
  select(
    -RNA.Raw.Counts
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model13, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = Intensity - .epred
  ) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
tidy_13_r_squared_df



# Create figure for mRNA and Protein expression for the simple brms model (venom_model13) using the stat_distribution method and epreds
venom_model13_plot <- ggplot(data = mRNA_protein_df, aes(x = mRNA.VST, y = Intensity, color = Venom.Family)) +
  
  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(mRNA_protein_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +
  
  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(mRNA_protein_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black
  
  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = tidy_13_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +
  
  # Add R² text
  geom_text_repel(
    data = tidy_13_r_squared_df,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +
  
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'Mass spectroscopy intensity',
    x = 'mRNA expression',
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
  )
venom_model13_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model13_plot_2024.11.13.png', plot = venom_model13_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
venom_model13_plot2 <- venom_model13_plot + facet_wrap(~ Venom.Family)
venom_model13_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model13_plot_individual_families_2024.11.13.png', plot = venom_model13_plot2, height = 15, width = 20, create.dir = T)


## Venom model 14
# Read venom model in
venom_model14 <- readRDS('Data/Models/ProtVSmRNA_with_miRNA/Gaussian/venom_model14_with_25000_iterations_from_brms_2024.11.13.rds')

# Model summary
venom_model14_summary <- summary(venom_model14)
venom_model14_summary

# See the tidy data
tidy_14 <- tidy(venom_model14)
glimpse(tidy_14)

# Create a data frame to contain epreds
tidy_14_epred_df <- mirna_mrna_protein_df %>%
  # Get epreds
  tidybayes::add_epred_draws(venom_model14, .draw = T)
glimpse(tidy_14_epred_df)

# Create a new data frame that contains the means of the epreds
tidy_14_mean_epred_df <- tidy_14_epred_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  dplyr::select(-.row) %>%
  mutate(residuals = Intensity - mean.epred)
glimpse(tidy_14_mean_epred_df)

# Create a data frame to hold the R^2 values
tidy_14_r_squared_df <- mirna_mrna_protein_df %>%
  select(
    -RNA.Raw.Counts
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model14, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = Intensity - .epred
  ) %>%
  ungroup() %>%
  group_by(Venom.Family) %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
tidy_14_r_squared_df



# Create figure for mRNA and Protein expression for the simple brms model (venom_model14) using the stat_distribution method and epreds
venom_model14_plot <- ggplot(data = mRNA_protein_df, aes(x = mRNA.VST, y = Intensity, color = Venom.Family)) +
  
  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(mRNA_protein_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5, show.legend = F
  ) +
  
  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(mRNA_protein_df, Protein.Detected == "Yes"),
    aes(fill = Venom.Family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black
  
  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = tidy_14_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +
  
  # Add R² text
  geom_text_repel(
    data = tidy_14_r_squared_df,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = T, show.legend = F
  ) +
  
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'Mass spectroscopy intensity',
    x = 'mRNA expression',
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
  )
venom_model14_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model14_plot_2024.11.13.png', plot = venom_model14_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
venom_model14_plot2 <- venom_model14_plot + facet_wrap(~ Venom.Family)
venom_model14_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Gaussian/venom_model14_plot_individual_families_2024.11.13.png', plot = venom_model14_plot2, height = 15, width = 20, create.dir = T)








# #### Bernoulli Modeling ####

# # Read in model
# venom_model20 <- readRDS('Data/Models/ProtVSmRNA_with_miRNA/Bernoulli/venom_model20_with_500000_iterations_from_brms_2024.11.13.rds')

# # Get model summary
# venom_model20_summary  <- summary(venom_model20)
# venom_model20_summary


# # File containing data
# binomial_model_data <- 'Data/miRNA/miRNA_mRNA_protein_data_for_binomial_model_2024.11.13.csv'

# # Read in data
# binomial_model_df  <- read.csv(binomial_model_data) %>% 
#   select(
#     -Total.Score, -Total.Energy
#   )
# glimpse(binomial_model_df)
# summary(binomial_model_df$miRNA.VST)
# summary(binomial_model_df$Protein.Detected)


# # Fixed effects visualization of the posterior distribution
# plot(venom_model20, pars = c("b_miRNA.VST", "b_Venom.Family", "b_miRNA.VST:Venom.Family"))


# # Plot the conditional effects of the posterior distribution
# conditional_effects(venom_model20, effects = "miRNA.VST:Venom.Family") %>%
#   plot()

# # Check model quality
# pp_check(venom_model20, nsamples = 100)


# # Plot the posterior distribution 
# # Generate and plot conditional effects for miRNA.VST
# miRNA_effects <- conditional_effects(venom_model20, effects = "miRNA.VST", re_formula = NULL)
# plot(miRNA_effects, plot = FALSE)[[1]] + 
#   ggplot2::labs(title = "Predicted Probability of Protein Detection by miRNA Expression",
#                 x = "miRNA Expression (VST)",
#                 y = "Predicted Probability of Protein Detection") +
#   ggplot2::theme_minimal()


# # Create a new data frame for prediction across miRNA.VST range
# miRNA_range_expansion_df <- binomial_model_df %>%
#   select(-miRNA.VST, -Venom.Family) %>%  # Remove existing columns to avoid conflict
#   expand_grid(
#     miRNA.VST = seq(min(binomial_model_df$miRNA.VST, na.rm = TRUE), 
#                     max(binomial_model_df$miRNA.VST, na.rm = TRUE), 
#                     length.out = 50),
#     Venom.Family = unique(binomial_model_df$Venom.Family)
#   )
# glimpse(miRNA_range_expansion_df)

# # Add predicted probabilities to the new data frame
# bernoulli_epred_df <- miRNA_range_expansion_df %>%
#   add_epred_draws(venom_model20, ndraws = 100) %>% 
#   select(
#     -.chain, -.iteration
#   ) # Remove rows to reduce memory consumption
# glimpse(bernoulli_epred_df)


# # Create color scheme for the venom genes
# venom_colors <- c(
#   SVMP = SVMP_color,
#   ADAM = ADAM_color,
#   SVSP = SVSP_color,
#   PLA2 = PLA2_color,
#   miRNA = miRNA_color,
#   VEGF = VEGF_color,
#   ohanin = ohanin_color,
#   myotoxin = myotoxin_color,
#   vQC = vQC_color,
#   CRISP = CRISP_color,
#   CTL = CTL_color,
#   EXO = EXO_color,
#   LAAO = LAAO_color,
#   BPP = BPP_color,
#   others = other_color
# )



# # Plot with ggplot2
# venom_model20_prob_distribution_plot <- ggplot(binomial_model_df, aes(x = miRNA.VST, y = Protein.Detected)) +
#   geom_point(
#     aes(
#       x = miRNA.VST,
#       y = Protein.Detected,
#       fill = Venom.Family
#     ),
#     shape = 21,          # Shape 21 allows for filled points
#     size = 2.5,
#     show.legend = TRUE
#   ) +
#   ggdist::stat_lineribbon(
#     data = bernoulli_epred_df,
#     aes(
#       x = miRNA.VST,
#       y = .epred, 
#       fill = Venom.Family
#     ),
#     .width = c(.99, .95, .8, .5),
#     alpha = 0.1,
#     show.legend = FALSE
#   ) +
#   scale_x_continuous() +  # Set x-axis to continuous scale
#   scale_y_continuous(limits = c(0, 1)) +  # Set y-axis to continuous scale
#   labs(
#     title = 'Probability of protein detection at given miRNA expression levels',
#     y = 'Probability',
#     x = 'miRNA expression (vst)',
#     fill = 'Venom Family'   # Only "fill" legend is needed
#   ) +
#   scale_fill_manual(values = venom_colors) +  # Apply color scheme for fill only
#   theme_linedraw() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 10, hjust = 0.5),
#     legend.position = 'bottom',
#     legend.title.position = 'top'
#   )
# # venom_model20_prob_distribution_plot
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Bernoulli/venom_model20_probability_distribution_2024.11.13.png', plot = venom_model20_prob_distribution_plot, height = 15, width = 20, create.dir = T)

# # Create seperate facets for each family
# venom_model20_prob_distribution_plot2 <- venom_model20_prob_distribution_plot + facet_wrap(~ Venom.Family)
# # venom_model20_prob_distribution_plot2
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_BRM/miRNA_modeling/Bernoulli/venom_model20_probability_distribution_individual_families_2024.11.13.png', plot = venom_model20_prob_distribution_plot2, height = 15, width = 20, create.dir = T)


#### Make Figure 6 ####

# Fuse expression plots together
expresion_plots <- plot_grid(
  miR_32_3p_svmp6_plot,
  miR_181c_5p_vegf1_plot,
  miR_34a_5p_myotoxin_plot,
  cluster_341_CTL4_plot,
  ncol = 1,
  nrow = 4,
  align = 'v'
)
expresion_plots <- ggdraw() +
  draw_label("c. Venom genes and highly effective miRNAs", fontface = 'bold', size = 14, x = 0.5, y = 0.99, hjust = 0.5) +
  draw_plot(expresion_plots, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
expresion_plots

# Add mRNA vs Protein plot with the bubble plot
bubble_plot_plus_expression_plots <- plot_grid(
  filtered_r_squared_binding_energy_bubble_plot,
  expresion_plots,
  align = 'v',
  axis = 'l',
  ncol = 2,
  nrow = 1
)
bubble_plot_plus_expression_plots

figure_6 <- plot_grid(
  all_gene_expression_epred_plot2,
  bubble_plot_plus_expression_plots,
  align = 'v',
  axis = 't',
  ncol = 1,
  nrow = 2,
  rel_widths = c(1, 3),
  rel_heights = c(1, 3)
)
figure_6
ggsave("Figures/1_Main_Figures/Figure_6/Figure_6_2024.11.13.png", plot = figure_6, width = 18, height = 35, dpi = 900, create.dir = T)


