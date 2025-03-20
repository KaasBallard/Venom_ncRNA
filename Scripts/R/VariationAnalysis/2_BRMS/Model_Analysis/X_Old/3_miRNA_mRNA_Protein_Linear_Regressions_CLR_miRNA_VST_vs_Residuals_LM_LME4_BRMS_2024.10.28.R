# Last Edited: 2024/10/28

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
library(caper)
library(phytools)
library(MuMIn)
library(ggtree)
library(tidytree)
library(tidybayes)
library(MuMIn)
library(broom.mixed)
library(lme4)
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
EXO_color <- '#49FFFF'
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

# CLR file
clr_data <- 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_2024.10.27.csv'

# Read file
all_samples_clr_df <- read.csv(clr_data)

# Also create a data frame that has only the detected proteins
all_samples_detected_proteins_only_df <- all_samples_clr_df %>% 
  filter(Protein.Detected == 'Yes')


# miRNA, mRNA, protein file
mirna_mrna_protein_file <- 'Data/miRNA/miRNA_mRNA_protein_data_for_regression_analysis_2024.10.27.csv'

# Read the file in for later regressions
all_samples_plot_miRNA_RNA_protein_df <- read.csv(mirna_mrna_protein_file)





#### Venom models ####



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


# Load the saved model
venom_model8 <- readRDS("Data/Models/venom_model8.1_from_brms_2024.10.27.rds")

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
# df       AIC
# venom_model1 24  974.4375
# venom_model2 24  914.0945
# venom_model3 25  916.0945
# venom_model4 29  916.8174
# venom_model5 28 1098.1736
# venom_model6 40  992.2261
# venom_model7 23 1091.3678
waic(venom_model8)
# > waic(venom_model8)
# 
# Computed from 45000 by 160 log-likelihood matrix.
# 
# Estimate   SE
# elpd_waic   -331.9 15.4
# p_waic        17.5  4.3
# waic         663.9 30.8
# 
# 8 (5.0%) p_waic estimates greater than 0.4. We recommend trying loo instead. 
# Warning message:
#   
#   8 (5.0%) p_waic estimates greater than 0.4. We recommend trying loo instead. 


## Extract lme4 venom_model2 information
# Summary
model2_summary <- summary(venom_model2)
print(model2_summary)

## Cool stuff from broom.mixed
# Get summary statistics
summary_stats_lme <- broom::tidy(venom_model2)
# Get individual observations data
augment_lme <- broom::augment(venom_model2) %>%
  rename_with(
    .cols = starts_with("."), .fn = ~ gsub("^\\.", "lme.", .)
  )

# Model stats
glance_lme <- broom::glance(venom_model2)

# Get predicted values
predictions_lme <- predict(venom_model2, newdata = all_samples_detected_proteins_only_df, re.form = NULL, se.fit = T) # The re.form = NULL means random effects are included, if you don't want to include them use re.form = NA which only uses fixed effects

# Extract fitted values and standard errors
fitted_values_lme <- predictions_lme$fit
se_fit_lme <- predictions_lme$se.fit

# Calculate confidence intervals
alpha <- 0.05  # For 95% confidence intervals
t_val_lme <- qt(1 - alpha / 2, df = df.residual(venom_model2))  # t-value for CI
ci_lower_lme <- fitted_values_lme - t_val_lme * se_fit_lme
ci_upper_lme <- fitted_values_lme + t_val_lme * se_fit_lme

# Create shared names set
shared_names_lme <- intersect(names(augment_lme), names(all_samples_detected_proteins_only_df))

# Add confidene intervals and fitted values
all_samples_clr_lme_df <- all_samples_detected_proteins_only_df %>%
  mutate(
    Fitted.lme = fitted_values_lme,
    se.fit.lme = se_fit_lme,
    CI.Lower.lme = ci_lower_lme,
    CI.Upper.lme = ci_upper_lme
  ) %>% left_join(
    augment_lme,
    by = shared_names_lme
  ) %>%
  mutate(
    Residuals.lme = CLR.Scaled.Protein.Expression - lme.fitted
  )

# Set a small tolerance for comparison (e.g., 1e-6)
tolerance <- 1e-6

# Get rows where they differ
differences_df <- all_samples_clr_lme_df %>% filter(abs(Fitted.lme - lme.fitted) > tolerance)
differences_df





## Extract info from the normal lm model: venom_model7
model7_summary <- summary(venom_model7)
print(model7_summary)

# Get summary statistics
summary_stats_lm <- broom::tidy(venom_model7)

# Get individual observations data
augment_lm <- broom::augment(venom_model7) %>%
  rename_with(
    .cols = starts_with("."), .fn = ~ gsub("^\\.", "lm.", .)
  )

# Model stats
glance_lm <- broom::glance(venom_model7)

# Get predicted values
predictions_lm <- predict(venom_model7, newdata = all_samples_detected_proteins_only_df, re.form = NA, se.fit = T, interval = 'confidence') # The re.form = NULL means random effects are included, if you don't want to include them use re.form = NA which only uses fixed effects

# Extract fitted values and standard errors
fitted_values_lm <- predictions_lm$fit[, 'fit']
se_fit_lm <- predictions_lm$se.fit
ci_lower_lm <- predictions_lm$fit[, 'lwr']
ci_upper_lm <- predictions_lm$fit[, 'upr']

# Create shared names set
shared_names_lm <- intersect(names(augment_lm), names(all_samples_detected_proteins_only_df))

# Add confidene intervals and fitted values
all_samples_clr_lm_df <- all_samples_detected_proteins_only_df %>%
  mutate(
    Fitted.lm = fitted_values_lm,
    se.fit.lm = se_fit_lm,
    CI.Lower.lm = ci_lower_lm,
    CI.Upper.lm = ci_upper_lm
  ) %>% left_join(
    augment_lm,
    by = shared_names_lm
  ) %>%
  mutate(
    Residuals.lm = CLR.Scaled.Protein.Expression - lm.fitted
  )



## Extract brms venom_model8 information
# Summary
model8_summary <- summary(venom_model8)
print(model8_summary)

# Cool stuff from broom.mixed
# Get summary statistics
summary_stats_brm <- broom.mixed::tidy(venom_model8)
# Get individual observations data
augment_brm <- broom.mixed::augment(venom_model8) %>%
  rename_with(
    .cols = starts_with("."), .fn = ~ gsub("^\\.", "brm.", .)
  )

# Model stats
glance_brm <- broom.mixed::glance(venom_model8)

# Get predicted values
predictions_brm <- predict(venom_model8, newdata = all_samples_detected_proteins_only_df, re.form = NULL, se.fit = T) # The re.form = NULL means random effects are included, if you don't want to include them use re.form = NA which only uses fixed effects

# Extract fitted values and standard errors
fitted_values_brm <- predictions_brm[, 'Estimate']
ci_upper_brm <- predictions_brm[, 'Q97.5']
ci_lower_brm <- predictions_brm[, 'Q2.5']


# Create shared names set
shared_names_brm <- intersect(names(augment_brm), names(all_samples_detected_proteins_only_df))

# Add confidene intervals and fitted values
all_samples_clr_brm_df <- all_samples_detected_proteins_only_df %>% 
  mutate(
    CLR.Scaled.Protein.Expression = round(CLR.Scaled.Protein.Expression, 6),
    CLR.Scaled.mRNA.Expression = round(CLR.Scaled.mRNA.Expression, 6)
  ) %>%
  left_join(
    augment_brm %>% 
      mutate(
        CLR.Scaled.Protein.Expression = round(CLR.Scaled.Protein.Expression, 6),
        CLR.Scaled.mRNA.Expression = round(CLR.Scaled.mRNA.Expression, 6)
      ),
    by = shared_names_brm
  ) %>%
  mutate(
    Fitted.brm = fitted_values_brm,
    CI.Upper.brm = ci_upper_brm,
    CI.Lower.brm = ci_lower_brm
  ) %>% 
  mutate(
    Residuals.brm = CLR.Scaled.Protein.Expression - brm.fitted
  ) # Check if the output residuals from the augment() function match what you would get by calculating it
# They don't match :(.






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



# Create figure for mRNA and Protein expression that compares all of the potential models
all_gene_expression_plot <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +
  
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
    size = 2.5, shape = 21, color = "black", stroke = 0.5
  ) +  # Outline with black
  
  # Add regression line from venom_model2
  # geom_line(data = all_samples_clr_lme_df, aes(y = lme.fitted), linewidth = 1, show.legend = F) +  # Add predicted lines
  geom_line(data = all_samples_clr_lme_df, aes(y = lme.fitted), linewidth = 0.5, linetype = 'dashed', color = 'blue') +  # Add predicted lines
  # Add credible intervals for venom_model2
  # geom_ribbon(data = all_samples_clr_lme_df, aes(ymin = CI.Lower.lme, ymax = CI.Upper.lme), fill = 'blue', color = NA, alpha = 0.2, show.legend = F) +  # Add credible intervals

  # Add regression line from venom_model7
  # geom_line(data = all_samples_clr_lm_df, aes(y = lm.fitted), linewidth = 1, show.legend = F) +  # Add predicted lines
  geom_line(data = all_samples_clr_lm_df, aes(y = lm.fitted), linewidth = 0.5, linetype = 'dashed', color = 'green') +  # Add predicted lines
  # Add credible intervals for venom_model7
  # geom_ribbon(data = all_samples_clr_lm_df, aes(ymin = CI.Lower.lm, ymax = CI.Upper.lm), fill = 'green', color = NA, alpha = 0.2, show.legend = F) +  # Add credible intervals

  # Add regression line from venom_model8
  # geom_line(data = all_samples_clr_brm_df, aes(y = brm.fitted), linewidth = 1, show.legend = F) +  # Add predicted lines
  geom_line(data = all_samples_clr_brm_df, aes(y = brm.fitted), linewidth = 0.5, linetype = 'dashed', color = 'black') +  # Add predicted lines
  # Add credible intervals for venom_model8
  # geom_ribbon(data = all_samples_clr_brm_df, aes(ymin = CI.Lower.brm, ymax = CI.Upper.brm), fill = 'black', color = NA, alpha = 0.2, show.legend = F) +  # Add credible intervals
  # geom_ribbon(data = all_samples_clr_brm_df, aes(ymin = CI.Lower.brm, ymax = CI.Upper.brm, fill = 'Venom.Family'), alpha = 0.2, show.legend = F) +  # Add credible intervals
  # geom_ribbon(data = all_samples_clr_brm_df, aes(ymin = CI.Lower.lme, ymax = CI.Upper.lme, fill = 'Venom.Family'), alpha = 0.2, show.legend = F) +  # Add credible intervals

  # Normal lm model (venom_model7)
  # geom_smooth(
  #   aes(group = Venom.Family, color = Venom.Family, fill = Venom.Family),
  #   method = 'lm',
  #   se = T,
  #   linetype = 'dashed',
  #   show.legend = F,
  #   formula = y ~ x
  # ) +
  # Normal lm model (venom_model7)
  # stat_poly_eq(
  #   aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
  #   formula = y ~ x
  # ) +
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
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_plot
# Save
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Model_Comparison_2024.10.28.pdf', plot = all_gene_expression_plot, height = 10, width = 10, create.dir = T)

# Add facet grid
all_gene_expression_plot2 <- all_gene_expression_plot + facet_wrap(~ Venom.Family)
all_gene_expression_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_Model_Comparison_2024.10.28.pdf', plot = all_gene_expression_plot2, height = 15, width = 20, create.dir = T)







# Create figure for mRNA and Protein expression for the simple lm model (venom_model7)
all_gene_expression_plot_lm <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +
  
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
    size = 2.5, shape = 21, color = "black", stroke = 0.5
  ) +  # Outline with black

  # Add regression line from venom_model7
  # geom_line(data = all_samples_clr_lm_df, aes(y = lm.fitted), linewidth = 1, show.legend = F) +  # Add predicted lines
  geom_line(data = all_samples_clr_lm_df, aes(y = lm.fitted), linewidth = 0.5, linetype = 'dashed', show.legend = F) +  # Add predicted lines
  # Add credible intervals for venom_model7
  geom_ribbon(data = all_samples_clr_lm_df, aes(ymin = CI.Lower.lm, ymax = CI.Upper.lm, fill = Venom.Family), color = NA, alpha = 0.2, show.legend = F) +  # Add credible intervals

  # # Normal lm model (venom_model7)
  # geom_smooth(
  #   aes(group = Venom.Family, color = Venom.Family, fill = Venom.Family),
  #   method = 'lm',
  #   se = T,
  #   linetype = 'dashed',
  #   show.legend = F,
  #   formula = y ~ x
  # ) +
  # Normal lm model (venom_model7)
  stat_poly_eq(
    data = subset(all_samples_clr_df, Protein.Detected == "Yes"),
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x
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
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_plot_lm
# Save
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_venom_model7_lm_2024.10.28.pdf', plot = all_gene_expression_plot_lm, height = 10, width = 10, create.dir = T)

# Add facet grid
all_gene_expression_plot_lm2 <- all_gene_expression_plot_lm + facet_wrap(~ Venom.Family)
all_gene_expression_plot_lm2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model7_lm_2024.10.28.pdf', plot = all_gene_expression_plot_lm2, height = 15, width = 20, create.dir = T)





# Create figure for mRNA and Protein expression for the simple lme4 model (venom_model2)
all_gene_expression_plot_lme <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +
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
    size = 2.5, shape = 21, color = "black", stroke = 0.5
  ) +  # Outline with black
  
  # Add regression line from venom_model2
  # geom_line(data = all_samples_clr_lme_df, aes(y = lme.fitted), linewidth = 1, show.legend = F) +  # Add predicted lines
  geom_line(data = all_samples_clr_lme_df, aes(y = lme.fitted), linewidth = 0.5, linetype = 'dashed', show.legend = F) +  # Add predicted lines
  # Add credible intervals for venom_model2
  geom_ribbon(data = all_samples_clr_lme_df, aes(ymin = CI.Lower.lme, ymax = CI.Upper.lme, fill = Venom.Family), color = NA, alpha = 0.2, show.legend = F) +  # Add credible intervals

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_plot_lme
# Save
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_venom_model2_lme_2024.10.28.pdf', plot = all_gene_expression_plot_lme, height = 10, width = 10, create.dir = T)

# Add facet grid
all_gene_expression_plot_lme2 <- all_gene_expression_plot_lme + facet_wrap(~ Venom.Family)
all_gene_expression_plot_lme2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model2_lme_2024.10.28.pdf', plot = all_gene_expression_plot_lme2, height = 15, width = 20, create.dir = T)






# Create figure for mRNA and Protein expression for the simple brms model (venom_model8) with data from the broom.mixed package
all_gene_expression_plot_brm <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +
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
    size = 2.5, shape = 21, color = "black", stroke = 0.5
  ) +  # Outline with black
  
  # Add regression line from venom_model8
  # geom_line(data = all_samples_clr_brm_df, aes(y = brm.fitted), linewidth = 1, show.legend = F) +  # Add predicted lines
  geom_line(data = all_samples_clr_brm_df, aes(y = brm.fitted), linewidth = 0.5, linetype = 'solid', show.legend = F) +  # Add predicted lines
  # Add credible intervals for venom_model8
  geom_ribbon(data = all_samples_clr_brm_df, aes(ymin = CI.Lower.brm, ymax = CI.Upper.brm, fill = Venom.Family), color = NA, alpha = 0.2, show.legend = F) +  # Add credible intervals

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_plot_brm
# Save
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_venom_model8_brm_2024.10.28.pdf', plot = all_gene_expression_plot_brm, height = 10, width = 10, create.dir = T)

# Add facet grid
all_gene_expression_plot_brm2 <- all_gene_expression_plot_brm + facet_wrap(~ Venom.Family)
all_gene_expression_plot_brm2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model8_brm_2024.10.28.pdf', plot = all_gene_expression_plot_brm2, height = 15, width = 20, create.dir = T)




## Extract important information form the brm modle using tidybayes (venom_model8)


# Set the number of draws
ndraws1 <- 100
ndraws2 <- 400


## Extract important information form the brm modle using tidybayes (venom_model8)

# Create a data frame to contain predictions (expected values of the response variable that do include observation noise) draws. I will use this for calculating residuals.
all_samples_clr_brm_prediction_df <- all_samples_clr_df %>%
  select(
    -RNA.VST, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Predictions are expected values of the response variable that do include observation noise
  tidybayes::add_predicted_draws(venom_model8, ndraws = ndraws2)


# Create a second data frame to contain mean predictions
all_samples_clr_brm_mean_prediction_df <- all_samples_clr_brm_prediction_df %>% 
  select(
    -.chain, -.iteration
  ) %>% 
  distinct() %>% 
  summarise(
    mean.prediction = mean(.prediction) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>% 
  select(-.row)


# Create figure for mRNA and Protein expression for the simple brms model (venom_model8)
all_gene_expression_plot_brm_pred <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +
  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5
  ) +
  
  # Add points that were used in fitting with a black outline
  geom_point(data = subset(all_samples_clr_df, Protein.Detected == "Yes"), 
             aes(fill = Venom.Family), 
             size = 2.5, shape = 21, color = "black", stroke = 0.5, show.legend = F) +  # Outline with black
  
  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = all_samples_clr_brm_prediction_df,
    aes(y = .prediction, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
  ) +
  
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme to points/lines
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_plot_brm_pred
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_venom_model8_brm_prediction_2024.10.28.pdf', plot = all_gene_expression_plot_brm_pred, height = 15, width = 20, create.dir = T)

# Add facet wrap
all_gene_expression_plot_brm_pred2  <- all_gene_expression_plot_brm_pred + facet_wrap(~ Venom.Family)
all_gene_expression_plot_brm_pred2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model8_brm_prediction_2024.10.28.pdf', plot = all_gene_expression_plot_brm_pred2, height = 15, width = 20, create.dir = T)


# Create another figure for mRNA and Protein expression for the simple brms model (venom_model8) with averaged predicted values
all_gene_expression_plot_brm_mean_pred <- ggplot(data = all_samples_clr_df, aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression, color = Venom.Family)) +
  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(all_samples_clr_df, Protein.Detected == "No"),
    aes(fill = Venom.Family),
    shape = 21, size = 2.5
  ) +
  
  # Add points that were used in fitting with a black outline
  geom_point(data = subset(all_samples_clr_df, Protein.Detected == "Yes"), 
             aes(fill = Venom.Family), 
             size = 2.5, shape = 21, color = "black", stroke = 0.5, show.legend = F) +  # Outline with black
  
  # Add mean predictions
  geom_line(data = all_samples_clr_brm_mean_prediction_df, aes(y = mean.prediction), linewidth = 1, linetype = 'solid', show.legend = F) +  # Add mean predicted lines
  
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
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_plot_brm_mean_pred
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_venom_model8_brm_mean_prediction_2024.10.28.pdf', plot = all_gene_expression_plot_brm_mean_pred, height = 15, width = 20, create.dir = T)

# Add facet wrap
all_gene_expression_plot_brm_mean_pred2  <- all_gene_expression_plot_brm_mean_pred + facet_wrap(~ Venom.Family)
all_gene_expression_plot_brm_mean_pred2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model8_brm_mean_prediction_2024.10.28.pdf', plot = all_gene_expression_plot_brm_mean_pred2, height = 15, width = 20, create.dir = T)









# Create a data frame to contain epred (expected values of the response variable that don't include observation noise) draws
all_samples_clr_brm_epred_df <- all_samples_clr_df %>%
  select(
    -RNA.VST, -Intensity, -Scaled.Protein.Expression, -Scaled.mRNA.Expression
  ) %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model8, ndraws = ndraws2)


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
    size = 2.5, shape = 21, color = "black", stroke = 0.5
  ) +  # Outline with black
  
  # Add regression line from venom_model8
  geom_line(
    data = all_samples_clr_brm_epred_df,
    aes(y = .epred, group = .draw),
    linewidth = 1,
    linetype = 'solid',
    show.legend = F,
    alpha = 1/4
  ) +  # Add epreds lines

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
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_plot_brm_epred
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_venom_model8_brm_epred_2024.10.28.pdf', plot = all_gene_expression_plot_brm_epred, height = 15, width = 20, create.dir = T)

# Add facet wrap
all_gene_expression_plot_brm_epred2  <- all_gene_expression_plot_brm_epred + facet_wrap(~ Venom.Family)
all_gene_expression_plot_brm_epred2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model8_brm_epred_2024.10.28.pdf', plot = all_gene_expression_plot_brm_epred2, height = 15, width = 20, create.dir = T)

# # Animate the plot because cool
# all_gene_expression_plot_brm_epred3 <- all_gene_expression_plot_brm_epred2 +
#   transition_states(.draw, 0, 1) +
#   shadow_mark(past = T, future = T, alpha = 1/20, color = 'gray50')
# 
# all_gene_expression_plot_brm_epred4 <- animate(
#   all_gene_expression_plot_brm_epred3, 
#   nframes = ndraws1, 
#   fps = 2.5, 
#   width = 10000, 
#   height = 10000, 
#   res = 600, 
#   renderer = gifski_renderer("all_gene_expression_animation.gif")
# )
# all_gene_expression_plot_brm_epred4
# anim_save('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_Individual_Families_venom_model8_brm_epred_animated_2024.10.28.gif', animation = all_gene_expression_plot_brm_epred4, create.dir = T)


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
    size = 2.5, shape = 21, color = "black", stroke = 0.5
  ) +  # Outline with black
  
  # Add regression line from venom_model8 with proper grouping and fill for Venom.Family
  ggdist::stat_lineribbon(
    data = all_samples_clr_brm_epred_df,
    aes(y = .epred, fill = Venom.Family),  # Color ribbons by Venom.Family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = F
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
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
all_gene_expression_epred_plot
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_With_BRM_Epred_Method_Venom_Model8_2024.10.28.pdf', plot = all_gene_expression_epred_plot, height = 15, width = 20, create.dir = T)

# Add facet wrap
all_gene_expression_epred_plot2 <- all_gene_expression_epred_plot + facet_wrap(~ Venom.Family)
all_gene_expression_epred_plot2
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_With_BRM_Epred_Method_Separated_by_Venom_Family_Venom_Model8_2024.10.28.pdf', plot = all_gene_expression_epred_plot, height = 15, width = 20, create.dir = T)











#### miRNA vs Protein for all samples ####

# Since zero values are not treated correctly by the CLR, I need to find the lowest non-zero expression level
min_nonzero_mirna <- all_samples_plot_miRNA_RNA_protein_df %>%
  filter(miRNA.VST > 0) %>%
  summarise(min_value = min(miRNA.VST)) %>%
  pull(min_value)
print(min_nonzero_mirna)

# Initialize miRNA pseudo count
miRNA_pseudo_count <- 1e-6

# Calculate CLR
all_samples_mirna_clr_df <- all_samples_plot_miRNA_RNA_protein_df %>%
  group_by(Sample.ID) %>%
  mutate(Scaled.Protein.Expression = Intensity / sum(Intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.miRNA.Expression = miRNA.VST / sum(miRNA.VST)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  # Replace zeros with the small constant before CLR
  mutate(Scaled.Protein.Expression = ifelse(Scaled.Protein.Expression == 0, protein_pseduo_count, Scaled.Protein.Expression)) %>%
  mutate(Scaled.miRNA.Expression = ifelse(Scaled.miRNA.Expression == 0, miRNA_pseudo_count, Scaled.miRNA.Expression)) %>%
  # Perform CLR transformation
  mutate(CLR.Scaled.Protein.Expression = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(CLR.Scaled.miRNA.Expression = as.numeric(compositions::clr(Scaled.miRNA.Expression)))


# Set variables for protein and RNA expression levels
all_miRNA_expression <- all_samples_mirna_clr_df$CLR.Scaled.miRNA.Expression
all_protein_expression <- all_samples_mirna_clr_df$CLR.Scaled.Protein.Expression
venom_family <- all_samples_mirna_clr_df$Venom.Family
# gene_label <- all_samples_mirna_clr_df$Genes


# Create color scheme for the venom genes
venom_colors <- c(
  SVMP = '#4A70B5',
  ADAM = '#9A70B5',
  SVSP = '#F0B830',
  PLA2 = '#7570B3',
  miRNA = '#8B0AA5',
  VEGF = '#74ADD1',
  ohanin = '#3A489C',
  myotoxin = '#B2182B',
  vQC = '#80BC50',
  CRISP = '#E7298A',
  CTL = '#F67E17',
  EXO = '#49FFFF',
  LAAO = '#B35806',
  BPP = '#1B9E77',
  others = '#666666'
)


# Create figure for mRNA and Protein expression
Venom_Family_Specific_expression_plot <- all_samples_mirna_clr_df %>%
  ggplot(aes(x = CLR.Scaled.miRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
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
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(miRNA expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
Venom_Family_Specific_expression_plot
# As expected, there is no correlation if you do all miRNAs against the protein expression for the genes they target
# Save
ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Linear_Regression_VST/All_Individuals_microRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = Venom_Family_Specific_expression_plot, create.dir = T)



#### Data frame for residuals and miRNA (VST) ####

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


# Create a data frame to contain predictions (expected values of the response variable that do include observation noise) draws. I will use this for calculating residuals.
residuals_df <- all_samples_clr_brm_prediction_df %>%
  # Get residuals by subtraction rather than add_residual_draws() because there is no way to get residual draws to match the predictions
  mutate(
    Residuals = CLR.Scaled.Protein.Expression - mean.prediction
  )



# Fit linear models for each Venom Family and extract residuals
residuals_df <- all_samples_clr_df %>%
  group_by(Venom.Family) %>%
  do({
    model <- lm(CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression, data = .)
    data.frame(Sample.ID = .$Sample.ID, Genes = .$Genes, Residuals = residuals(model), Fitted.lme = fitted(model))
  })

# Create name intersection for data
names <- intersect(names(residuals_df), names(all_samples_clr_df))
names2 <- intersect(names(all_samples_plot_miRNA_RNA_protein_df), names(all_samples_clr_df))

# Combine the data frames
all_samples_residuals_miRNA_df <- left_join(
  all_samples_clr_df,
  residuals_df,
  by = names
) %>%
  mutate(
    Res = CLR.Scaled.Protein.Expression - Fitted.lme
  ) %>%
  left_join(
    all_samples_plot_miRNA_RNA_protein_df,
    by = names2
  ) %>%
  # Shrink size to save ram
  select(
    -contains('Scaled')
  ) %>%
  distinct() %>%
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove, collapse = "|"), '')) %>%
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove2, collapse = "|"), ''))

wider_all_samples_residuals_miRNA_df <- all_samples_residuals_miRNA_df %>%
  pivot_wider(
    names_from = 'miRNA.Cluster',
    values_from = 'miRNA.VST'
  )




#### miRNA vs Residuals Graphs (based on VST) ####

# Get the number of available cores
num_cores <- detectCores()


source('Scripts/R/Protein-miRNA_RNA_Analysis/Functions/Residuals_Correlation_Function.R')

# Create a list of venom genes for a function to iterate through
venom_genes <- wider_all_samples_residuals_miRNA_df %>% distinct(Genes)
venom_genes <- venom_genes$Genes
# Create path for the plots to go in
path <- 'Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Residuals_Plots'

# # Loop the function
# for (gene in venom_genes) {
#  # Run the function in a loop
#  plots <- residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path, '2024.10.28', filter_r_squared = F)
# }


# Add only the plots with R squared higher than
# Create a second path
path2 <- 'Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Residuals_Plots/Good_R_squared'

# # Loop through the function
# for (gene in venom_genes) {
#  # Run the function in a loop
#  plots <- residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path2, '2024.10.28', filter_r_squared = T)
# }



# # Define a helper function to call residuals_vs_miRNA_plot2 with different arguments
# run_plot <- function(gene, filter_r_squared, path) {
#   residuals_vs_miRNA_plot2(
#     data = wider_all_samples_residuals_miRNA_df,
#     Gene = gene,
#     parent_directory = path,
#     Date = '2024.10.28',
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
    residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path, '2024.10.28', filter_r_squared = FALSE)
    gc()
  }
}, mc.cores = num_cores)

# # Run the function in parallel for the second case (filtered by R-squared)
# mclapply(venom_genes, run_plot, filter_r_squared = TRUE, path = path2, mc.cores = num_cores)

mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot2(wider_all_samples_residuals_miRNA_df, gene, path2, '2024.10.28', filter_r_squared = TRUE)
    gc()
  }
}, mc.cores = num_cores)





# Create a dataframe just for SVSP7 and filter at any missing data
svsp7_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_SVSP7') %>%
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
colors <- setNames(svsp7_df$Color, svsp7_df$Venom.Family)


# Create graph of values
cluster_1395_svsp7_plot <- svsp7_df %>%
  ggplot(aes(x = Cluster_1395, y = Residuals)) +
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
    title = 'Cluster_1395 miRNA VST vs SVSP7 protein expression',
    y = 'Intensity (log scaled)',
    x = 'miRNA VST',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
cluster_1395_svsp7_plot



# Create a dataframe just for SVSP7 and filter at any missing data
svmp7_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_SVMP7') %>%
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
colors <- setNames(svmp7_df$Color, svmp7_df$Venom.Family)


# Create graph of values
cluster_1395_svmp7_plot <- svmp7_df %>%
  ggplot(aes(x = Cluster_1395, y = Log.Intensity)) +
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
    title = 'Cluster_1395 miRNA VST vs SVMP7 protein expression',
    y = 'Intensity (log scaled)',
    x = 'miRNA VST',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
cluster_1395_svmp7_plot


pla2b1_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_PLA2B1') %>%
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
colors <- setNames(pla2b1_df$Color, pla2b1_df$Venom.Family)
# Create graph of values
let_7c_5p_2_PLA2B1_plot <- pla2b1_df %>%
  ggplot(aes(x = `cvi-let-7c-5p-2`, y = Log.Intensity)) +
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
    title = 'cvi-let-7c-5p-2 miRNA VST vs PLA2B1 protein expression',
    y = 'Intensity (log scaled)',
    x = 'miRNA VST',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
let_7c_5p_2_PLA2B1_plot


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
  ggplot(aes(x = `cvi-miR-34a-5p`, y = Log.Intensity)) +
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
    title = 'cvi-miR-34a-5p miRNA VST vs myotoxin protein expression',
    y = 'Intensity (log scaled)',
    x = 'miRNA VST',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
miR_34a_5p_myotoxin_plot


svmp4_df <- wider_all_samples_residuals_miRNA_df %>%
  filter(Genes == 'Venom_SVMP4') %>%
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
colors <- setNames(svmp4_df$Color, svmp4_df$Venom.Family)


# Create graph of values
cluster_919_svmp4_plot <- svmp4_df %>%
  ggplot(aes(x = Cluster_919, y = Log.Intensity)) +
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
    title = 'Cluster_919 miRNA VST vs SVMP4 protein expression',
    y = 'Intensity (log scaled)',
    x = 'miRNA VST',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5)),
    legend.position = 'none',
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10, hjust = 0.5),
    # legend.position = 'bottom',
    # legend.title.position = 'top'
  )
cluster_919_svmp4_plot




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
    Sample.ID, Genes, Venom.Family, Residuals, miRNA.Cluster, Total.Score, Total.Energy, miRNA.VST, Origin
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
  model <- lm(Residuals ~ miRNA.VST, data = data)
  r_squared <- summary(model)$r.squared
  slope <- coef(model)['miRNA.VST'] # Extract the slope (coefficient for miRNA.VST)

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
    -Sample.ID, -Residuals, -miRNA.VST
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
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Energy',
    y = 'R Squared',
    title = 'Binding Energy vs R Squared'
  )
be_vs_r_squared_plot
# They do not correlate at all
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/miRNA_Binding_Energy_vs_RSquared_Regressions_2024.10.28.pdf", plot = be_vs_r_squared_plot, width = 10, height = 10, dpi = 900, create.dir = T)

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
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Energy',
    y = 'R Squared',
    title = 'Binding Energy vs R Squared'
  )
be_vs_r_squared_plot2
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/Labeled/Labeled_miRNA_Binding_Energy_vs_RSquared_Regressions_2024.10.28.pdf", plot = be_vs_r_squared_plot2, width = 25, height = 25, dpi = 900, limitsize = F, create.dir = T)



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
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Score',
    y = 'R Squared',
    title = 'Binding Score vs R Squared'
  )
bs_vs_r_squared_plot
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/miRNA_Binding_Score_vs_RSquared_Regressions_2024.10.28.pdf", plot = bs_vs_r_squared_plot, width = 10, height = 10, dpi = 900, create.dir = T)

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
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) +
  labs(
    x = 'Binding Score',
    y = 'R Squared',
    title = 'Binding Score vs R Squared'
  )
bs_vs_r_squared_plot2
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Binding_Data_vs_RSquared/Labeled/Labeled_miRNA_Binding_Score_vs_RSquared_Regressions_2024.10.28.pdf", plot = bs_vs_r_squared_plot2, width = 25, height = 25, dpi = 900, create.dir = T)


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
  # theme_classic2() +
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
r_squared_binding_score_bubble_plot
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/RSquared_and_Binding_Score_2024.10.28.pdf", plot = r_squared_binding_score_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = T)


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
  # theme_classic2() +
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
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/RSquared_and_Binding_Energy_2024.10.28.pdf", plot = r_squared_binding_energy_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = T)


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
  # theme_classic2() +
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
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Filtered_RSquared_and_Binding_Energy_2024.10.28.pdf", plot = filtered_r_squared_binding_energy_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)


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
  # theme_classic2() +
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
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Filtered_RSquared_and_Binding_Score_2024.10.28.pdf", plot = filtered_r_squared_binding_score_bubble_plot, width = 10, height = 15, dpi = 900, create.dir = T)




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
  # theme_classic2() +
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
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Heatmaps/Filtered_RSquared_Heatmap_(BS=160_BE=-20)_2024.10.28.pdf", plot = R_squared_heatmap, width = 10, height = 15, dpi = 900, create.dir = T)


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
write_csv3(max_scores_df2, file = 'Tables/Expression_Tables/Venom_Family_Specific_LM_BRM_LME/Top30_Best_miRNA-Gene_Pairs_2024.10.28.csv')

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
  # theme_classic2() +
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
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Heatmaps/Filtered_RSquared_Heatmap_Best_30_2024.10.28.pdf", plot = R_squared_heatmap_best_30, width = 10, height = 15, dpi = 900, create.dir = T)



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
  # theme_classic2() +
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
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Top30_RSquared_and_Binding_Score_2024.10.28.pdf", plot = r_squared_binding_score_bubble_plot_top_30, width = 10, height = 15, dpi = 900, create.dir = T)


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
  # theme_classic2() +
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
ggsave("Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/Residuals/Bubble_Plot/Top30_RSquared_and_Binding_Energy_2024.10.28.pdf", plot = r_squared_binding_energy_bubble_plot_top_30, width = 10, height = 15, dpi = 900, create.dir = T)




#### mRNA vs Protein for all C. viridis ####

# Create new data frame for only the viridis samples
viridis_plot_df <- all_samples_clr_df %>%
  filter(str_detect(Sample.ID, '.viridis.'))


# Set variables for protein and RNA expression levels
viridis_mRNA_expression <- viridis_plot_df$CLR.Scaled.mRNA.Expression
viridis_protein_expression <- viridis_plot_df$CLR.Scaled.Protein.Expression
venom_family <- viridis_plot_df$Venom.Family


# Create figure for mRNA and Protein expression
viridis_only_plot <- viridis_plot_df %>%
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
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
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  # geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
viridis_only_plot

# Save
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/Viridis_Only_Individuals_mRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = viridis_only_plot, create.dir = T)



#### mRNA vs Protein for CV1081 (viridis) ####

# Create new data frame for only the viridis samples
CV1081_df <- all_samples_clr_df %>%
  filter(str_detect(Sample.ID, 'LVG.2.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV1081_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV1081_df$CLR.Scaled.Protein.Expression
venom_family <- CV1081_df$Venom.Family


# Create figure for mRNA and Protein expression
CV1081_plot <- CV1081_df %>%
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
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
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  # geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
CV1081_plot

# Save
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/CV1081_viridis_mRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = CV1081_plot, create.dir = T)


#### mRNA vs Protein for CV0857 (viridis) ####

# Create new data frame for only the viridis samples
CV0857_df <- all_samples_clr_df %>%
  filter(str_detect(Sample.ID, 'LVG.4.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV0857_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV0857_df$CLR.Scaled.Protein.Expression
venom_family <- CV0857_df$Venom.Family


# Create figure for mRNA and Protein expression
CV0857_plot <- CV0857_df %>%
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
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
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  # geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
CV0857_plot

# Save
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/CV0857_viridis_mRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = CV0857_plot, create.dir = T)



#### mRNA vs Protein for CV1086 (viridis) ####

# Create new data frame for only the viridis samples
CV1086_df <- all_samples_clr_df %>%
  filter(str_detect(Sample.ID, 'LVG.9.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV1086_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV1086_df$CLR.Scaled.Protein.Expression
venom_family <- CV1086_df$Venom.Family


# Create figure for mRNA and Protein expression
CV1086_plot <- CV1086_df %>%
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
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
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  # geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
CV1086_plot

# Save
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/CV1086_viridis_mRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = CV1086_plot, create.dir = T)



#### mRNA vs Protein for CV1087 (viridis) ####

# Create new data frame for only the viridis samples
CV1087_df <- all_samples_clr_df %>%
  filter(str_detect(Sample.ID, 'RVG.5S.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV1087_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV1087_df$CLR.Scaled.Protein.Expression
venom_family <- CV1087_df$Venom.Family


# Create figure for mRNA and Protein expression
CV1087_plot <- CV1087_df %>%
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
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
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  # geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
CV1087_plot

# Save
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/CV1087_viridis_mRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = CV1087_plot, create.dir = T)


#### mRNA vs Protein for CV0987 (lutosus) ####

# Create new data frame for only the viridis samples
CV0987_df <- all_samples_clr_df %>%
  filter(str_detect(Sample.ID, 'RVG.6S.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV0987_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV0987_df$CLR.Scaled.Protein.Expression
venom_family <- CV0987_df$Venom.Family


# Create figure for mRNA and Protein expression
CV0987_plot <- CV0987_df %>%
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
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
  # xlim(0, NA) +  # Set the desired x-axis limits
  # ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
CV0987_plot

# Save
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/CV0987_lutosus_mRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = CV0987_plot, create.dir = T)



#### mRNA vs Protein for CV0985 (concolor) ####

# Create new data frame for only the viridis samples
CV0985_df <- all_samples_clr_df %>%
  filter(str_detect(Sample.ID, 'RVG.7S.'))


# Set variables for protein and RNA expression levels
mRNA_expression <- CV0985_df$CLR.Scaled.mRNA.Expression
protein_expression <- CV0985_df$CLR.Scaled.Protein.Expression
venom_family <- CV0985_df$Venom.Family


# Create figure for mRNA and Protein expression
CV0985_plot <- CV0985_df %>%
  ggplot(aes(x = CLR.Scaled.mRNA.Expression, y = CLR.Scaled.Protein.Expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
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
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  # geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
CV0985_plot

# Save
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/CV0985_concolor_mRNAvsProtein_Expression_Plot_2024.10.28.pdf', plot = CV0985_plot, create.dir = T)





#### miRNA vs mRNA expression ####

# Let's check to see how miRNA and mRNA expression levels are correlated with each other
# First thing is to include mRNAs that didn't have any protein
mi_df2 <- read.table(file = miRNA_mRNA_protein_data, header = T) %>%
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
  # Change the name of the venom adam genes so that they are correct
  mutate(
    Genes = ifelse(Genes == "Venom_ADAM28_1", 'Venom_SVMP12', Genes),
    Genes = ifelse(Genes == 'Venom_ADAM28_2', 'ADAM28', Genes)
  )

# Create a data frame for venom genes only
venom_genes_vst_df2 <- mi_df2 %>%
  filter(
    str_starts(Genes, 'Venom_')
  ) %>%
  dplyr::select(
    miRNA.Cluster,
    Genes,
    contains('miRNA.Counts.'),
    contains('RNA.VST.'),
    contains('Intensity.'),
    contains('miRNA.VST.')
  ) %>%
  mutate(Venom.Family = case_when(grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
                                  grepl('VEGF', Genes) ~ 'VEGF',
                                  grepl('ohanin', Genes) ~ 'ohanin',
                                  grepl('vQC', Genes) ~ 'vQC',
                                  grepl('SVSP', Genes) ~ 'SVSP',
                                  grepl('PLA2', Genes) ~ 'PLA2',
                                  grepl('ADAM', Genes) ~ 'ADAM',
                                  grepl('CRISP', Genes) ~ 'CRISP',
                                  grepl('CTL', Genes) ~ 'CTL',
                                  grepl('EXO', Genes) ~ 'EXO',
                                  grepl('LAAO', Genes) ~ 'LAAO',
                                  grepl('myotoxin', Genes) ~ 'myotoxin',
                                  grepl('BPP', Genes) ~ 'BPP',
                                  TRUE ~ 'others')) %>%
  distinct()

# Pivot data frame so that all samples are in a single column
all_samples_plot_protein_df2 <- venom_genes_vst_df2 %>%
  dplyr::select(-contains('miRNA.'), -contains('Intensity')) %>%
  distinct() %>%
  pivot_longer(
    cols = contains('RNA.VST.'),
    names_to = 'Sample.ID',
    values_to = 'RNA.VST'
  ) %>%
  mutate(Sample.ID = gsub('RNA.VST.', '', Sample.ID))

# Pivot this one too
all_samples_plot_mRNA_df2 <- venom_genes_vst_df2 %>%
  dplyr::select(-contains('miRNA.'), -contains('RNA.')) %>%
  distinct() %>%
  pivot_longer(
    cols = contains('Intensity.'),
    names_to = 'Sample.ID',
    values_to = 'Intensity'
  ) %>%
  mutate(Sample.ID = gsub('Intensity.', '', Sample.ID))

# Pivot data that all the miRNA data is the same column
all_sampls_plot_miRNA_df2 <- venom_genes_vst_df2 %>%
  dplyr::select(miRNA.Cluster, Genes, Venom.Family, contains('miRNA.VST')) %>%
  distinct() %>%
  pivot_longer(
    cols = contains('miRNA.VST'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.VST'
  ) %>%
  mutate(Sample.ID = gsub('miRNA.VST.', '', Sample.ID))


# Fuse RNA and protein data back together
all_samples_plot_RNA_protein_df2 <- left_join(
  all_samples_plot_protein_df2,
  all_samples_plot_mRNA_df2,
  by = c('Sample.ID', 'Genes', 'Venom.Family')
)

# Fuse miRNA to the above as well
all_samples_plot_miRNA_RNA_protein_df2 <- left_join(
  all_samples_plot_RNA_protein_df2,
  all_sampls_plot_miRNA_df2,
  by = c('Sample.ID', 'Genes', 'Venom.Family')
) %>%
  distinct()


# Create plot
miRNA_vs_mRNA_all_samples_plot <- ggplot(
  all_samples_plot_miRNA_RNA_protein_df2,
  aes(x = miRNA.VST, y = RNA.VST)
) +
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
    title = 'miRNA vs mRNA',
    y = 'mRNA exression',
    x = 'miRNA expression',
    color = 'Venom Family'
  ) +
  scale_color_manual(values = venom_colors) +  # Apply color scheme
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
miRNA_vs_mRNA_all_samples_plot
# Didn't work... at all
# Zero correlation
# ggsave('Figures/Expression_Plots/Venom_Family_Specific_LM_BRM_LME/miRNA_vs_Protein/Linear_Regression_VST/All_Individuals_mRNAvsmicroRNA_Expression_Plot_2024.10.28.pdf', plot = miRNA_vs_mRNA_all_samples_plot, create.dir = T)





#### Make Figure 6 ####

# Fuse expression plots together
expresion_plots <- plot_grid(
  cluster_1395_svsp7_plot,
  cluster_1395_svmp7_plot,
  let_7c_5p_2_PLA2B1_plot,
  miR_34a_5p_myotoxin_plot,
  cluster_919_svmp4_plot,
  ncol = 1,
  nrow = 5,
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
  ncol = 2,
  nrow = 1
)
bubble_plot_plus_expression_plots

figure_6 <- plot_grid(
  all_gene_expression_plot,
  bubble_plot_plus_expression_plots,
  align = 'v',
  ncol = 1,
  nrow = 2,
  rel_widths = c(1, 3),
  rel_heights = c(1, 3)
)
figure_6
# ggsave("Figures/1_Main_Figures/Figure_6/Figure_6_2024.10.28.pdf", plot = figure_6, width = 18, height = 30, dpi = 2400, create.dir = T)


