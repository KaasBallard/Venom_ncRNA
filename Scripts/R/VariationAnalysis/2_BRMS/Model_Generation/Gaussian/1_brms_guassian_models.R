# Last Edited: 2025/03/13

# Set up and Read Data ----

## Load packages ----

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
library(readxl)
library(ggpmisc)
library(gganimate)
library(gifski)
library(modelr)
library(scales)
# Load parallel package for the loop later on
library(parallel)
# Load libraries that are important for the tree model
library(brms)
library(caper)
library(visreg)
library(ape)
library(plyr)
library(phytools)
library(ggtree)
library(tidytree)
library(tidybayes)
library(bayesplot)
library(MuMIn)
library(broom.mixed)
library(lme4)
library(cmdstanr)
library(arrow)
# Set working directory
setwd("/home/administrator/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA")

## Set data ----

# Expression data for venom genes and their miRNAs
mRNA_protein_data <- 'Data/Merged/mRNA_Protein_Combined_Data_2025.01.22.parquet'
file.exists(mRNA_protein_data)

## Read in expression data ----
mRNA_protein_df <- read_parquet(file = mRNA_protein_data) %>%
  dplyr::filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    # Remove any instances where a protein was not observed
    protein.observed == 'Yes',
    # Filter out any proteins that truely were not observed, as in 0 intensity
    intensity > 0
  ) %>%
  dplyr::mutate(
    log.protein = log(intensity + 1)
  ) %>%
  select(
    genes, sample.id, mRNA.vst, mRNA.ntd, venom.family, log.protein
  )
glimpse(mRNA_protein_df)

## Tree data ----
# Tree file for the 12 snakes samples
tree_file <- 'Data/Tree/11Snake_Tree_ViridisPolytomy.phy'

# Read tree file in
snake_tree <- ape::read.tree(tree_file)

## Colors ----
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

# Format the phylogeny and get the matrix ----

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
  "CV0857_viridis_North_M" = 'CV0857_viridis',
  "CV1081_viridis_Mid_M" = 'CV1081_viridis',
  "CV1087_viridis_North_F" = 'CV1087_viridis',
  "CV1086_viridis_South_M" = 'CV1086_viridis',
  "CV0985_concolor_Other_F" = 'CV0985_concolor',
  "CV0987_lutosus_Other_F" = 'CV0987_lutosus'
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
# tree_plot

# Get tip labels from the phylogenetic tree
tree_tips <- snake_tree$tip.label

# Generate the phylogenetic covariance matrix from the tree
phylo_cov_matrix <- ape::vcv.phylo(snake_tree)
eigen(phylo_cov_matrix)$values  # All eigenvalues should be positive

# Check for mismatched samples between the tree and data
# Check for mismatches between tree tips and sample.id
mismatched_samples <- setdiff(mRNA_protein_df$sample.id, tree_tips)

# Output mismatches (if any)
mismatched_samples

# Check for missing or infinite values in the covariance matrix
any_na_inf <- any(is.na(phylo_cov_matrix) | is.infinite(phylo_cov_matrix))

# Output the result
any_na_inf

# Check the structure of the tree
str(snake_tree)

# Ensure there are valid branch lengths
summary(snake_tree$edge.length)

## Functions ----
# This function creates a set of diagnostics for the model it is given
model_diagnostics <- function(model, directory, model_number) {
  # Check if needed packages are installed and load them
  if (!require("brms", quietly = TRUE)) {
    install.packages("brms")
    library(brms)
  }
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if (!require("bayesplot", quietly = TRUE)) {
    install.packages("bayesplot")
    library(bayesplot)
  }

  # Create directory for current model
  model_dir <- file.path(directory, sprintf("Model%s_Diagnostics", model_number))
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

  # Model Diagnostics Title
  cat(sprintf("\n\n--- Model %s Diagnostics ---\n", model_number))

  # Model Summary
  model_summary_text <- capture.output(summary(model))
  cat("\nModel Summary:\n")
  print(summary(model))
  # Create a file for summary
  file <- paste0(model_dir, "/summary.txt")
  writeLines(model_summary_text, file)

  # Posterior Summary
  model_posterior_summary_text <- capture.output(posterior_summary(model))
  cat("\nPosterior Summary:\n")
  print(posterior_summary(model))
  # Create a file for Posterior Summary
  file <- paste0(model_dir, "/posterior_summary.txt")
  writeLines(model_posterior_summary_text, file)

  # Print the rhat values
  model_Rhat_text <- capture.output(brms::rhat(model)) # Call the correct function
  cat("\nRhat values:\n")
  print(brms::rhat(model))
  # Create a file for Rhat
  file <- paste0(model_dir, "/rhat.txt")
  writeLines(model_Rhat_text, file)

  # Calculate the effective sample size
  model_neff_ratio_text <- capture.output(neff_ratio(model))
  cat("\nEffective Sample Size:\n")
  print(neff_ratio(model))
  # Create a file for ESS
  file <- paste0(model_dir, "/neff_ratio.txt")
  writeLines(model_neff_ratio_text, file)

  # Posterior exploration
  pdf_file_plot <- file.path(model_dir, "plot.pdf")
  pdf(file = pdf_file_plot, width = 8.5, height = 11)
  print(plot(model))
  dev.off()

  # Posterior predictive checks
  pdf_file_pp_check <- file.path(model_dir, "pp_check.pdf")
  pdf(file = pdf_file_pp_check, width = 8.5, height = 11)
  print(pp_check(model))
  dev.off()

  # Check the mcmc info
  pdf_file_parcoord <- file.path(model_dir, "mcmc_parcoord.pdf")
  pdf(file = pdf_file_parcoord, width = 8.5, height = 11)
  print(bayesplot::mcmc_parcoord(as.array(model), pars = vars(starts_with("b"))))
  dev.off()
  pdf_file_trace <- file.path(model_dir, "mcmc_trace.pdf")
  pdf(file = pdf_file_trace, width = 8.5, height = 11)
  print(bayesplot::mcmc_trace(as.array(model), pars = vars(starts_with("b"))))
  dev.off()

  pdf_file_dens <- file.path(model_dir, "mcmc_dens.pdf")
  pdf(file = pdf_file_dens, width = 8.5, height = 11)
  print(bayesplot::mcmc_dens(as.array(model), pars = vars(starts_with("b"))))
  dev.off()

  # Check conditional effects
  pdf_file_conditionl <- file.path(model_dir, "conditional_effects.pdf")
  pdf(file = pdf_file_conditionl, width = 8.5, height = 11)
  print(conditional_effects(model))
  dev.off()
}


# BRMS modeling ----

## Set brms model parameters and priors ----
warmup_percentage <- 0.25
chain <- 6
thin <- 500
iterations <- 500000
iterations2 <- 25000 # Temporarily reduce iterations to decrease run time for debugging purposes
warmup <- iterations * warmup_percentage
warmup2 <- iterations2 * warmup_percentage

# Set priors
priors <- c(
  # For the normal() functions, the first argument is the mean and the second is the standard deviation
  # For the cauchy() functions, the first argument is the location and the second is the scale
  # Fixed effects
  prior(normal(0, 10), class = "b"),

  # Intercept
  prior(normal(0, 10), class = "Intercept"),

  # Random effects standard deviation
  prior(cauchy(0, 2), class = "sd"),

  # Residual standard deviation, allowing for high variability
  prior(cauchy(0, 5), class = "sigma")
)

# Set a second set of priors for the student-t distribution model
priors_t <- c(
  # For the normal() functions, the first argument is the mean and the second is the standard deviation
  # For the cauchy() functions, the first argument is the location and the second is the scale
  # Fixed effects
  prior(normal(0, 10), class = "b"),

  # Intercept
  prior(normal(0, 10), class = "Intercept"),

  # Random effects standard deviation
  prior(cauchy(0, 2), class = "sd"),

  # Residual standard deviation, allowing for high variability
  prior(cauchy(0, 5), class = "sigma"),

  prior(gamma(2, 0.1), class = "nu")  # Prior for degrees of freedom parameter in t-distribution
)

# Priors for heteroscedastic model
priors_het <- c(
  prior(normal(0, 10), class = "b"),
  prior(normal(0, 10), class = "Intercept"),
  prior(cauchy(0, 2), class = "sd"),
  prior(normal(0, 3), class = "b", dpar = "sigma")  # Prior for the sigma coefficients
)

# Diagnostics directory
diagnostics_path <- sprintf("Scripts/R/VariationAnalysis/2_BRMS/Model_Diagnostics/Gaussian/Inter_%d_2025.03.13", iterations)

## Guassian models ----
### Venom model #1: mRNA vs Protein by family ----
#### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
model1 <- brm(
  log.protein ~
    mRNA.vst * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = mRNA_protein_df,
  family = gaussian(link = 'identity'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.95), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = TRUE),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model1 with ', iterations, ' iterations in ', timer))

##### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model1,
  directory = diagnostics_path,
  model_number = '1'
)

#### Save the model ----
# Create a filename for the model
model1_filename <- sprintf('Data/Models/mRNA_vs_Protein/Gaussian/Model1_with_%d_iterations_2025.03.13.rds', iterations)
# Save the model
saveRDS(model1, file = model1_filename)
rm(model1)
gc()


### Venom model #2: mRNA vs Protein with student-t distribution ----
#### Run the model ----
# Get time
t1 <- Sys.time()

# Run the model
model2 <- brm(
  log.protein ~
    mRNA.vst * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = mRNA_protein_df,
  family = student(link = 'identity'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.95), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors_t,
  save_pars = save_pars(all = TRUE),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model2 with ', iterations, ' iterations in ', timer))

##### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model2,
  directory = diagnostics_path,
  model_number = '2'
)

#### Save the model ----
# Create a filename for the model
model2_filename <- sprintf('Data/Models/mRNA_vs_Protein/Gaussian/Model2_with_%d_iterations_2025.03.13.rds', iterations)
# Save the model
saveRDS(model2, file = model2_filename)
rm(model2)
gc()


### Venom model #3: mRNA vs Protein allowing heteroscedasticity ----
#### Run the model ----
# Get time
t1 <- Sys.time()

# Run the model
model3 <- brm(
  bf(log.protein ~ mRNA.vst * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)), sigma ~ 0 + venom.family),  # Allows different variances for each venom family
  data = mRNA_protein_df,
  family = gaussian(link = 'identity'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.95), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors_het,
  save_pars = save_pars(all = TRUE),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model3 with ', iterations, ' iterations in ', timer))

##### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model3,
  directory = diagnostics_path,
  model_number = '3'
)

#### Save the model ----
# Create a filename for the model
model3_filename <- sprintf('Data/Models/mRNA_vs_Protein/Gaussian/Model3_with_%d_iterations_2025.03.13.rds', iterations)
# Save the model
saveRDS(model3, file = model3_filename)
rm(model3)
gc()
