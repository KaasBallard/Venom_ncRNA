# Last Edited: 2025/03/02

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

# Number of miRNAs per gene data
miRNA_number_data <- 'Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filtered_miRNA_number_data <- 'Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
no_feature_type_data <- 'Data/Merged/No_feature_type_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filt_no_feature_type_data <- 'Data/Merged/No_feature_type_filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'

## Read miRNA number data ----
miRNA_num_df <- read_parquet(file = miRNA_number_data) %>%
  filter(
    # !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    feature.type == 'three_prime_utr'
  ) %>%
  mutate(
    mRNA.binary = ifelse(mRNA.counts > 0, 1, 0),
    protein.binary = ifelse(intensity > 0, 1, 0)
  )
glimpse(miRNA_num_df)

# Create a variable to stor the mean value for the mRNA counts
mean_mRNA <- mean(miRNA_num_df$mRNA.counts)

# Create a variable to store the mean value for the protein intensity
mean_protein <- mean(miRNA_num_df$intensity)

# Create new vectors in the data that represent miRNA and protein expression as binary values
miRNA_num_df <- miRNA_num_df %>%
  mutate(
    mRNA.level = ifelse(mRNA.counts > mean_mRNA, 1, 0),
    protein.level = ifelse(intensity > mean_protein, 1, 0),
    mRNA.level2 = if_else(mRNA.counts > 1000, 1, 0)
  )

## Read miRNA number data (no feature type version) ----
no_ft_miRNA_num_df <- read_parquet(file = no_feature_type_data) %>%
  filter(
    # !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_')
  ) %>%
  mutate(
    mRNA.binary = ifelse(mRNA.counts > 0, 1, 0),
    protein.binary = ifelse(intensity > 0, 1, 0)
  )

# Create a variable to stor the mean value for the mRNA counts
mean_mRNA <- mean(no_ft_miRNA_num_df$mRNA.counts)

# Create a variable to store the mean value for the protein intensity
mean_protein <- mean(no_ft_miRNA_num_df$intensity)

# Create new vectors in the data that represent miRNA and protein expression as binary values
no_ft_miRNA_num_df <- no_ft_miRNA_num_df %>%
  mutate(
    mRNA.level = ifelse(mRNA.counts > mean_mRNA, 1, 0),
    protein.level = ifelse(intensity > mean_protein, 1, 0),
    mRNA.level2 = if_else(mRNA.counts > 1000, 1, 0)
  )

## Read filtered miRNA number data ----
filt_miRNA_num_df <- read_parquet(file = filtered_miRNA_number_data) %>%
  filter(
    # !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    feature.type == 'three_prime_utr'
  ) %>%
  mutate(
    mRNA.binary = ifelse(mRNA.counts > 0, 1, 0),
    protein.binary = ifelse(intensity > 0, 1, 0)
  )

# Create a variable to stor the mean value for the mRNA counts
mean_mRNA <- mean(filt_miRNA_num_df$mRNA.counts)

# Create a variable to store the mean value for the protein intensity
mean_protein <- mean(filt_miRNA_num_df$intensity)

# Create new vectors in the data that represent miRNA and protein expression as binary values
filt_miRNA_num_df <- filt_miRNA_num_df %>%
  mutate(
    mRNA.level = ifelse(mRNA.counts > mean_mRNA, 1, 0),
    protein.level = ifelse(intensity > mean_protein, 1, 0),
    mRNA.level2 = if_else(mRNA.counts > 1000, 1, 0)
  )

## Read filtered miRNA number data (no feature type version) ----
no_ft_filt_miRNA_num_df <- read_parquet(file = filt_no_feature_type_data) %>%
  filter(
    # !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_')
  ) %>%
  mutate(
    mRNA.binary = ifelse(mRNA.counts > 0, 1, 0),
    protein.binary = ifelse(intensity > 0, 1, 0)
  )

# Create a variable to stor the mean value for the mRNA counts
mean_mRNA <- mean(no_ft_filt_miRNA_num_df$mRNA.counts)

# Create a variable to store the mean value for the protein intensity
mean_protein <- mean(no_ft_filt_miRNA_num_df$intensity)

# Create new vectors in the data that represent miRNA and protein expression as binary values
no_ft_filt_miRNA_num_df <- no_ft_filt_miRNA_num_df %>%
  mutate(
    mRNA.level = ifelse(mRNA.counts > mean_mRNA, 1, 0),
    protein.level = ifelse(intensity > mean_protein, 1, 0),
    mRNA.level2 = if_else(mRNA.counts > 1000, 1, 0)
  )

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
mismatched_samples <- setdiff(miRNA_num_df$sample.id, tree_tips)

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
  print(conditional_effects(model, ask = FALSE))
  dev.off()
}


# BRMS modeling ----
## Set brms model parameters and priors ----
warmup_percentage <- 0.25
chain <- 6
# thin <- 500
thin <- 500
iterations <- 500000
iterations2 <- 25000 # Temporarily reduce iterations to decrease run time for debugging purposes
warmup <- iterations * warmup_percentage
warmup2 <- iterations2 * warmup_percentage

# Set priors
priors <- c(
  # Set prior for effect of number of miRNAs
  # The negative normal(-1, 2) means I suspect that the number of miRNAs will have a negative effect on mRNA and protein expression
  # The -1 (negative because the correlation will be negative) represents the mean of the prior and the 2 represents the standard deviation
  # The first number represents the mean of the prior and the second number represents the standard deviation

  # Fixed effects for the mu portion of the model
  # The negative normal(-1, 2) means I suspect that the number of miRNAs will have a negative effect on mRNA and protein expression
  prior(normal(0, 5), class = "b"),

  # Set the prior for the intercept (the mean of the prior is 0 and the standard deviation is 3)
  prior(normal(0, 5), class = "Intercept"),

  # Set the prior for random effects (first number is the location parameter and the second number is the scale parameter)
  prior(cauchy(0, 2), class = "sd") # Random effects standard deviation
)

# Diagnostics directory
diagnostics_path <- sprintf("Figures/Model_Diagnostics/Bernoulli/Inter_%d_2025.03.02", iterations)

## Bernoulli models ----
### Venom model #1: mRNA vs miRNA numbers ----
#### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
model1 <- brm(
  mRNA.binary ~
    number.of.miRNAs + (1 | r(sample.id, cov = phylo_cov_matrix)),
  data = miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model1 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model1,
  directory = diagnostics_path,
  model_number = '1'
)

#### Save the model ----
# Create a filename for the model
model1_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model1_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model1, file = model1_filename)
rm(model1)
gc()

### Venom model #2: mRNA vs filtered miRNA numbers ----
#### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
model2 <- brm(
  mRNA.binary ~
    number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = filt_miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model2 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model2,
  directory = diagnostics_path,
  model_number = '2'
)

#### Save the model ----
# Create a filename for the model
model2_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model2_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model2, file = model2_filename)
rm(model2)
gc()

### Venom model #3: number of miRNAs vs Protein ----
#### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
model3 <- brm(
  protein.binary ~
    number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model3 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model3,
  directory = diagnostics_path,
  model_number = '3'
)

#### Save the model ----
# Create a filename for the model
model3_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model3_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model3, file = model3_filename)
rm(model3)
gc()

### Venom model #4: Protein vs filtered miRNA numbers ----
#### Run the model ----
# Get time
t1 <- Sys.time()

# Run the model
model4 <- brm(
  protein.binary ~
    number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = filt_miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model4 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model4,
  directory = diagnostics_path,
  model_number = '4'
)

#### Save the model ----
# Create a filename for the model
model4_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model4_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model4, file = model4_filename)
rm(model4)
gc()



### Venom model #5: mRNA vs miRNA numbers ----
#### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
model5 <- brm(
  mRNA.level ~
    number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model5 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model5,
  directory = diagnostics_path,
  model_number = '5'
)

#### Save the model ----
# Create a filename for the model
model5_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model5_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model5, file = model5_filename)
rm(model5)
gc()

### Venom model #6: mRNA vs filtered miRNA numbers ----
#### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
model6 <- brm(
  mRNA.level ~
    number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = filt_miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model6 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model6,
  directory = diagnostics_path,
  model_number = '6'
)

#### Save the model ----
# Create a filename for the model
model6_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model6_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model6, file = model6_filename)
rm(model6)
gc()

### Venom model #7: number of miRNAs vs Protein ----
#### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
model7 <- brm(
  protein.level ~
    number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model7 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model7,
  directory = diagnostics_path,
  model_number = '7'
)

#### Save the model ----
# Create a filename for the model
model7_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model7_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model7, file = model7_filename)
rm(model7)
gc()

### Venom model #8: Protein vs filtered miRNA numbers ----
#### Run the model ----
# Get time
t1 <- Sys.time()

# Run the model
model8 <- brm(
  protein.level ~
    number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)),
  data = filt_miRNA_num_df,
  family = bernoulli(link = 'logit'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
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
# Get time
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed model8 with ', iterations, ' iterations in ', timer))

#### Model Diagnostics ----

# Run the diagnostics function
model_diagnostics(
  model = model8,
  directory = diagnostics_path,
  model_number = '8'
)

#### Save the model ----
# Create a filename for the model
model8_filename <- sprintf('Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model8_with_%d_iterations_2025.03.02.rds', iterations)
# Save the model
saveRDS(model8, file = model8_filename)
rm(model8)
gc()