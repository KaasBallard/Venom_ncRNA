# Last Edited: 2025/01/30

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
library(ggrepel)
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
library(MuMIn)
library(broom.mixed)
library(lme4)
library(cmdstanr)
library(arrow)
# Source my functions
# source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# Set working directory
setwd("/home/administrator/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA")

## Set data ----

# CLR Data
clr_data <- 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_venom_2025.01.30.parquet'
# Expression data for venom genes and their miRNAs
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_Venom_2025.01.22.parquet'

# Number of miRNAs per gene data
miRNA_number_data <- 'Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filtered_miRNA_number_data <- 'Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
no_feature_type_data <- 'Data/Merged/No_feature_type_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filt_no_feature_type_data <- 'Data/Merged/No_feature_type_filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'


## Read reference data ----

# Read reference data in as a DataFrame
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  dplyr::filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_')
  ) %>% 
  dplyr::select(
    sample.id, genes, venom.family, mRNA.counts, mRNA.vst, intensity, protein.observed, in.library, total.energy, total.score, feature.type, miRNA.cluster, miRNA.counts, miRNA.rpm, miRNA.vst
  )
rm(miRNA_mRNA_protein_df)

## Read miRNA number data ----
miRNA_num_df <- read_parquet(file = miRNA_number_data) %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    feature.type == 'three_prime_utr'
  )

## Read miRNA number data (no feature type version) ----
no_ft_miRNA_num_df <- read_parquet(file = no_feature_type_data) %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_')
  )


## Read filtered miRNA number data ----
filt_miRNA_num_df <- read_parquet(file = filtered_miRNA_number_data) %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    feature.type == 'three_prime_utr'
  )

## Read filtered miRNA number data (no feature type version) ----
no_ft_filt_miRNA_num_df <- read_parquet(file = filt_no_feature_type_data) %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_')
  )


## Read CLR data ----
clr_df <- read_parquet(file = clr_data)

# Remove any undetected proteins from the data for fitting purposes
clr_observed_prot_df <- clr_df %>%
  filter(protein.observed == 'Yes')
glimpse(clr_observed_prot_df)


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

# Check for mismatched samples between the tree and data ----
# Check for mismatches between tree tips and sample.id
mismatched_samples <- setdiff(clr_observed_prot_df$sample.id, tree_tips)

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


# BRMS modeling ----

## Venom model #1: mRNA vs Protein ----

# Set brms model parameters
warmup_percentage <- 0.25
chain <- 6
thin <- 500
iterations <- 500000
iterations2 <- 25000 # Temporarily reduce iterations to decrease run time for debugging purposes
warmup <- iterations * warmup_percentage
warmup2 <- iterations2 * warmup_percentage

# Set priors
priors <- c(
  prior(normal(0, 10), class = "b"), # Fixed effects
  prior(normal(0, 10), class = "Intercept"), # Intercept
  prior(cauchy(0, 1), class = "sd"), # Random effects standard deviation
  prior(cauchy(0, 2), class = "sigma") # Residual standard deviation, allowing for high variability
)

# Get time
t1 <- Sys.time()

# TODO: Run this as a hurdle model to account for zeros (run it as an additional model)
# TODO: Increase adapt_delta to above 0.9
# Run the model
venom_model1 <- brm(
  protein.clr ~
    mRNA.clr + (1|gr(sample.id, cov = phylo_cov_matrix)),
  data = clr_df,
  family = gaussian(link = 'identity'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model1 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model1, file = 'Data/Models/mRNA_vs_Protein/Gaussian/venom_model1_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model1)
gc()


## Venom model #2: mRNA vs Protein by family ----

# Get time
t1 <- Sys.time()

# TODO: set moment match to TRUE (moment_match = TRUE)
# Run the model
venom_model2 <- brm(
  protein.clr ~
    mRNA.clr * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix)),
  data = clr_df,
  family = gaussian(link = 'identity'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model2 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model2, file = 'Data/Models/mRNA_vs_Protein/Gaussian/venom_model2_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model2)
gc()


## Venom model #3: mRNA vs Protein + miRNA expression ----

### Format the data ----

# Get shared columns
shared_columns <- intersect(names(mi_df), names(clr_df))

# Fuse miRNA data to the clr data
miRNA_clr_df <- left_join(
  mi_df,
  clr_df,
  by = shared_columns
) %>% 
  select(
    -feature.type, -total.score, -total.energy
  ) %>% 
  distinct()
# any(is.na(miRNA_clr_df))
glimpse(miRNA_clr_df)


### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
venom_model3 <- brm(
  protein.clr ~
    mRNA.clr * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix)) + miRNA.rpm,
  data = miRNA_clr_df,
  family = gaussian(link = 'identity'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model3 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model3, file = 'Data/Models/mRNA_vs_Protein/Gaussian/venom_model3_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model3)
gc()


## Venom model #4: mRNA vs Protein + miRNA expression, filtered by score and energy ----

### Format the data ----

# Get shared columns
shared_columns <- intersect(names(mi_df), names(clr_df))

# Fuse miRNA data to the clr data
filt_miRNA_clr_df <- left_join(
  mi_df,
  clr_df,
  by = shared_columns
) %>%
  filter(
    total.energy <= -7,
    total.score >= 155
  ) %>% 
  select(
    -feature.type, -total.score, -total.energy
  ) %>% 
  distinct()
# any(is.na(filt_miRNA_clr_df))
glimpse(filt_miRNA_clr_df)


### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
venom_model4 <- brm(
  protein.clr ~
    mRNA.clr * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix)) + miRNA.rpm,
  data = filt_miRNA_clr_df,
  family = gaussian(link = 'identity'),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4), 
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model4 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model4, file = 'Data/Models/mRNA_vs_Protein/Gaussian/venom_model4_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model4)
gc()



## Venom model #5: mRNA vs miRNA numbers ----

### Format the data ----

# Get shared columns
shared_columns <- intersect(names(miRNA_num_df), names(clr_df))

# Fuse miRNA data to the clr data
miRNA_num_clr_df <- left_join(
  miRNA_num_df,
  clr_df,
  by = shared_columns
) %>%
  select(
    -feature.type
  ) %>% 
  distinct()
any(is.na(miRNA_num_clr_df))
glimpse(miRNA_num_clr_df)


### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
venom_model5 <- brm(
  bf(
    mRNA.ntd ~
      number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix)),
    hu ~ number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix))
  ),
  data = miRNA_num_clr_df,
  family = hurdle_lognormal(),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model5 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model5, file = 'Data/Models/mRNA_vs_Number_of_miRNAs/Hurdle_Lognormal/venom_model5_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model5)
gc()



## Venom model #6: mRNA vs number of miRNAs + filtered miRNA numbers ----

### Format the data ----

# Get shared columns
shared_columns <- intersect(names(filt_miRNA_num_df), names(clr_df))

# Fuse miRNA data to the clr data
filt_miRNA_num_clr_df <- left_join(
  filt_miRNA_num_df,
  clr_df,
  by = shared_columns
) %>%
  select(
    -feature.type
  ) %>% 
  distinct()
any(is.na(miRNA_num_clr_df))
glimpse(miRNA_num_clr_df)


### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
venom_model6 <- brm(
  bf(
    mRNA.ntd ~
      number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix)),
    hu ~ number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix))
  ),
  data = filt_miRNA_num_clr_df,
  family = hurdle_lognormal(),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model6 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model6, file = 'Data/Models/mRNA_vs_Number_of_miRNAs/Hurdle_Lognormal/venom_model6_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model6)
gc()

## Venom model #7: number of miRNAs vs Protein ----

### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
venom_model7 <- brm(
  bf(
    intensity ~
      number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix)),
    hu ~ number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix))
  ),
  data = miRNA_num_clr_df,
  family = hurdle_lognormal(),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model7 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model7, file = 'Data/Models/mRNA_vs_Number_of_miRNAs/Hurdle_Lognormal/venom_model7_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model7)
gc()


## Venom model #8: Protein vs number of miRNAs + filtered miRNA numbers ----

### Run the model ----

# Get time
t1 <- Sys.time()

# Run the model
venom_model8 <- brm(
  bf(
    intensity ~
      number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix)),
    hu ~ number.of.miRNAs * venom.family + (1|gr(sample.id, cov = phylo_cov_matrix))
  ),
  data = filt_miRNA_num_clr_df,
  family = hurdle_lognormal(),
  data2 = list(phylo_cov_matrix = phylo_cov_matrix),
  control = list(max_treedepth = 12, adapt_delta = 0.9), # Increase tree depth
  chains = chain,
  iter = iterations,
  thin = thin,
  warmup = warmup,
  prior = priors,
  save_pars = save_pars(all = T),
  cores = chain,
  threads = threading(4),
  # Use cmdstanr because it is supposedly faster
  backend = "cmdstanr"
)
t2 <- Sys.time()
timer <- t2 - t1
print(paste('Completed venom_model8 with ', iterations, ' iterations in ', timer))

# Save the model
saveRDS(venom_model8, file = 'Data/Models/mRNA_vs_Number_of_miRNAs/Hurdle_Lognormal/venom_model8_with_500000_iterations_from_brms_2025.01.30.rds')
rm(venom_model8)
gc()
