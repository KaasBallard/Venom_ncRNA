# Last Edited: 2024/10/25

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
library(cmdstanr)
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






#### Gene family plots data frame for Venom Gene Dot Plots (VST) ####

# Create a data frame for venom genes only
venom_genes_vst_df <- mi_df %>%
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




#### Build data frame for future analysis ####


# Pivot data frame so that all samples are in a single column
all_samples_plot_protein_df <- venom_genes_vst_df %>%
  dplyr::select(-contains('miRNA.'), -contains('Intensity')) %>%
  distinct() %>%
  pivot_longer(
    cols = contains('RNA.VST.'),
    names_to = 'Sample.ID',
    values_to = 'RNA.VST'
  ) %>%
  mutate(Sample.ID = gsub('RNA.VST.', '', Sample.ID)) %>% 
  # Add species so that the shapes can be changed in the plot
  mutate(
    Species = case_when(
      grepl('.viridis.', Sample.ID) ~ 'C.viridis',
      grepl('.lutosus.', Sample.ID) ~ 'C.lutosus',
      grepl('.concolor.', Sample.ID) ~ 'C.concolor'
    )
  )

# Pivot this one too
all_samples_plot_mRNA_df <- venom_genes_vst_df %>%
  dplyr::select(-contains('miRNA.'), -contains('RNA.')) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('Intensity.'),
    names_to = 'Sample.ID',
    values_to = 'Intensity'
  ) %>% 
  mutate(Sample.ID = gsub('Intensity.', '', Sample.ID)) %>% 
  # Add species so that the shapes can be changed in the plot
  mutate(
    Species = case_when(
      grepl('.viridis.', Sample.ID) ~ 'C.viridis',
      grepl('.lutosus.', Sample.ID) ~ 'C.lutosus',
      grepl('.concolor.', Sample.ID) ~ 'C.concolor'
    )
  )


# Pivot data that all the miRNA data is the same column
all_sampls_plot_miRNA_df <- venom_genes_vst_df %>% 
  dplyr::select(miRNA.Cluster, Genes, Venom.Family, contains('miRNA.VST')) %>%
  distinct() %>%
  pivot_longer(
    cols = contains('miRNA.VST'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.VST'
  ) %>%
  mutate(Sample.ID = gsub('miRNA.VST.', '', Sample.ID)) %>% 
  # Add species so that the shapes can be changed in the plot
  mutate(
    Species = case_when(
      grepl('.viridis.', Sample.ID) ~ 'C.viridis',
      grepl('.lutosus.', Sample.ID) ~ 'C.lutosus',
      grepl('.concolor.', Sample.ID) ~ 'C.concolor'
    )
  )



# Fuse RNA and protein data back together
all_samples_plot_RNA_protein_df <- left_join(
  all_samples_plot_protein_df,
  all_samples_plot_mRNA_df,
  by = c('Sample.ID', 'Genes', 'Venom.Family', 'Species')
) %>% 
  # Add a column to indicate their was zero protein expression detected for this data point
  mutate(
    Protein.Detected = case_when(
      Intensity == 0 ~ 'No',
      Intensity > 0 ~ 'Yes'
    )
  )

# Fuse miRNA to the above as well
all_samples_plot_miRNA_RNA_protein_df <- left_join(
  all_samples_plot_RNA_protein_df,
  all_sampls_plot_miRNA_df,
  by = c('Sample.ID', 'Genes', 'Venom.Family', 'Species')
)
# Save data for future use
write_csv(all_samples_plot_miRNA_RNA_protein_df, file = 'Data/miRNA/miRNA_mRNA_protein_data_for_regression_analysis_2024.10.27.csv')


# Apparently zero values are not treated properly by the clr transform, so I am going to create a pseudo-count to fix it.
# First I need to find the lowest value in each set of data though
# Find the lowest non-zero Intensity value
min_nonzero_intensity <- all_samples_plot_RNA_protein_df %>%
  filter(Intensity > 0) %>%
  summarise(min_value = min(Intensity)) %>%
  pull(min_value)
print(min_nonzero_intensity)
# Find the lowest non-zero RNA.VST value
min_nonzero_rna_vst <- all_samples_plot_RNA_protein_df %>%
  filter(RNA.VST > 0) %>%
  summarise(min_value = min(RNA.VST)) %>%
  pull(min_value)
print(min_nonzero_rna_vst)

# Initialize pseudo-count for mRNA
mRNA_pseudo_count <- 1e-6
# Initialize pseudo-count for protein
protein_pseduo_count <- 1e-6





# Apparently zero values are not treated properly by the clr transform, so I have already filtered out proteins that were not been observed.
# Calculate CLR
all_samples_clr_df <- all_samples_plot_RNA_protein_df %>%
  # Replace zeros with the small constant before CLR
  mutate(Intensity = ifelse(Intensity == 0, protein_pseduo_count, Intensity)) %>%
  mutate(RNA.VST = ifelse(RNA.VST == 0, mRNA_pseudo_count, RNA.VST)) %>%
  group_by(Sample.ID) %>%
  mutate(Scaled.Protein.Expression = Intensity / sum(Intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.mRNA.Expression = RNA.VST / sum(RNA.VST)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  # Perform CLR transformation
  mutate(CLR.Scaled.Protein.Expression = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(CLR.Scaled.mRNA.Expression = as.numeric(compositions::clr(Scaled.mRNA.Expression)))
# Save the above as a data frame to be read in another script
write_csv(all_samples_clr_df, file = 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_2024.10.27.csv')


# Remove any undetected proteins from the data for fitting purposes
all_samples_detected_proteins_only_df <- all_samples_clr_df %>% 
  filter(Protein.Detected == 'Yes')


# Check for mismatches between tree tips and Sample.ID
mismatched_samples <- setdiff(all_samples_detected_proteins_only_df$Sample.ID, tree_tips)

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




#### Venom modeling ####

# Set brms model parameters
# Run 3 or 4 times
warmup_percentage <- 0.25 # Set the percent to be kept by the model
chain <- 6
iter <- 5000000 # Add a zero
thin <- 500 # make this at least 100, maybe 500
warmup <- iter * warmup_percentage # Keep 50%

priors <- c(
  prior(normal(0, 10), class = "b"),       # Fixed effects
  prior(normal(0, 10), class = "Intercept"), # Intercept
  prior(cauchy(0, 1), class = "sd"),       # Random effects standard deviation
  prior(cauchy(0, 2), class = "sigma")     # Residual standard deviation, allowing for high variability
)

# Number of runs
n_runs <- 10

# Loop to run the model 10 times
for (i in 1:n_runs) {
  t1 <- Sys.time() # Start timer
  
  # Run the model
  venom_model8 <- brm(
    CLR.Scaled.Protein.Expression ~ CLR.Scaled.mRNA.Expression * Venom.Family + (1|gr(Sample.ID, cov = phylo_cov_matrix)),
    data = all_samples_detected_proteins_only_df,
    family = gaussian(link = 'identity'),
    data2 = list(phylo_cov_matrix = phylo_cov_matrix),
    control = list(max_treedepth = 12), # Increase tree depth
    chains = chain,
    iter = iter,
    thin = thin,
    warmup = warmup,
    prior = priors,
    save_pars = save_pars(all = T),
    cores = parallel::detectCores(),
    # Use cmdstanr because it is supposedly faster
    backend = "cmdstanr"
  )
  
  # End timer and calculate elapsed time
  t2 <- Sys.time() # End timer
  timer <- t2 - t1
  print(paste('Run', i, 'completed in', timer))
  
  # Save each model with a unique file name
  saveRDS(venom_model8, file = paste0("Data/Models/venom_model8.", i, "_from_brms_2024.10.28.rds"))
}



