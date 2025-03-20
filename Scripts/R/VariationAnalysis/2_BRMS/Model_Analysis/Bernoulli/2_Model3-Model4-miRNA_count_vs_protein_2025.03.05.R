# Last Edited: 2025/03/05

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

## Set data ----
# Set working directory
setwd("/home/administrator/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA")

# Set path for the main data
miRNA_number_data <- "Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet"
filtered_miRNA_number_data <- "Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet"

## Read miRNA number data ----
miRNA_num_df <- read_parquet(file = miRNA_number_data) %>%
  filter(
    !str_detect(genes, "maker-scaffold|augustus|XP_|ADAM28"),
    in.library == "Yes",
    str_detect(genes, "Venom_"),
    feature.type == "three_prime_utr"
  ) %>%
  # Create the binary expression columns
  mutate(
    mRNA.binary = ifelse(mRNA.counts > 0, 1, 0),
    protein.binary = ifelse(intensity > 0, 1, 0),
    # Create a second binary column to see if it improves the model
    mRNA.binary2 = ifelse(mRNA.counts > 1000, 1, 0)
  )

# Create a variable to stor the mean value for the mRNA counts
mean_mRNA <- mean(miRNA_num_df$mRNA.counts)

# Create a variable to store the mean value for the protein intensity
mean_protein <- mean(miRNA_num_df$intensity)

# Create new vectors in the data that represent miRNA and protein expression as binary values
miRNA_num_df <- miRNA_num_df %>%
  mutate(
    mRNA.level = ifelse(mRNA.counts > mean_mRNA, 1, 0),
    protein.level = ifelse(intensity > mean_protein, 1, 0)
  )


## Read filtered miRNA number data ----
filt_miRNA_num_df <- read_parquet(file = filtered_miRNA_number_data) %>%
  filter(
    !sample.id == "CV1082_viridis",
    !str_detect(genes, "maker-scaffold|augustus|XP_|ADAM28"),
    in.library == "Yes",
    str_detect(genes, "Venom_"),
    feature.type == "three_prime_utr"
  ) %>%
  # Create the binary expression columns
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
    protein.level = ifelse(intensity > mean_protein, 1, 0)
  )


## Create color scheme for the venom genes ----
SVMP_color <- "#4A70B5"
ADAM_color <- "#9A70B5"
SVSP_color <- "#F0B830"
PLA2_color <- "#7570B3"
miRNA_color <- "#8B0AA5"
VEGF_color <- "#74ADD1"
ohanin_color <- "#3A489C"
myotoxin_color <- "#B2182B"
vQC_color <- "#80BC50"
CRISP_color <- "#E7298A"
CTL_color <- "#F67E17"
EXO_color <- "#005824"
LAAO_color <- "#B35806"
BPP_color <- "#1B9E77"
other_color <- "#666666"
three_prime_color <- "black"
# five_prime_color <- '#0072b2'
five_prime_color <- "#1B9E77"
# cds_color <- '#d55e00'
cds_color <- "#4A70B5"

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
# protein vs miRNA numbers model
model3 <- readRDS(
  "Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model3_with_500000_iterations_2025.03.02.rds"
)

# mRNA vs miRNA numbers model, filtered
model4 <- readRDS(
  "Data/Models/mRNA_vs_Number_of_miRNAs/Bernoulli/Model4_with_500000_iterations_2025.03.02.rds"
)

### Summarize the models ----
# Summarize model3
model3_summary <- summary(model3)
model3_summary
#  Family: bernoulli 
#   Links: mu = logit 
# Formula: protein.binary ~ number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)) 
#    Data: miRNA_num_df (Number of observations: 185) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.31      0.31     0.01     1.13 1.00     4597     4289

# Regression Coefficients:
#                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            2.39      0.41     1.63     3.20 1.00     4481     4327
# number.of.miRNAs    -0.17      0.03    -0.24    -0.11 1.00     4198     4533

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# Summarize model4
model4_summary <- summary(model4)
model4_summary
#  Family: bernoulli 
#   Links: mu = logit 
# Formula: protein.binary ~ number.of.miRNAs + (1 | gr(sample.id, cov = phylo_cov_matrix)) 
#    Data: filt_miRNA_num_df (Number of observations: 116) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.44      0.40     0.02     1.47 1.00     4494     4707

# Regression Coefficients:
#                  Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept            3.74      0.71     2.47     5.25 1.00     4694     4488
# number.of.miRNAs    -0.79      0.15    -1.10    -0.52 1.00     4722     4620

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


# Plot and examine venom model 3 ----
## Format the input data ----

## Get tidy dat for the model ----
tidy3 <- broom.mixed::tidy(model3)
tidy3


## Get epreds ----
model3_df <- miRNA_num_df %>%
  # Get epreds
  tidybayes::add_epred_draws(
    model3,
    .draw = TRUE,
    .chain = TRUE,
    .iteration = TRUE,
    # Return all draws with ndraws = NULL
    ndraws = NULL
  )
glimpse(model3_df)

## Plot conditional effects ----

# Create a plot for the conditional effects of the entire model
model3_conditional_effects <- conditional_effects(model3)
model3_conditional_effects


## Plot effects with ggplot ----

# Create a plot with tidybayes output for the model
model3_plot <- ggplot(
  model3_df,
  aes(
    x = number.of.miRNAs,
    y = protein.binary
  )
) +
  # Add the points
  geom_point(
    aes(y = protein.binary, color = venom.family),
    shape = 16
  ) +
  # Get rid of the gray box around the venom family legend points
  guides(color = guide_legend(override.aes = list(fill = NA, shape = 16))) +
  # Add the tidybayes ribbon
  ggdist::stat_lineribbon(
    aes(y = .epred),
    .width = c(0.99, 0.95, 0.8, 0.5),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  scale_color_manual(values = venom_colors) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  scale_y_continuous() +
  labs(
    title = "Probability of venom protein expression",
    x = "Number of targeting miRNAs",
    y = "Probability of venom protein expression",
    color = "Venom Family"
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    legend.title.position = "top",
    axis.title = element_text(face = "italic")
  )
model3_plot
ggsave("Figures/Expression_Plots/BRMS/miRNA_numbers_vs_mRNA/Bernoulli/model3_plot_2025.03.05.png", plot = model3_plot, width = 10, height = 10, create.dir = TRUE)


## Try plotting with a regular glm ----
# Create a plot with tidybayes output for the model
glm_model3_plot <- ggplot(
  miRNA_num_df,
  aes(
    x = number.of.miRNAs,
    y = protein.binary
  )
) +
  geom_point(aes(y = protein.binary, color = venom.family)) +
  geom_smooth(
    formula = y ~ x,
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE,
    color = "blue",
    fill = "lightblue"
  ) +
  scale_color_manual(values = venom_colors) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(
    title = "Probability of venom protein expression",
    x = "Number of targeting miRNAs",
    y = "Probability of venom protein expression",
    color = "Venom Family"
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    legend.title.position = "top",
    axis.title = element_text(face = "italic")
  )
glm_model3_plot
ggsave("Figures/Expression_Plots/GLM/miRNA_numbers_vs_mRNA/Bernoulli/model3_plot_2025.03.05.png", plot = glm_model3_plot, width = 10, height = 10, create.dir = TRUE)



# Plot and examine venom model 4 ----
## Format the input data ----

## Get tidy dat for the model ----
tidy4 <- broom.mixed::tidy(model4)
tidy4

## Get epreds ----
model4_df <- filt_miRNA_num_df %>%
  # Get epreds
  tidybayes::add_epred_draws(
    model4,
    .draw = TRUE,
    .chain = TRUE,
    .iteration = TRUE,
    # Return all draws with ndraws = NULL
    ndraws = NULL
  )
glimpse(model4_df)

## Plot conditional effects ----

# Create a plot for the conditional effects of the entire model
model4_conditional_effects <- conditional_effects(model4)
model4_conditional_effects


## Plot effects with ggplot ----

# Create a plot with tidybayes output for the model
model4_plot <- ggplot(
  model4_df,
  aes(
    x = number.of.miRNAs,
    y = protein.level
  )
) +
  # Add the points
  geom_point(
    aes(y = protein.level, color = venom.family),
    shape = 16
  ) +
  # Get rid of the gray box around the venom family legend points
  guides(color = guide_legend(override.aes = list(fill = NA, shape = 16))) +
  # Add the tidybayes ribbon
  ggdist::stat_lineribbon(
    aes(y = .epred),
    .width = c(0.99, 0.95, 0.8, 0.5),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  scale_color_manual(values = venom_colors) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 12)) +
  scale_y_continuous() +
  labs(
    title = "Probability of venom protein expression",
    x = "Number of targeting miRNAs",
    y = "Probability of venom protein expression",
    color = "Venom Family"
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    legend.title.position = "top",
    axis.title = element_text(face = "italic")
  )
model4_plot
ggsave("Figures/Expression_Plots/BRMS/miRNA_numbers_vs_mRNA/Bernoulli/model4_plot_2025.03.05.png", plot = model4_plot, width = 10, height = 10, create.dir = TRUE)

## Try plotting with a regular glm ----
# Create a plot with tidybayes output for the model
glm_model4_plot <- ggplot(
  filt_miRNA_num_df,
  aes(
    x = number.of.miRNAs,
    y = protein.level
  )
) +
  geom_point(aes(y = protein.level, color = venom.family)) +
  geom_smooth(
    formula = y ~ x,
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE,
    color = "blue",
    fill = "lightblue"
  ) +
  scale_color_manual(values = venom_colors) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(
    title = "Probability of venom protein expression",
    x = "Number of targeting miRNAs",
    y = "Probability of venom protein expression",
    color = "Venom Family"
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "bottom",
    legend.title.position = "top",
    axis.title = element_text(face = "italic")
  )
glm_model4_plot
ggsave("Figures/Expression_Plots/GLM/miRNA_numbers_vs_mRNA/Bernoulli/model4_plot_2025.03.05.png", plot = glm_model4_plot, width = 10, height = 10, create.dir = TRUE)
