# Last Edited: 2025/02/13

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

# CLR Data
clr_data <- 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_venom_2025.01.30.parquet'

# Expression data for venom genes and their miRNAs
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_Venom_2025.01.22.parquet'

# Number of miRNAs per gene data
miRNA_number_data <- 'Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filtered_miRNA_number_data <- 'Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
no_feature_type_data <- 'Data/Merged/No_feature_type_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filt_no_feature_type_data <- 'Data/Merged/No_feature_type_filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'

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
glimpse(mi_df)


# Read in the CLR data
clr_df <- read_parquet(file = clr_data)

# Remove any undetected proteins from the data for fitting purposes
clr_observed_prot_df <- clr_df %>%
  filter(protein.observed == 'Yes')
glimpse(clr_observed_prot_df)

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
# mRNA vs Protein model
venom_model1 <- readRDS('Data/Models/mRNA_vs_Protein/Gaussian/venom_model1_with_500000_iterations_from_brms_2025.01.30.rds')

# mRNA vs Protein model with venom family as an interaction term
venom_model2 <- readRDS('Data/Models/mRNA_vs_Protein/Gaussian/venom_model2_with_500000_iterations_from_brms_2025.01.30.rds')

# mRNA vs Protein model with venom family as an interaction term and miRNA expression as a fixed effect
venom_model3 <- readRDS('Data/Models/mRNA_vs_Protein/Gaussian/venom_model3_with_500000_iterations_from_brms_2025.01.30.rds')

# mRNA vs Protein model with venom family as an interaction term and miRNA expression as a fixed effect, but miRNAs are filtered by score and energy
venom_model4 <- readRDS('Data/Models/mRNA_vs_Protein/Gaussian/venom_model4_with_500000_iterations_from_brms_2025.01.30.rds')

### Compare the venom models ----

# Run waic to see model complexity
waic(venom_model1, venom_model2)

# Output of model 'venom_model1':

# Computed from 4500 by 288 log-likelihood matrix.

#           Estimate   SE
# elpd_waic  -1171.9  8.4
# p_waic         2.9  0.2
# waic        2343.9 16.8

# Output of model 'venom_model2':

# Computed from 4500 by 288 log-likelihood matrix.

#           Estimate   SE
# elpd_waic  -1069.7 16.0
# p_waic        24.2  3.2
# waic        2139.4 31.9

# 13 (4.5%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

# Model comparisons:
#              elpd_diff se_diff
# venom_model2    0.0       0.0 
# venom_model1 -102.2      15.5 
# Warning message:

# 13 (4.5%) p_waic estimates greater than 0.4. We recommend trying loo instead.

# Add loo
add_criterion(
  venom_model1,
  criterion = 'loo'
)

add_criterion(
  venom_model1,
  criterion = 'bayes_R2'
)

add_criterion(
  venom_model2,
  criterion = 'loo'
)

add_criterion(
  venom_model2,
  criterion = 'bayes_R2'
)

# Run loo to see model complexity
loo(venom_model1, venom_model2)
# Output of model 'venom_model1':

# Computed from 4500 by 288 log-likelihood matrix.

#          Estimate   SE
# elpd_loo  -1171.9  8.4
# p_loo         2.9  0.2
# looic      2343.9 16.8
# ------
# MCSE of elpd_loo is 0.0.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.9, 1.0]).

# All Pareto k estimates are good (k < 0.7).
# See help('pareto-k-diagnostic') for details.

# Output of model 'venom_model2':

# Computed from 4500 by 288 log-likelihood matrix.

#          Estimate   SE
# elpd_loo  -1070.4 16.1
# p_loo        24.9  3.4
# looic      2140.8 32.1
# ------
# MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume MCMC draws (r_eff in [0.8, 1.1]).

# Pareto k diagnostic values:
#                          Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     286   99.3%   669     
#    (0.7, 1]   (bad)        2    0.7%   <NA>    
#    (1, Inf)   (very bad)   0    0.0%   <NA>    
# See help('pareto-k-diagnostic') for details.

# Model comparisons:
#              elpd_diff se_diff
# venom_model2    0.0       0.0 
# venom_model1 -101.5      15.6 
# Warning message:
# Found 2 observations with a pareto_k > 0.7 in model 'venom_model2'. We recommend to set 'moment_match = TRUE' in order to perform moment matching for problematic observations.

# Summarize model 1
venom_model1_summary <- summary(venom_model1)
venom_model1_summary
#  Family: gaussian
#   Links: mu = identity; sigma = identity 
# Formula: protein.clr ~ mRNA.clr + (1 | gr(sample.id, cov = phylo_cov_matrix)) 
#    Data: clr_df (Number of observations: 288) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.69      0.66     0.02     2.35 1.00     4372     4490

# Regression Coefficients:
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept     0.02      1.00    -1.95     1.94 1.00     4438     4573
# mRNA.clr      1.69      0.11     1.48     1.91 1.00     4548     4576

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma    14.07      0.59    12.99    15.26 1.00     4303     4447

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# Summarize model 2
venom_model2_summary <- summary(venom_model2)
venom_model2_summary
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: protein.clr ~ mRNA.clr * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)) 
#    Data: clr_df (Number of observations: 288) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.61      0.60     0.02     2.27 1.00     4306     4406

# Regression Coefficients:
#                               Estimate Est.Error l-95% CI u-95% CI Rhat
# Intercept                       -10.11      3.31   -16.48    -3.68 1.00
# mRNA.clr                          3.54      0.63     2.30     4.76 1.00
# venom.familyCRISP                -1.28      3.72    -8.59     5.91 1.00
# venom.familyCTL                   2.48      3.88    -4.96    10.01 1.00
# venom.familyEXO                 -11.72      4.28   -20.02    -3.27 1.00
# venom.familyLAAO                  5.26      3.84    -2.25    12.79 1.00
# venom.familymyotoxin             -5.46      6.99   -19.22     8.01 1.00
# venom.familyohanin                3.05      8.87   -14.52    20.14 1.00
# venom.familyPLA2                  2.24      4.23    -6.28    10.12 1.00
# venom.familySVMP                  8.43      3.50     1.57    15.30 1.00
# venom.familySVSP                 15.23      4.72     5.78    24.54 1.00
# venom.familyVEGF                 -2.89      4.13   -11.05     4.97 1.00
# venom.familyvQC                 -15.15      4.98   -24.87    -5.33 1.00
# mRNA.clr:venom.familyCRISP       -2.66      0.68    -3.98    -1.32 1.00
# mRNA.clr:venom.familyCTL         -2.73      0.65    -3.98    -1.41 1.00
# mRNA.clr:venom.familyEXO         -1.73      3.80    -9.28     5.72 1.00
# mRNA.clr:venom.familyLAAO        -1.83      0.69    -3.16    -0.47 1.00
# mRNA.clr:venom.familymyotoxin     1.07      1.26    -1.34     3.59 1.00
# mRNA.clr:venom.familyohanin       1.75      2.83    -3.82     7.21 1.00
# mRNA.clr:venom.familyPLA2         0.39      0.90    -1.35     2.19 1.00
# mRNA.clr:venom.familySVMP         0.77      0.74    -0.68     2.21 1.00
# mRNA.clr:venom.familySVSP        -1.03      1.14    -3.28     1.21 1.00
# mRNA.clr:venom.familyVEGF         3.51      1.52     0.48     6.51 1.00
# mRNA.clr:venom.familyvQC          8.79      2.08     4.74    12.78 1.00
#                               Bulk_ESS Tail_ESS
# Intercept                         4486     4486
# mRNA.clr                          4547     4445
# venom.familyCRISP                 4553     4295
# venom.familyCTL                   4314     4523
# venom.familyEXO                   4465     4187
# venom.familyLAAO                  4646     4445
# venom.familymyotoxin              4596     3793
# venom.familyohanin                4334     4244
# venom.familyPLA2                  4526     4548
# venom.familySVMP                  4411     4329
# venom.familySVSP                  4576     3871
# venom.familyVEGF                  4545     4390
# venom.familyvQC                   4712     4367
# mRNA.clr:venom.familyCRISP        4543     4347
# mRNA.clr:venom.familyCTL          4494     4492
# mRNA.clr:venom.familyEXO          4615     4619
# mRNA.clr:venom.familyLAAO         4383     4263
# mRNA.clr:venom.familymyotoxin     4109     4162
# mRNA.clr:venom.familyohanin       4591     4447
# mRNA.clr:venom.familyPLA2         4560     4343
# mRNA.clr:venom.familySVMP         4449     4365
# mRNA.clr:venom.familySVSP         4638     4436
# mRNA.clr:venom.familyVEGF         4383     4144
# mRNA.clr:venom.familyvQC          4757     4436

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     9.46      0.41     8.70    10.30 1.00     4405     4318

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

### Compare the models that use miRNAs ----

# Compare the models
waic(venom_model3)
waic(venom_model4)

# Summarize model 3
venom_model3_summary <- summary(venom_model3)
venom_model3_summary
# Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: protein.clr ~ mRNA.clr * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)) + miRNA.rpm 
#    Data: miRNA_clr_df (Number of observations: 6834) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     3.43      1.26     1.88     6.61 1.00     4442     4331

# Regression Coefficients:
#                               Estimate Est.Error l-95% CI u-95% CI Rhat
# Intercept                       -10.94      3.65   -18.09    -3.80 1.00
# mRNA.clr                          3.65      0.40     2.84     4.43 1.00
# venom.familyCRISP                -7.23      3.07   -13.44    -1.30 1.00
# venom.familyCTL                  -5.28      3.07   -11.58     0.69 1.00
# venom.familyEXO                 -12.48      3.09   -18.85    -6.57 1.00
# venom.familyLAAO                  5.70      3.10    -0.55    11.69 1.00
# venom.familymyotoxin             -6.91      3.35   -13.45    -0.37 1.00
# venom.familyohanin               14.22      5.87     2.61    25.50 1.00
# venom.familyPLA2                  4.32      3.15    -1.94    10.38 1.00
# venom.familySVMP                  6.86      3.06     0.56    12.71 1.00
# venom.familySVSP                 18.97      3.10    12.78    25.01 1.00
# venom.familyVEGF                 -2.29      3.10    -8.52     3.64 1.00
# venom.familyvQC                 -16.01      3.28   -22.53    -9.53 1.00
# miRNA.rpm                         0.02      0.01     0.00     0.04 1.00
# mRNA.clr:venom.familyCRISP       -3.26      0.41    -4.04    -2.44 1.00
# mRNA.clr:venom.familyCTL         -3.23      0.40    -4.01    -2.41 1.00
# mRNA.clr:venom.familyEXO         -3.38      0.77    -4.86    -1.87 1.00
# mRNA.clr:venom.familyLAAO        -1.99      0.41    -2.79    -1.18 1.00
# mRNA.clr:venom.familymyotoxin     1.22      0.46     0.31     2.11 1.00
# mRNA.clr:venom.familyohanin      -1.36      1.62    -4.54     1.80 1.00
# mRNA.clr:venom.familyPLA2         0.19      0.44    -0.66     1.06 1.00
# mRNA.clr:venom.familySVMP         1.36      0.41     0.57     2.17 1.00
# mRNA.clr:venom.familySVSP        -1.92      0.44    -2.78    -1.07 1.00
# mRNA.clr:venom.familyVEGF         3.51      0.48     2.59     4.46 1.00
# mRNA.clr:venom.familyvQC          9.83      0.65     8.57    11.14 1.00
#                               Bulk_ESS Tail_ESS
# Intercept                         4747     4367
# mRNA.clr                          4591     4490
# venom.familyCRISP                 4562     4405
# venom.familyCTL                   4589     4503
# venom.familyEXO                   4517     4406
# venom.familyLAAO                  4590     4566
# venom.familymyotoxin              4529     4534
# venom.familyohanin                4640     4225
# venom.familyPLA2                  4561     4245
# venom.familySVMP                  4593     4616
# venom.familySVSP                  4514     4531
# venom.familyVEGF                  4597     4488
# venom.familyvQC                   4564     4403
# miRNA.rpm                         4402     4476
# mRNA.clr:venom.familyCRISP        4603     4483
# mRNA.clr:venom.familyCTL          4583     4449
# mRNA.clr:venom.familyEXO          4497     4547
# mRNA.clr:venom.familyLAAO         4613     4622
# mRNA.clr:venom.familymyotoxin     4503     4441
# mRNA.clr:venom.familyohanin       4443     4273
# mRNA.clr:venom.familyPLA2         4566     4331
# mRNA.clr:venom.familySVMP         4606     4577
# mRNA.clr:venom.familySVSP         4454     4530
# mRNA.clr:venom.familyVEGF         4690     4515
# mRNA.clr:venom.familyvQC          4409     4448

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     8.17      0.07     8.03     8.31 1.00     4616     4290

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# Summarize model 3
venom_model4_summary <- summary(venom_model4)
venom_model4_summary
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: protein.clr ~ mRNA.clr * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)) + miRNA.rpm 
#    Data: filt_miRNA_clr_df (Number of observations: 1446) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     3.47      1.41     1.67     7.09 1.00     4180     4407

# Regression Coefficients:
#                               Estimate Est.Error l-95% CI u-95% CI Rhat
# Intercept                       -12.55      3.74   -20.12    -5.29 1.00
# mRNA.clr                          3.88      0.50     2.89     4.83 1.00
# venom.familyCRISP                -6.99      3.30   -13.30    -0.31 1.00
# venom.familyCTL                   1.84      3.22    -4.34     8.44 1.00
# venom.familyEXO                 -10.66      3.34   -17.10    -3.98 1.00
# venom.familyLAAO                  6.48      3.48    -0.23    13.38 1.00
# venom.familymyotoxin             -4.63      4.29   -12.80     3.75 1.00
# venom.familyohanin                4.12      8.95   -13.60    21.72 1.00
# venom.familyPLA2                 -4.73      4.49   -13.46     4.24 1.00
# venom.familySVMP                  7.32      3.19     1.24    13.79 1.00
# venom.familySVSP                 20.63      3.68    13.47    28.06 1.00
# venom.familyVEGF                 -0.63      3.32    -7.04     6.12 1.00
# venom.familyvQC                 -12.58      4.59   -21.49    -3.39 1.00
# miRNA.rpm                         0.09      0.03     0.03     0.15 1.00
# mRNA.clr:venom.familyCRISP       -3.59      0.51    -4.60    -2.57 1.00
# mRNA.clr:venom.familyCTL         -2.87      0.51    -3.86    -1.87 1.00
# mRNA.clr:venom.familyEXO         -3.33      1.66    -6.59    -0.04 1.00
# mRNA.clr:venom.familyLAAO        -2.22      0.53    -3.27    -1.21 1.00
# mRNA.clr:venom.familymyotoxin     0.92      0.70    -0.47     2.27 1.00
# mRNA.clr:venom.familyohanin       1.80      2.79    -3.64     7.25 1.00
# mRNA.clr:venom.familyPLA2         1.87      0.80     0.31     3.45 1.00
# mRNA.clr:venom.familySVMP         1.45      0.52     0.45     2.46 1.00
# mRNA.clr:venom.familySVSP        -2.04      0.72    -3.47    -0.63 1.00
# mRNA.clr:venom.familyVEGF         3.32      0.73     1.91     4.76 1.00
# mRNA.clr:venom.familyvQC          8.87      1.63     5.62    12.03 1.00
#                               Bulk_ESS Tail_ESS
# Intercept                         4556     4341
# mRNA.clr                          4588     4331
# venom.familyCRISP                 4716     4330
# venom.familyCTL                   4586     4450
# venom.familyEXO                   4616     4464
# venom.familyLAAO                  4675     4189
# venom.familymyotoxin              4618     4422
# venom.familyohanin                4050     3952
# venom.familyPLA2                  4624     4083
# venom.familySVMP                  4621     4408
# venom.familySVSP                  4480     4491
# venom.familyVEGF                  4541     4407
# venom.familyvQC                   4580     4450
# miRNA.rpm                         4462     4492
# mRNA.clr:venom.familyCRISP        4567     4254
# mRNA.clr:venom.familyCTL          4542     4451
# mRNA.clr:venom.familyEXO          4574     4195
# mRNA.clr:venom.familyLAAO         4615     4410
# mRNA.clr:venom.familymyotoxin     4862     4661
# mRNA.clr:venom.familyohanin       4464     4338
# mRNA.clr:venom.familyPLA2         4725     4621
# mRNA.clr:venom.familySVMP         4554     4446
# mRNA.clr:venom.familySVSP         4311     4664
# mRNA.clr:venom.familyVEGF         4662     4481
# mRNA.clr:venom.familyvQC          4601     4441

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     8.80      0.16     8.48     9.12 1.00     4552     4400

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


# Plot and examine venom model 1 ----

## Extract important information form the brm modle using tidybayes (venom_model1) ----

# Extract the data
tidy_1 <- tidy(venom_model1)

# Create a data frame to contain epred (expected values of the response variable that don't include observation noise) draws
model1_df <- clr_df %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(
    venom_model1, 
    .draw = TRUE, 
    .chain = TRUE, 
    .iteration = TRUE
    # ndraws = 50
  )
glimpse(model1_df)

# Create a new data frame that contains the means of the epreds
model1_mean_epred_df <- model1_df %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = protein.clr - mean.epred)
glimpse(model1_mean_epred_df)

# Create a data frame to hold the R^2 values
r2_model1 <- clr_df %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model1, .draw = TRUE) %>%
  select(
    -.chain, -.iteration
  ) %>%
  distinct() %>%
  summarise(
    mean.epred = mean(.epred) # summarise all of the predicted values to obtain a single point estimate for each observed value - predicted value pair
  ) %>%
  select(-.row) %>%
  mutate(residuals = protein.clr - mean.epred) %>%
  ungroup() %>%
  summarise(
    R2 = 1 - sum((protein.clr - mean.epred)^2) / sum((protein.clr - mean(protein.clr))^2)
  )
glimpse(r2_model1)

# Create a data frame to hold the R^2 values
r2_model1_2 <- clr_df %>%
  # Epreds are expected values of the response variable that don't include observation noise
  tidybayes::add_epred_draws(venom_model1, .draw = TRUE) %>%
  distinct() %>%
  mutate(
    residual.draws = protein.clr - .epred
  ) %>%
  ungroup() %>%
  summarise(
    var.fit = var(.epred),
    var.res = var(residual.draws)
  ) %>%
  mutate(
    R2 = (var.fit) / (var.fit + var.res)
  )
r2_model1_2


## Create Plots ----

# Create figure for mRNA and Protein expression for the simple brms model 
mRNA_vs_protein_plot <- ggplot(data = clr_df, aes(x = mRNA.clr, y = protein.clr)) +
  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(clr_df, protein.observed == "No"),
    aes(fill = venom.family),
    shape = 21, size = 2.5, show.legend = FALSE
  ) +
  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(clr_df, protein.observed == "Yes"),
    aes(fill = venom.family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line
  geom_line(
    data = model1_df,
    aes(y = .epred, group = .draw),
    linewidth = 1,
    linetype = 'solid',
    show.legend = FALSE,
    alpha = 1/4
  ) +  # Add epreds lines

  # Add R² text
  geom_text_repel(
    data = r2_model1_2,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = TRUE, show.legend = FALSE
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
mRNA_vs_protein_plot
ggsave('Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Main/Venom_model1_mRNA_vs_protein_clr_expression_line_cluster.2025.02.13.png', plot = mRNA_vs_protein_plot, height = 15, width = 20, create.dir = TRUE)

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
mRNA_vs_protein_plot_ribbon <- ggplot(data = clr_df, aes(x = mRNA.clr, y = protein.clr)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(clr_df, protein.observed == "No"),
    aes(fill = venom.family),
    shape = 21, size = 2.5, show.legend = FALSE
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(clr_df, protein.observed == "Yes"),
    aes(fill = venom.family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add regression line from venom_model8 with proper grouping and fill for venom.family
  ggdist::stat_lineribbon(
    data = model1_df,
    aes(y = .epred),  # Color ribbons by venom.family
    .width = c(.99, .95, .8, .5),
    alpha = 0.1, # Set transparency for the ribbon
    show.legend = FALSE
  ) +

  # Add R² text
  geom_text_repel(
    data = r2_model1_2,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = TRUE, show.legend = FALSE
  ) +

  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  labs(
    title = 'a. mRNA vs Protein',
    y = 'clr(peak intensity)',
    x = 'clr(gene expression)',
    fill = 'Venom Family'
  ) +
  # scale_color_manual(values = venom_colors) +  # Apply color scheme
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
mRNA_vs_protein_plot_ribbon
ggsave('Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Main/Venom_model1_mRNA_vs_protein_clr_expression_ribbon.2025.02.13.png', plot = mRNA_vs_protein_plot_ribbon, height = 15, width = 20, create.dir = TRUE)


# Create a mean epred plot
mRNA_vs_protein_plot_mean <- ggplot(data = clr_df, aes(x = mRNA.clr, y = protein.clr)) +

  # Plot all points that were not used for line fitting
  geom_point(
    data = subset(clr_df, protein.observed == "No"),
    aes(fill = venom.family),
    shape = 21, size = 2.5, show.legend = FALSE
  ) +

  # Add points that were used in fitting with a black outline
  geom_point(
    data = subset(clr_df, protein.observed == "Yes"),
    aes(fill = venom.family),
    size = 2.5, shape = 21, color = "black", stroke = 1
  ) +  # Outline with black

  # Add mean predictions
  geom_line(data = model1_mean_epred_df, aes(y = mean.epred), linewidth = 1, linetype = 'solid', show.legend = FALSE) +  # Add mean predicted lines

  # Add R² text
  geom_text_repel(
    data = r2_model1_2,
    aes(x = -Inf, y = Inf, label = paste("R² =", round(R2, 3))),
    hjust = 1.2, vjust = 1.2, size = 3, # color = "black",
    inherit.aes = TRUE, show.legend = FALSE
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
  # scale_color_manual(values = venom_colors) +  # Apply color scheme to points/lines
  scale_fill_manual(values = venom_colors) +  # Apply color scheme for ribbons
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
mRNA_vs_protein_plot_mean
ggsave('Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Main/Venom_model1_mRNA_vs_protein_clr_expression_mean.2025.02.13.png', plot = mRNA_vs_protein_plot_mean, height = 15, width = 20, create.dir = TRUE)



## Format residuals data ----

# Create name intersection for data
names <- intersect(names(model1_mean_epred_df), names(clr_df))
names2 <- intersect(names(mi_df), names(clr_df))

# Combine the data frames
model1_residuals_df <- left_join(
  clr_df,
  model1_mean_epred_df,
  by = names
) %>%
  left_join(
    mi_df,
    by = names2
  ) %>%
  select(
    -feature.type, -miRNA.counts, -miRNA.vst, -total.score, -total.energy
  ) %>%
  distinct()
glimpse(model1_residuals_df)


## miRNA vs residual protein expression Graphs (based on VST) ----

# Get the number of available cores
num_cores <- detectCores()

# Source plotting functions
source('Scripts/R/Functions/Residuals_Correlation_Function.R')

# Create a list of venom genes for a function to iterate through
venom_genes <- model1_residuals_df %>% distinct(genes)
venom_genes <- venom_genes$genes

# Create path for the plots to go in
path <- 'Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Residual_Plots/All'

# Add only the plots with R squared higher than
# Create a second path
path2 <- 'Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Residual_Plots/Good_R2'

# Define your batch size
batch_size <- 5  # Adjust this depending on your system's memory limits

# Split genes into batches
gene_batches <- split(venom_genes, ceiling(seq_along(venom_genes)/batch_size))

# Parallelize the execution of batches and do all plots
mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot3(model1_residuals_df, gene, path, '2025.02.13', filter_r_squared = FALSE)
    gc()
  }
}, mc.cores = num_cores)

# For good R squared values
mclapply(gene_batches, function(gene_batch) {
  for (gene in gene_batch) {
    # Run your plot function for each gene
    residuals_vs_miRNA_plot3(model1_residuals_df, gene, path2, '2025.02.13', filter_r_squared = TRUE)
    gc()
  }
}, mc.cores = num_cores)


## Faceted residuals plots ----

# # Create sub-data frames for a subsets of genes
# svsp_df  <- model1_residuals_df %>%
#   filter(venom.family %in% c('SVSP'))

# # Create a figure where each gene is a facet, and the miRNAs are plotted on the same graph
# protein_residuals_plot <- ggplot(svsp_df, aes(x = miRNA.rpm, y = residuals, group = miRNA.cluster)) +
#   geom_point(aes(color = miRNA.cluster)) +
#   geom_smooth(
#     aes(color = miRNA.cluster),
#     method = 'lm',
#     se = FALSE,
#     linetype = 'dashed',
#     formula = y ~ x,
#     show.legend = FALSE
#   ) +
#   stat_poly_eq(
#     aes(
#       color = miRNA.cluster,
#       label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")
#     ),
#     formula = y ~ x
#   ) +
#   scale_x_continuous() +  # Set x-axis to continuous scale
#   scale_y_continuous() +  # Set y-axis to continuous scale
#   labs(
#     title = 'miRNA expression vs Protein Residuals',
#     y = 'Protein Residuals',
#     x = 'miRNA RPM',
#     color = 'Venom Family'
#   ) +
#   theme_linedraw() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 10, hjust = 0.5),
#     legend.position = 'none',
#     legend.title.position = 'none'
#   ) +
#   facet_wrap( ~ genes)
# protein_residuals_plot
# # This will be impossible to parse, unfortunately, time to abbandon it

## Bubble plots ----

### Format data for the bubble plots ----
# Calculate the covariance
covariance_df <- model1_residuals_df %>%
  ungroup() %>%
  select(sample.id, genes, venom.family, residuals, miRNA.cluster, miRNA.rpm) %>%
  group_by(venom.family, genes, miRNA.cluster) %>%
  summarise(
    pearson.cov = cov(miRNA.rpm, residuals, method = 'pearson'),
    kendall.cov = cov(miRNA.rpm, residuals, method = 'kendall'),
    spearman.cov = cov(miRNA.rpm, residuals, method = 'spearman'),
    pearson.cor = cor(miRNA.rpm, residuals, method = 'pearson')
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

# Create a data frame that contains residauls and other important data
miRNA_covariance_df <- left_join(
  covariance_df,
  mi_df2,
  by = shared_column_names,
  relationship = 'many-to-many'
) %>%
  dplyr::select(
    genes, venom.family, contains('cov'), pearson.cor, miRNA.cluster, total.score, total.energy, feature.type
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
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Covariance/miRNA_expression_covariance_against_protein_residuals_vs_binding_energy.2025.02.14.pdf", plot = be_vs_covariance_plot, width = 10, height = 10, dpi = 900, create.dir = TRUE)

# Same as the above, but labeled with the miRNA cluster and gene the dot represents
be_vs_covariance_plot_labeled <- be_vs_covariance_plot +
  geom_text(
    aes(label = paste(genes, miRNA.cluster, sep = ": ")), # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  )
be_vs_covariance_plot_labeled
# They do not correlate at all
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Covariance/Labeled_miRNA_expression_covariance_against_protein_residuals_vs_binding_energy.2025.02.14.pdf", plot = be_vs_covariance_plot_labeled, width = 10, height = 10, dpi = 900, create.dir = TRUE)


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
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Covariance/miRNA_expression_covariance_against_protein_residuals_vs_binding_score.2025.02.14.pdf", plot = bs_vs_covariance_plot, width = 10, height = 10, dpi = 900, create.dir = TRUE)


# Same as the above, but labeled with the miRNA cluster and gene the dot represents
bs_vs_covariance_plot_labeled <- bs_vs_covariance_plot +
  geom_text(
    aes(label = paste(genes, miRNA.cluster, sep = ": ")), # Combine Genes and miRNA.Cluster in label
    hjust = 0, vjust = -1,  # Adjust the position of the labels
    check_overlap = TRUE  # Avoid overlapping labels
  )
bs_vs_covariance_plot_labeled
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Covariance/Labeled_miRNA_expression_covariance_against_protein_residuals_vs_binding_score.2025.02.14.pdf", plot = bs_vs_covariance_plot_labeled, width = 10, height = 10, dpi = 900, create.dir = TRUE)


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
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Correlation/miRNA_residauls_correlation_and_bs_bubble_plot_2025.02.14.pdf", plot = r_and_binding_score_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)


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
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Correlation/miRNA_residauls_correlation_and_be_bubble_plot_2025.02.14.pdf", plot = r_and_binding_energy_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)

### Format data for the filtered bubble plots ----

# Set limits for grey and red color scale
min_correlation = 0.5
max_total_energy = -7
min_total_score = 155

miRNA_covariance_df2 <- miRNA_covariance_df %>%
  filter(
    abs(pearson.cor) > min_correlation,
    total.energy <= max_total_energy,
    min_total_score >= 155
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
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Correlation/filtered_miRNA_residauls_correlation_and_bs_bubble_plot_2025.02.14.pdf", plot = filt_r_and_binding_score_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)


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
ggsave("Figures/Expression_Plots/BRMS/mRNA_vs_Protein/Guassian/Residuals/Venom_Model1/Correlation/filtered_miRNA_residauls_correlation_and_be_bubble_plot_2025.02.14.pdf", plot = filt_r_and_binding_energy_bubble_plot, width = 15, height = 28, dpi = 900, create.dir = TRUE)
