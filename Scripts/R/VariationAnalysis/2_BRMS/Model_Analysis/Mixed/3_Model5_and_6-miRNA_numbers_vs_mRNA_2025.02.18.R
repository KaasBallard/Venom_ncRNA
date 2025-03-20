# Last Edited: 2025/02/18

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
miRNA_number_data <- 'Data/Merged/miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'
filtered_miRNA_number_data <- 'Data/Merged/filtered_miRNA_numbers_per_sample_mRNA-Protein.2025.01.22.parquet'

# Set path for the clr transformed data
clr_data <- 'Data/CLR_transformed_data/CLR_mRNA_vs_protein_expression_venom_2025.01.30.parquet'

## Read miRNA number data ----
miRNA_num_df <- read_parquet(file = miRNA_number_data) %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM28'),
    in.library == 'Yes',
    str_detect(genes, 'Venom_'),
    feature.type == 'three_prime_utr'
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

## Read in clr data ----
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
# mRNA vs miRNA numbers model
venom_model5 <- readRDS('Data/Models/mRNA_vs_Number_of_miRNAs/Hurdle_Lognormal/venom_model5_with_500000_iterations_from_brms_2025.01.30.rds')

# mRNA vs miRNA numbers model, filtered
venom_model6 <- readRDS('Data/Models/mRNA_vs_Number_of_miRNAs/Hurdle_Lognormal/venom_model6_with_500000_iterations_from_brms_2025.01.30.rds')

### Summarize the models ----
# Summarize model5
model5_summary <- summary(venom_model5)
model5_summary
# Family: hurdle_lognormal 
#   Links: mu = identity; sigma = identity; hu = logit 
# Formula: mRNA.ntd ~ number.of.miRNAs * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)) 
#          hu ~ number.of.miRNAs * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix))
#    Data: miRNA_num_clr_df (Number of observations: 185) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#                   Estimate Est.Error  l-95% CI  u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)         0.08      0.06      0.01      0.24 1.20       23       50
# sd(hu_Intercept) 262194.19 103685.41 121167.73 490108.00 1.81        9       32

# Regression Coefficients:
#                                             Estimate  Est.Error    l-95% CI
# Intercept                                       2.61       0.39        1.85
# hu_Intercept                               373626.22  117127.44   212474.15
# number.of.miRNAs                               -0.08       0.03       -0.14
# venom.familyCTL                                -2.94       0.67       -4.05
# venom.familyEXO                                 0.71       3.55       -5.76
# venom.familymyotoxin                            2.65       9.46      -13.47
# venom.familyohanin                             -0.25       9.45      -17.35
# venom.familyPLA2                               -0.31       0.39       -1.18
# venom.familySVMP                               -0.06       0.39       -0.80
# venom.familySVSP                                0.07       0.42       -0.71
# venom.familyVEGF                               -1.34       8.97      -16.22
# venom.familyvQC                                 0.08       0.72       -1.39
# number.of.miRNAs:venom.familyCTL                0.13       0.04        0.04
# number.of.miRNAs:venom.familyEXO                0.05       0.08       -0.10
# number.of.miRNAs:venom.familymyotoxin          -0.17       0.95       -2.35
# number.of.miRNAs:venom.familyohanin             0.15       2.38       -5.78
# number.of.miRNAs:venom.familyPLA2               0.13       0.04        0.06
# number.of.miRNAs:venom.familySVMP               0.06       0.03       -0.01
# number.of.miRNAs:venom.familySVSP               0.08       0.04       -0.01
# number.of.miRNAs:venom.familyVEGF               0.14       0.47       -0.86
# number.of.miRNAs:venom.familyvQC                0.07       0.07       -0.04
# hu_number.of.miRNAs                            -0.24       0.31       -0.88
# hu_venom.familyCTL                         122406.57   74782.89    18336.56
# hu_venom.familyEXO                        -625131.54 1376697.99 -3943490.75
# hu_venom.familymyotoxin                   -998205.80 1021387.29 -3306485.50
# hu_venom.familyohanin                    -1261319.03 1488632.26 -5530776.25
# hu_venom.familyPLA2                       -190716.57  231777.77  -699602.68
# hu_venom.familySVMP                        -46488.04   48121.91  -164203.92
# hu_venom.familySVSP                       -191301.34  110258.53  -435284.40
# hu_venom.familyVEGF                      -1114956.92 1293487.28 -3634553.25
# hu_venom.familyvQC                       -2108851.00 2605695.98 -8659412.50
# hu_number.of.miRNAs:venom.familyCTL         -4761.91    2990.42   -13528.90
# hu_number.of.miRNAs:venom.familyEXO        -15227.22   42927.57  -103564.67
# hu_number.of.miRNAs:venom.familymyotoxin    24306.09   89917.36  -106214.52
# hu_number.of.miRNAs:venom.familyohanin      53460.29  336344.28  -488199.10
# hu_number.of.miRNAs:venom.familyPLA2       -87635.38   71559.70  -222534.92
# hu_number.of.miRNAs:venom.familySVMP       -11098.41   15388.70   -67451.04
# hu_number.of.miRNAs:venom.familySVSP       -40823.65   48356.03  -178526.02
# hu_number.of.miRNAs:venom.familyVEGF        -3831.28   56575.09  -149773.90
# hu_number.of.miRNAs:venom.familyvQC          5172.80  278527.92  -745715.82
#                                            u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                      3.28 2.39        7       15
# hu_Intercept                              661158.95 1.63       10       19
# number.of.miRNAs                              -0.01 2.75        7       12
# venom.familyCTL                               -1.41 1.47       12       26
# venom.familyEXO                                6.75 1.38       14       42
# venom.familymyotoxin                          24.45 1.25       19       20
# venom.familyohanin                            23.18 1.62       10       30
# venom.familyPLA2                               0.45 2.10        8       12
# venom.familySVMP                               0.71 2.68        7       15
# venom.familySVSP                               0.96 2.46        7       13
# venom.familyVEGF                              17.69 1.49       11       23
# venom.familyvQC                                1.22 1.48       12       32
# number.of.miRNAs:venom.familyCTL               0.19 1.97        8       29
# number.of.miRNAs:venom.familyEXO               0.19 1.31       17       45
# number.of.miRNAs:venom.familymyotoxin          1.44 1.25       19       21
# number.of.miRNAs:venom.familyohanin            4.39 1.63       10       30
# number.of.miRNAs:venom.familyPLA2              0.22 1.44       12       24
# number.of.miRNAs:venom.familySVMP              0.12 2.72        7       12
# number.of.miRNAs:venom.familySVSP              0.16 2.18        8       14
# number.of.miRNAs:venom.familyVEGF              0.92 1.49       11       24
# number.of.miRNAs:venom.familyvQC               0.21 1.44       12       34
# hu_number.of.miRNAs                            0.30 1.12       41      160
# hu_venom.familyCTL                        347016.20 2.33        7       11
# hu_venom.familyEXO                       1459786.50 1.98        8       49
# hu_venom.familymyotoxin                   117323.32 2.68        7       19
# hu_venom.familyohanin                     325843.92 2.38        8       12
# hu_venom.familyPLA2                       136326.35 1.75        9       19
# hu_venom.familySVMP                         9385.33 1.62       10       42
# hu_venom.familySVSP                          565.44 1.41       13       26
# hu_venom.familyVEGF                      1649238.50 1.88        9       16
# hu_venom.familyvQC                       1009892.25 2.86        7       11
# hu_number.of.miRNAs:venom.familyCTL         -698.65 2.26        8       11
# hu_number.of.miRNAs:venom.familyEXO        71705.67 1.76        9       21
# hu_number.of.miRNAs:venom.familymyotoxin  212686.77 1.84        9       20
# hu_number.of.miRNAs:venom.familyohanin    839618.00 1.98        8       11
# hu_number.of.miRNAs:venom.familyPLA2       46025.76 1.51       12       21
# hu_number.of.miRNAs:venom.familySVMP         582.46 1.76        9       11
# hu_number.of.miRNAs:venom.familySVSP       21383.70 2.03        8       11
# hu_number.of.miRNAs:venom.familyVEGF       90435.56 1.73        9       20
# hu_number.of.miRNAs:venom.familyvQC       623476.20 3.14        7       11

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.38      0.02     0.35     0.43 1.12       39      126

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# Summarize model5
model6_summary <- summary(venom_model6)
model6_summary
#  Family: hurdle_lognormal 
#   Links: mu = identity; sigma = identity; hu = logit 
# Formula: mRNA.ntd ~ number.of.miRNAs * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix)) 
#          hu ~ number.of.miRNAs * venom.family + (1 | gr(sample.id, cov = phylo_cov_matrix))
#    Data: filt_miRNA_num_clr_df (Number of observations: 116) 
#   Draws: 6 chains, each with iter = 5e+05; warmup = 125000; thin = 500;
#          total post-warmup draws = 4500

# Multilevel Hyperparameters:
# ~sample.id (Number of levels: 6) 
#                     Estimate   Est.Error   l-95% CI    u-95% CI Rhat Bulk_ESS
# sd(Intercept)           0.09        0.06       0.01        0.25 1.02      248
# sd(hu_Intercept) 18555819.41 20375533.24 1169275.75 77533795.00 2.42        7
#                  Tail_ESS
# sd(Intercept)         323
# sd(hu_Intercept)       16

# Regression Coefficients:
#                                            Estimate   Est.Error      l-95% CI
# Intercept                                     -0.61        0.36         -1.28
# hu_Intercept                            28365106.46 28428049.87    1954431.75
# number.of.miRNAs                               0.48        0.09          0.32
# venom.familyCTL                               -0.57       10.08        -20.15
# venom.familyEXO                                3.10        1.20          0.57
# venom.familyohanin                             1.57        6.69        -11.58
# venom.familyPLA2                               3.42        0.54          2.40
# venom.familySVMP                               3.27        0.36          2.57
# venom.familySVSP                               1.17        6.89        -12.42
# venom.familyVEGF                              -0.67        9.55        -18.23
# venom.familyvQC                                0.41        8.97        -16.88
# number.of.miRNAs:venom.familyCTL              -0.24        1.01         -2.27
# number.of.miRNAs:venom.familyEXO              -0.54        0.22         -0.95
# number.of.miRNAs:venom.familyohanin            1.18        6.69        -11.92
# number.of.miRNAs:venom.familyPLA2             -0.50        0.25         -0.98
# number.of.miRNAs:venom.familySVMP             -0.52        0.09         -0.69
# number.of.miRNAs:venom.familySVSP              1.64        6.90        -11.82
# number.of.miRNAs:venom.familyVEGF              0.15        1.59         -2.88
# number.of.miRNAs:venom.familyvQC               0.91        4.49         -7.80
# hu_number.of.miRNAs                     -2648769.40  3083045.22  -12116877.50
# hu_venom.familyCTL                     -22148165.74 61513551.19 -178569175.00
# hu_venom.familyEXO                     -43529863.92 75271542.26 -230034150.00
# hu_venom.familyohanin                  -18568304.98 31757219.99 -101432900.00
# hu_venom.familyPLA2                     -9432920.76 22592260.86  -72005487.50
# hu_venom.familySVMP                    -20006898.27 30535862.43 -118211900.00
# hu_venom.familySVSP                    -10550186.27 25576877.03  -72110222.50
# hu_venom.familyVEGF                    -16397472.70 28644919.14  -78851540.00
# hu_venom.familyvQC                      -3059238.26 15451094.10  -31877427.50
# hu_number.of.miRNAs:venom.familyCTL       202137.22  4422959.56   -7986276.00
# hu_number.of.miRNAs:venom.familyEXO      2226611.78 10192642.04   -9862944.25
# hu_number.of.miRNAs:venom.familyohanin  -6336723.90 27728430.72  -59753250.00
# hu_number.of.miRNAs:venom.familyPLA2   -19482818.77 25327939.26  -91931897.50
# hu_number.of.miRNAs:venom.familySVMP      763069.29  2901349.08   -2996892.25
# hu_number.of.miRNAs:venom.familySVSP    -7742578.85 30143207.67  -85971962.50
# hu_number.of.miRNAs:venom.familyVEGF    -2373293.34  5758590.72  -18957160.00
# hu_number.of.miRNAs:venom.familyvQC    -12085383.65 17188117.24  -64007807.50
#                                            u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept                                      0.08 1.05       97      276
# hu_Intercept                           104107200.00 3.11        7       12
# number.of.miRNAs                               0.65 1.06       79      198
# venom.familyCTL                               19.81 1.04      186      633
# venom.familyEXO                                5.34 1.02      339      449
# venom.familyohanin                            14.71 1.01      372      601
# venom.familyPLA2                               4.50 1.05      101      279
# venom.familySVMP                               3.98 1.05       85      236
# venom.familySVSP                              14.69 1.02      266      674
# venom.familyVEGF                              17.60 1.02      244      810
# venom.familyvQC                               17.77 1.01      284      451
# number.of.miRNAs:venom.familyCTL               1.73 1.04      180      587
# number.of.miRNAs:venom.familyEXO              -0.06 1.02      323      432
# number.of.miRNAs:venom.familyohanin           14.28 1.01      372      676
# number.of.miRNAs:venom.familyPLA2             -0.02 1.03      208      601
# number.of.miRNAs:venom.familySVMP             -0.36 1.06       80      186
# number.of.miRNAs:venom.familySVSP             15.22 1.02      261      666
# number.of.miRNAs:venom.familyVEGF              3.06 1.02      244      842
# number.of.miRNAs:venom.familyvQC               9.58 1.01      283      455
# hu_number.of.miRNAs                       -86390.41 2.83        7       12
# hu_venom.familyCTL                      70284595.00 2.02        8       11
# hu_venom.familyEXO                      36195202.50 1.76        9       15
# hu_venom.familyohanin                   14268725.00 2.17        8       15
# hu_venom.familyPLA2                     17609512.50 2.52        7       14
# hu_venom.familySVMP                      -125390.15 3.11        7       14
# hu_venom.familySVSP                     24478365.00 2.50        7       26
# hu_venom.familyVEGF                     29543732.50 1.47       12       32
# hu_venom.familyvQC                      39162337.50 1.77        9       11
# hu_number.of.miRNAs:venom.familyCTL      9649594.25 2.21        8       21
# hu_number.of.miRNAs:venom.familyEXO     31094335.00 1.74        9       28
# hu_number.of.miRNAs:venom.familyohanin  58666572.50 1.96        8       13
# hu_number.of.miRNAs:venom.familyPLA2     4426435.50 2.24        8       12
# hu_number.of.miRNAs:venom.familySVMP     8133741.25 1.80       11       39
# hu_number.of.miRNAs:venom.familySVSP    41837835.00 2.81        7       18
# hu_number.of.miRNAs:venom.familyVEGF     6203747.00 1.84        9       15
# hu_number.of.miRNAs:venom.familyvQC      4787004.75 2.03        8       15

# Further Distributional Parameters:
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma     0.26      0.02     0.22     0.30 1.01      474     1072

# Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).


# Plot and examine venom model 6 ----
## Format the input data ----

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
any(is.na(filt_miRNA_num_clr_df))
glimpse(filt_miRNA_num_clr_df)


## Get tidy dat for the model ----
tidy <- broom.mixed::tidy(venom_model6)
tidy


## Get epreds ----
model6_df <- filt_miRNA_num_clr_df %>%
  # Get epreds
  tidybayes::add_epred_draws(
    venom_model6,
    .draw = TRUE,
    .chain = TRUE,
    .iteration = TRUE,
    dpar = 'hu'
  )
glimpse(model6_df)

## Get probability of expression for the data ----
expression_prob_df <- filt_miRNA_num_clr_df %>%
  # Create a expresion probability column based on whether the mRNA is expressed
  mutate(
    mRNA.expressed = if_else(mRNA.counts > 1000, 1, 0),
    protein.expressed = if_else(intensity > 0, 1, 0)
  )
glimpse(expression_prob_df)


## Plot conditional effects ----

# Create a plot for the conditional effects of the entire model
model6_conditional_effects <- conditional_effects(venom_model6)
model6_conditional_effects

# Create a plot for the conditional effects for the hu part of the model
model6_hu <- conditional_effects(venom_model6, dpar = 'hu')
model6_hu

# Creat a plot for the conditional effects for the mu part of the model
model6_mu <- conditional_effects(venom_model6, dpar = 'mu')
model6_mu

## Plot effects with ggplot ----

# Create a plot with tidybayes output for the model
model6_hu_plot <- ggplot(
  model6_df,
  aes(
    x = number.of.miRNAs,
    y = hu
  )
) +
  ggdist::stat_lineribbon(
    .width = c(0.99, 0.95, 0.8, 0.5),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(
    title = 'Probability of mRNA expression',
    x = 'Number of targeting miRNAs',
    y = 'Probability of mRNA expression'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
model6_hu_plot
ggsave('Figures/Expression_Plots/BRMS/miRNA_numbers_vs_mRNA/hurdle_lognormal/venom_model6_hu_plot_2025.02.18.png', width = 10, height = 10, plot = model6_hu_plot, create.dir = TRUE)

### Try plotting with a regular glm ----
# Create a plot with tidybayes output for the model
logistic_regression_plot_mRNA <- ggplot(
  expression_prob_df,
  aes(
    x = number.of.miRNAs,
    y = mRNA.expressed
  )
) +
  geom_smooth(
    formula = y ~ x,
    method = 'glm',
    method.args = list(family = 'binomial'),
    se = TRUE,
    color = 'blue',
    fill = 'grey'
  ) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(
    title = 'Probability of mRNA expression',
    x = 'Number of targeting miRNAs',
    y = 'Probability of mRNA expression'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
logistic_regression_plot_mRNA
ggsave('Figures/Expression_Plots/GLM/miRNA_numbers_vs_mRNA/Bernoulli/miRNA_numbers_vs_mRNA_expression_2025.02.18.png', width = 10, height = 10, plot = logistic_regression_plot_mRNA, create.dir = TRUE)


# Create a plot with tidybayes output for the model
logistic_regression_plot_protein <- ggplot(
  expression_prob_df,
  aes(
    x = number.of.miRNAs,
    y = protein.expressed
  )
) +
  geom_smooth(
    formula = y ~ x,
    method = 'glm',
    method.args = list(family = 'binomial'),
    se = TRUE,
    color = 'blue',
    fill = 'grey'
  ) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(
    title = 'Probability of protein expression',
    x = 'Number of targeting miRNAs',
    y = 'Probability of protein expression'
  ) +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold', margin = margin(b = 5, t = 5), size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.position = 'bottom',
    legend.title.position = 'top'
  )
logistic_regression_plot_protein
ggsave('Figures/Expression_Plots/GLM/miRNA_numbers_vs_Protein/Bernoulli/miRNA_numbers_vs_Protein_expression_2025.02.18.png', width = 10, height = 10, plot = logistic_regression_plot_protein, create.dir = TRUE)
