# Last Edited: 2025/01/23

# Set up and Read Data ----

# Load in packages
library(compositions)
library(cowplot)
library(tidyverse)
library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(readxl)
library(ggpmisc)
library(ggplot2)
library(pheatmap)
# Load libraries that are important for the tree model
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
library(arrow)
# Source my functions
source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')

# Create variable for the fused dataset.
# If this is in Dropbox the file needs to be on the machine
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read both in as data frames
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)
glimpse(miRNA_mRNA_protein_df)

# Create a list of columns to remove
removed_columns <- c(
  'miRNA.target.chrom', 'miRNA.target.start', 'miRNA.target.end', 'miRNA.target.length', 
  'strand', 'best.miRNA.ortholog', 'miRNA.sequence.chrom', 'miRNA.start', 'miRNA.end', 'miRNA.strandedness', 
  'miRNA.length', 'miRNA.name.probability', 'blast.percent.identity', 'E.value', 'bit.score', 
  'miRNA.sequence', 'hairpin.sequence'
)

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    !sample.id == 'CV1082_viridis'
  ) %>% # This should get rid of all of the weirdly annotated genes
  dplyr::select(
    -all_of(removed_columns) # Remove columns to save memory
  ) %>% 
  distinct()
rm(miRNA_mRNA_protein_df)

# Tree file for the 12 snakes samples
tree_file <- 'Data/Tree/11Snake_Tree_ViridisPolytomy.phy'


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
five_prime_color <- '#0072b2'
cds_color <- '#d55e00'
viridis_color <- '#2D8A5C'
lutosus_color <- '#8C4720'
concolor_color <- '#E0A229'



# Set up and format tree ----

# Read tree file in
snake_tree <- ape::read.tree(tree_file)

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
  "CV0857_viridis_North_M" = 'C.viridis_CV0857', 
  "CV1081_viridis_Mid_M" = 'C.viridis_CV1081', 
  "CV1087_viridis_North_F" = 'C.viridis_CV1087', 
  "CV1086_viridis_South_M" = 'C.viridis_CV1086', 
  "CV0985_concolor_Other_F" = 'C.concolor_CV0985', 
  "CV0987_lutosus_Other_F" = 'C.lutosus_CV0987'
)


# Rename tips correctly using a match
for (i in seq_along(tip_names)) {
  old_name <- names(tip_names)[i]
  new_name <- tip_names[i]
  snake_tree$tip.label[snake_tree$tip.label == old_name] <- new_name
}

# Check tree
snake_tree$tip.label
# Replace underscores with newlines
snake_tree$tip.label <- gsub("_", "\n", snake_tree$tip.label)

# Make tree plot
tree_plot <- ggtree(snake_tree, layout = 'rectangular', size = 1) +
  geom_tiplab(size = 3, angle = 0) +
  geom_text2(
    aes(label = round(branch.length, 2), subset = !isTip),
    vjust = -0.5,
    hjust = 1.25, size = 3
  ) +
  # theme_tree2() +
  xlim(-0.1, 1.1)  # Expand x-axis range by 20%
tree_plot
ggsave('Figures/Trees/Crotalus_Tree_2025.01.23.pdf', device = 'pdf', plot = tree_plot, height = 4, width = 6, dpi = 900, create.dir = T)



# Venom Transcriptome Profile Pie Charts ----

# Remove the following characters from the sample.id. column
characters_to_remove <- c(
  '_viridis',
  '_viridis',
  '_viridis',
  '_viridis',
  '_concolor',
  '_lutosus'
)

# Read mRNA data
mRNA_data <- 'Data/mRNA/mRNA_raw_counts-vst-rlog-norm_counts.2025.01.22.parquet'

mRNA_df <- read_parquet(file = mRNA_data) %>% 
  mutate(sample.id = str_replace(sample.id, paste(characters_to_remove, collapse = "|"), '')) %>% 
  # Remove unnecessary expression info
  select(
    sample.id, genes, venom.family, mRNA.vst
  ) %>% 
  pivot_wider(
    names_from = sample.id,
    values_from = mRNA.vst
  ) %>% 
  filter(str_starts(genes, 'Venom_')) %>% 
  filter(!str_starts(genes, 'Venom_ADAM')) %>%
  # Change the venom.family column so that only 
  mutate(
    venom.family = case_when(
      grepl('SVMP', genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('SVSP', genes) ~ 'SVSP',
      grepl('PLA2', genes) ~ 'PLA2',
      grepl('CRISP', genes) ~ 'CRISP',
      grepl('LAAO', genes) ~ 'LAAO',
      grepl('myotoxin', genes) ~ 'myotoxin',
      grepl('BPP', genes) ~ 'BPP',
      # grepl('CTL', genes) ~ 'CTL',
      # grepl('EXO', genes) ~ 'EXO',
      # grepl('VEGF', genes) ~ 'VEGF',
      # grepl('ohanin', genes) ~ 'ohanin',
      # grepl('vQC', genes) ~ 'vQC',
      TRUE ~ 'others'
    )
  ) %>% 
  mutate(
    color = case_when(
      grepl('SVMP', genes) ~ SVMP_color,
      grepl('SVSP', genes) ~ SVSP_color,
      grepl('PLA2', genes) ~ PLA2_color,
      grepl('CRISP', genes) ~ CRISP_color,
      grepl('LAAO', genes) ~ LAAO_color,
      grepl('myotoxin', genes) ~ myotoxin_color,
      grepl('BPP', genes) ~ BPP_color,
      # grepl('VEGF', genes) ~ VEGF_color,
      # grepl('ohanin', genes) ~ ohanin_color,
      # grepl('vQC', genes) ~ vQC_color,
      # grepl('CTL', genes) ~ CTL_color,
      # grepl('EXO', genes) ~ EXO_color,
      TRUE ~ other_color
    )
  )

## Limit the data to CV0857 ----
mRNA_CV0857_df <- mRNA_df %>% 
  dplyr::select(venom.family, CV0857, color) %>% 
  group_by(venom.family) %>% 
  mutate(expression = mean(CV0857)) %>% 
  dplyr::select(-CV0857) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((expression / sum(expression)) * 100, 2))  # round to 1 decimal place


# Create a vector for colors
colors <- setNames(mRNA_CV0857_df$color, mRNA_CV0857_df$venom.family)


# Create pie chart for the mRNA data for CV0857
venom_mRNA_CV0857_pie_chart <- mRNA_CV0857_df %>% 
  ggplot(aes(x = "", y = expression, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    # title = expression(paste('Venom Transcriptome in Northern ', italic('C. viridis'), ' (CV0857)'))
    title = expression(paste(italic('C. viridis'), ' (CV0857)'))
  )
venom_mRNA_CV0857_pie_chart
ggsave('Figures/Pie_Charts/Transcriptome/Transcripome_viridis_CV0857_2025.01.23.pdf', venom_mRNA_CV0857_pie_chart, create.dir = T)


## Limit the data to CV1087 ----
mRNA_CV1087_df <- mRNA_df %>% 
  dplyr::select(venom.family, CV1087, color) %>% 
  group_by(venom.family) %>% 
  mutate(expression = mean(CV1087)) %>% 
  dplyr::select(-CV1087) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((expression / sum(expression)) * 100, 2))  # round to 1 decimal place

# Create a vector for colors
colors <- setNames(mRNA_CV1087_df$color, mRNA_CV1087_df$venom.family)

# Create pie chart for the mRNA data for CV1087
venom_mRNA_CV1087_pie_chart <- mRNA_CV1087_df %>% 
  ggplot(aes(x = "", y = expression, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    # title = expression(paste('Venom Transcriptome in Northern ', italic('C. viridis'), ' (CV1087)'))
    title = expression(paste(italic('C. viridis'), ' (CV1087)'))
  )
venom_mRNA_CV1087_pie_chart
ggsave('Figures/Pie_Charts/Transcriptome/Transcripome_viridis_CV1087_2025.01.23.pdf', venom_mRNA_CV1087_pie_chart, create.dir = T)



## Limit the data to CV1081 ----
mRNA_CV1081_df <- mRNA_df %>% 
  dplyr::select(venom.family, CV1081, color) %>% 
  group_by(venom.family) %>% 
  mutate(expression = mean(CV1081)) %>% 
  dplyr::select(-CV1081) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((expression / sum(expression)) * 100, 2))  # round to 1 decimal place

# Create a vector for colors
colors <- setNames(mRNA_CV1081_df$color, mRNA_CV1081_df$venom.family)

# Create pie chart for the mRNA data for CV1081
venom_mRNA_CV1081_pie_chart <- mRNA_CV1081_df %>% 
  ggplot(aes(x = "", y = expression, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    title = expression(paste(italic('C. viridis'), ' (CV1081)'))
  )
venom_mRNA_CV1081_pie_chart
ggsave('Figures/Pie_Charts/Transcriptome/Transcripome_viridis_CV1081_2025.01.23.pdf', venom_mRNA_CV1081_pie_chart, create.dir = T)


## Limit the data to CV1086 ----
mRNA_CV1086_df <- mRNA_df %>% 
  dplyr::select(venom.family, CV1086, color) %>% 
  group_by(venom.family) %>% 
  mutate(expression = mean(CV1086)) %>% 
  dplyr::select(-CV1086) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((expression / sum(expression)) * 100, 2))  # round to 1 decimal place

# Create a vector for colors
colors <- setNames(mRNA_CV1086_df$color, mRNA_CV1086_df$venom.family)

# Create pie chart for the mRNA data for CV1086
venom_mRNA_CV1086_pie_chart <- mRNA_CV1086_df %>% 
  ggplot(aes(x = "", y = expression, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    title = expression(paste(italic('C. viridis'), ' (CV1086)'))
  )
venom_mRNA_CV1086_pie_chart
ggsave('Figures/Pie_Charts/Transcriptome/Transcripome_viridis_CV1086_2025.01.23.pdf', venom_mRNA_CV1086_pie_chart, create.dir = T)


## Limit the data to CV0987 ----
mRNA_CV0987_df <- mRNA_df %>% 
  dplyr::select(venom.family, CV0987, color) %>% 
  group_by(venom.family) %>% 
  mutate(expression = mean(CV0987)) %>% 
  dplyr::select(-CV0987) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((expression / sum(expression)) * 100, 2))  # round to 1 decimal place

# Create a vector for colors
colors <- setNames(mRNA_CV0987_df$color, mRNA_CV0987_df$venom.family)

# Create pie chart for the mRNA data for CV0987
venom_mRNA_CV0987_pie_chart <- mRNA_CV0987_df %>% 
  ggplot(aes(x = "", y = expression, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    title = expression(paste(italic('C. lutosus'), ' (CV0987)'))
  )
venom_mRNA_CV0987_pie_chart
ggsave('Figures/Pie_Charts/Transcriptome/Transcripome_lutosus_CV0987_2025.01.23.pdf', venom_mRNA_CV0987_pie_chart, create.dir = T)


## Limit the data to CV0985 ----
mRNA_CV0985_df <- mRNA_df %>% 
  dplyr::select(venom.family, CV0985, color) %>% 
  group_by(venom.family) %>% 
  mutate(expression = mean(CV0985)) %>% 
  dplyr::select(-CV0985) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((expression / sum(expression)) * 100, 2))  # round to 1 decimal place

# Create a vector for colors
colors <- setNames(mRNA_CV0985_df$color, mRNA_CV0985_df$venom.family)


# Create pie chart for the mRNA data for CV0985
venom_mRNA_CV0985_pie_chart <- mRNA_CV0985_df %>% 
  ggplot(aes(x = "", y = expression, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    title = expression(paste(italic('C. concolor'), ' (CV0985)'))
  )
venom_mRNA_CV0985_pie_chart
ggsave('Figures/Pie_Charts/Transcriptome/Transcripome_concolor_CV0985_2025.01.23.pdf', venom_mRNA_CV0985_pie_chart, create.dir = T)


## Remove legends for the following pie charts ----
venom_mRNA_CV0857_pie_chart2 <- venom_mRNA_CV0857_pie_chart +
  theme(legend.position = 'none')
venom_mRNA_CV1087_pie_chart2 <- venom_mRNA_CV1087_pie_chart +
  theme(legend.position = 'none')
venom_mRNA_CV1081_pie_chart2 <- venom_mRNA_CV1081_pie_chart +
  theme(legend.position = 'none')
venom_mRNA_CV1086_pie_chart2 <- venom_mRNA_CV1086_pie_chart +
  theme(legend.position = 'none')
venom_mRNA_CV0985_pie_chart2 <- venom_mRNA_CV0985_pie_chart +
  theme(legend.position = 'none')
venom_mRNA_CV0987_pie_chart2 <- venom_mRNA_CV0987_pie_chart +
  theme(legend.position = 'none')



# Create a set of transcriptome figures
transcriptome_six_samples <- plot_grid(
  venom_mRNA_CV0857_pie_chart2, 
  venom_mRNA_CV1087_pie_chart2, 
  venom_mRNA_CV1081_pie_chart2, 
  venom_mRNA_CV1086_pie_chart2,
  venom_mRNA_CV0985_pie_chart2,
  venom_mRNA_CV0987_pie_chart,
  ncol = 6,
  align = 'v'
)
transcriptome_six_samples
ggsave('Figures/Pie_Charts/Transcriptome/All_Transcriptomes_2025.01.23.pdf', plot = transcriptome_six_samples, width = 30, height = 8, dpi = 900, create.dir = T) # Make sure to change the size of the plots window when saving this

# Create a smaller set
transcriptome_four_samples <- plot_grid(
  venom_mRNA_CV1087_pie_chart2,
  venom_mRNA_CV0857_pie_chart2,
  venom_mRNA_CV0985_pie_chart2,
  venom_mRNA_CV0987_pie_chart,
  ncol = 4,
  align = 'v'
)
transcriptome_four_samples
ggsave('Figures/Pie_Charts/Transcriptome/Four_Transcriptomes_2025.01.23.pdf', plot = transcriptome_four_samples, width = 24, height = 8, dpi = 900, create.dir = T) # Make sure to change the size of the plots window when saving this

# Create a smaller set
transcriptome_three_samples <- plot_grid(
  venom_mRNA_CV1087_pie_chart2,
  venom_mRNA_CV0985_pie_chart2,
  venom_mRNA_CV0987_pie_chart,
  ncol = 3,
  align = 'v'
)
transcriptome_three_samples
# Create a title
transcriptome_three_samples <- ggdraw() +
  draw_label("a. Transcriptome profiles of three species", fontface = 'bold', size = 14, x = 0.5, y = 0.95, hjust = 0.5) +
  draw_plot(transcriptome_three_samples, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
transcriptome_three_samples
ggsave('Figures/Pie_Charts/Transcriptome/Three_Transcriptomes_2025.01.23.pdf', plot = transcriptome_three_samples, width = 24, height = 8, dpi = 900, create.dir = T) # Make sure to change the size of the plots window when saving this



# Venom Proteome Profile Pie Charts ----

# Remove the following characters from the sample.id. column
characters_to_remove <- c(
  '_viridis',
  '_viridis',
  '_viridis',
  '_viridis',
  '_concolor',
  '_lutosus'
)

# Set and read the data
protein_data <- 'Data/Protein/protein_intensity_data.2025.01.22.parquet'
protein_df <- read_parquet(file = protein_data) %>% 
  mutate(sample.id = str_replace(sample.id, paste(characters_to_remove, collapse = "|"), '')) %>% 
  filter(in.library == 'Yes') %>% # Add this filtering step so as not to include proteins that were not in the library sent to Anthony
  dplyr::select(intensity, sample.id, genes) %>% 
  filter(str_starts(genes, 'Venom_')) %>% 
  filter(!str_detect(genes, 'ADAM')) %>% 
  pivot_wider(
    names_from = sample.id,
    values_from = intensity
  ) %>% 
  mutate(
    venom.family = case_when(
      grepl('SVMP', genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('SVSP', genes) ~ 'SVSP',
      grepl('PLA2', genes) ~ 'PLA2',
      grepl('CRISP', genes) ~ 'CRISP',
      grepl('LAAO', genes) ~ 'LAAO',
      grepl('myotoxin', genes) ~ 'myotoxin',
      grepl('BPP', genes) ~ 'BPP',
      # grepl('CTL', genes) ~ 'CTL',
      # grepl('EXO', genes) ~ 'EXO',
      # grepl('ADAM', genes) ~ 'ADAM',
      # grepl('VEGF', genes) ~ 'VEGF',
      # grepl('ohanin', genes) ~ 'ohanin',
      # grepl('vQC', genes) ~ 'vQC',
      TRUE ~ 'others'
    )
  ) %>% 
  mutate(
    color = case_when(
      grepl('SVMP', genes) ~ SVMP_color,
      grepl('SVSP', genes) ~ SVSP_color,
      grepl('PLA2', genes) ~ PLA2_color,
      grepl('CRISP', genes) ~ CRISP_color,
      grepl('LAAO', genes) ~ LAAO_color,
      grepl('myotoxin', genes) ~ myotoxin_color,
      grepl('BPP', genes) ~ BPP_color,
      # grepl('VEGF', genes) ~ VEGF_color,
      # grepl('ohanin', genes) ~ ohanin_color,
      # grepl('vQC', genes) ~ vQC_color,
      # grepl('CTL', genes) ~ CTL_color,
      # grepl('EXO', genes) ~ EXO_color,
      # grepl('ADAM', genes) ~ ADAM_color,
      TRUE ~ other_color
    )
  )

## Create a new data frame for only CV0987 ----
prot_CV0987_df <- protein_df %>% 
  dplyr::select(venom.family, CV0987, color) %>% 
  group_by(venom.family) %>% 
  mutate(total.intensity = mean(CV0987)) %>% 
  dplyr::select(-CV0987) %>% 
  # filter(total.intensity > 0) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((total.intensity / sum(total.intensity)) * 100, 2))  # round to 2 decimal place

# Create a vector for colors
colors <- setNames(prot_CV0987_df$color, prot_CV0987_df$venom.family)


# Create plot
venom_protein_CV0987_pie_chart <- prot_CV0987_df %>% 
  ggplot(aes(x = "", y = total.intensity, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    # title = expression(paste('Venom Proteome in ', italic('C. lutosus'), ' (CV0987)'))
    title = expression(paste(italic('C. lutosus'), ' (CV0987)'))
  )
venom_protein_CV0987_pie_chart
ggsave('Figures/Pie_Charts/Proteome/Proteome_lutosus_CV0987_2025.01.23.pdf', venom_protein_CV0987_pie_chart, create.dir = T)


## Create a new data frame for only CV0985 ----
prot_CV0985_df <- protein_df %>% 
  dplyr::select(venom.family, CV0985, color) %>% 
  group_by(venom.family) %>% 
  mutate(total.intensity = mean(CV0985)) %>% 
  dplyr::select(-CV0985) %>% 
  # filter(total.intensity > 0) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((total.intensity / sum(total.intensity)) * 100, 2))  # round to 2 decimal place

# Create a vector for colors
colors <- setNames(prot_CV0985_df$color, prot_CV0985_df$venom.family)


# Create plot
venom_protein_CV0985_pie_chart <- prot_CV0985_df %>% 
  ggplot(aes(x = "", y = total.intensity, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    # title = expression(paste('Venom Proteome in ', italic('C. concolor'), ' (CV0985)'))
    title = expression(paste(italic('C. concolor'), ' (CV0985)'))
  )
venom_protein_CV0985_pie_chart
ggsave('Figures/Pie_Charts/Proteome/Proteome_concolor_CV0985_2025.01.23.pdf', venom_protein_CV0985_pie_chart, create.dir = T)

## Create a new data frame for only CV1087 ----
prot_CV1087_df <- protein_df %>% 
  dplyr::select(venom.family, CV1087, color) %>% 
  group_by(venom.family) %>% 
  mutate(total.intensity = mean(CV1087)) %>% 
  dplyr::select(-CV1087) %>% 
  # filter(total.intensity > 0) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((total.intensity / sum(total.intensity)) * 100, 2))  # round to 2 decimal place

# Create a vector for colors
colors <- setNames(prot_CV1087_df$color, prot_CV1087_df$venom.family)


# Create plot
venom_protein_CV1087_pie_chart <- prot_CV1087_df %>% 
  ggplot(aes(x = "", y = total.intensity, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    # title = expression(paste('Venom Proteome in Northern ', italic('C. viridis'), ' (CV1087)'))
    title = expression(paste(italic('C. viridis'), ' (CV1087)'))
  )
venom_protein_CV1087_pie_chart
ggsave('Figures/Pie_Charts/Proteome/Proteome_viridis_CV1087_2025.01.23.pdf', venom_protein_CV1087_pie_chart, create.dir = T)



## Create a new data frame for only CV1086 ----
prot_CV1086_df <- protein_df %>% 
  dplyr::select(venom.family, CV1086, color) %>% 
  group_by(venom.family) %>% 
  mutate(total.intensity = mean(CV1086)) %>% 
  dplyr::select(-CV1086) %>% 
  # filter(total.intensity > 0) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((total.intensity / sum(total.intensity)) * 100, 2))  # round to 2 decimal place

# Create a vector for colors
colors <- setNames(prot_CV1086_df$color, prot_CV1086_df$venom.family)


# Create plot
venom_protein_CV1086_pie_chart <- prot_CV1086_df %>% 
  ggplot(aes(x = "", y = total.intensity, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    # title = expression(paste('Venom Proteome in Southern ', italic('C. viridis'), ' (CV1086)'))
    title = expression(paste(italic('C. viridis'), ' (CV1086)'))
  )
venom_protein_CV1086_pie_chart
ggsave('Figures/Pie_Charts/Proteome/Proteome_viridis_CV1086_2025.01.23.pdf', venom_protein_CV1086_pie_chart, create.dir = T)


## Create a new data frame for only CV1081 ----
prot_CV1081_df <- protein_df %>% 
  dplyr::select(venom.family, CV1081, color) %>% 
  group_by(venom.family) %>% 
  mutate(total.intensity = mean(CV1081)) %>% 
  dplyr::select(-CV1081) %>% 
  # filter(total.intensity > 0) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((total.intensity / sum(total.intensity)) * 100, 2))  # round to 2 decimal place

# Create a vector for colors
colors <- setNames(prot_CV1081_df$color, prot_CV1081_df$venom.family)


# Create plot
venom_protein_CV1081_pie_chart <- prot_CV1081_df %>% 
  ggplot(aes(x = "", y = total.intensity, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family',
    # title = expression(paste('Venom Proteome in Central ', italic('C. viridis'), ' (CV1081)'))
    title = expression(paste(italic('C. viridis'), ' (CV1081)'))
    
  )
venom_protein_CV1081_pie_chart
ggsave('Figures/Pie_Charts/Proteome/Proteome_viridis_CV1081_2025.01.23.pdf', venom_protein_CV1081_pie_chart, create.dir = T)


## Create a new data frame for only CV0857 ----
prot_CV0857_df <- protein_df %>% 
  dplyr::select(venom.family, CV0857, color) %>% 
  group_by(venom.family) %>% 
  mutate(total.intensity = mean(CV0857)) %>% 
  dplyr::select(-CV0857) %>% 
  # filter(total.intensity > 0) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(percentage = round((total.intensity / sum(total.intensity)) * 100, 2))  # round to 2 decimal place

# Create a vector for colors
colors <- setNames(prot_CV0857_df$color, prot_CV0857_df$venom.family)


# Create plot
venom_protein_CV0857_pie_chart <- prot_CV0857_df %>% 
  ggplot(aes(x = "", y = total.intensity, fill = venom.family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  ggtitle(expression(paste('Venom Proteome in Northern ', italic('C. viridis')))) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 15),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Venom Family', 
    # title = expression(paste('Venom Proteome in Northern ', italic('C. viridis'), ' (CV0857)'))
    title = expression(paste(italic('C. viridis'), ' (CV0857)'))
  )
venom_protein_CV0857_pie_chart
ggsave('Figures/Pie_Charts/Proteome/Proteome_viridis_CV0857_2025.01.23.pdf', venom_protein_CV0857_pie_chart, create.dir = T)

# Remove legends for the following pie charts
venom_protein_CV0857_pie_chart2 <- venom_protein_CV0857_pie_chart +
  theme(legend.position = 'none')
venom_protein_CV1087_pie_chart2 <- venom_protein_CV1087_pie_chart +
  theme(legend.position = 'none')
venom_protein_CV1081_pie_chart2 <- venom_protein_CV1081_pie_chart +
  theme(legend.position = 'none')
venom_protein_CV1086_pie_chart2 <- venom_protein_CV1086_pie_chart +
  theme(legend.position = 'none')
venom_protein_CV0985_pie_chart2 <- venom_protein_CV0985_pie_chart +
  theme(legend.position = 'none')
venom_protein_CV0987_pie_chart2 <- venom_protein_CV0987_pie_chart +
  theme(legend.position = 'none')




# Create a set of proteome figures
proteome_six_samples <- plot_grid(
  venom_protein_CV0857_pie_chart2, 
  venom_protein_CV1087_pie_chart2, 
  venom_protein_CV1081_pie_chart2, 
  venom_protein_CV1086_pie_chart2,
  venom_protein_CV0985_pie_chart2,
  venom_protein_CV0987_pie_chart,
  ncol = 6,
  align = 'v'
  )
proteome_six_samples
# # Add a title to the combined plot and adjust the title position
# proteome_six_samples_with_title <- ggdraw() + 
#   draw_label("Proteome Comparison of Six Samples", fontface = 'bold', size = 14, x = 0.5, y = 0.70, hjust = 0.5) + 
#   draw_plot(proteome_six_samples, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
# 
# proteome_six_samples_with_title
ggsave('Figures/Pie_Charts/Proteome/All_Proteomes_2025.01.23.pdf', plot = proteome_six_samples, width = 30, height = 8, dpi = 900, create.dir = T) # Make sure to change the size of the plots window when saving this

# Create a smaller set
proteome_four_samples <- plot_grid(
  venom_protein_CV1087_pie_chart2,
  venom_protein_CV0857_pie_chart2,
  venom_protein_CV0985_pie_chart2,
  venom_protein_CV0987_pie_chart,
  ncol = 4,
  align = 'v'
)
proteome_four_samples
ggsave('Figures/Pie_Charts/Proteome/Four_Proteomes_2025.01.23.pdf', plot = proteome_four_samples, width = 24, height = 8, dpi = 900, create.dir = T) # Make sure to change the size of the plots window when saving this

# Create a smaller set
proteome_three_samples <- plot_grid(
  venom_protein_CV1087_pie_chart2,
  venom_protein_CV0985_pie_chart2,
  venom_protein_CV0987_pie_chart,
  ncol = 3,
  align = 'v'
)
proteome_three_samples
# Create a title
proteome_three_samples <- ggdraw() +
  draw_label("b. Proteom profiles of three species", fontface = 'bold', size = 14, x = 0.5, y = 0.95, hjust = 0.5) +
  draw_plot(proteome_three_samples, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
proteome_three_samples
ggsave('Figures/Pie_Charts/Proteome/Three_Proteomes_2025.01.23.pdf', plot = proteome_three_samples, width = 24, height = 8, dpi = 900, create.dir = T) # Make sure to change the size of the plots window when saving this





# Figure a, Transcriptomes and Proteomes ----

## Change mRNA plot titles ----
# Create another figure for figure one with proteome and transcriptome put together for each sample
# Change the title for the following pie charts
venom_mRNA_CV0857_pie_chart3 <- venom_mRNA_CV0857_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV0857) transcriptome'))
  )

venom_mRNA_CV1087_pie_chart3 <- venom_mRNA_CV1087_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV1087) transcriptome'))
  )

venom_mRNA_CV1081_pie_chart3 <- venom_mRNA_CV1081_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV1081) transcriptome'))
  )

venom_mRNA_CV1086_pie_chart3 <- venom_mRNA_CV1086_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV1086) transcriptome'))
  )

venom_mRNA_CV0985_pie_chart3 <- venom_mRNA_CV0985_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. concolor'), ' (CV0985) transcriptome'))
  )

venom_mRNA_CV0987_pie_chart3 <- venom_mRNA_CV0987_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. lutosus'), ' (CV0987) transcriptome'))
  )

## Change Protein plot titles ----
# Remove legends for the following pie charts and change title
venom_protein_CV0857_pie_chart3 <- venom_protein_CV0857_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV0857) proteome'))
  )

venom_protein_CV1087_pie_chart3 <- venom_protein_CV1087_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV1087) proteome'))
  )

venom_protein_CV1081_pie_chart3 <- venom_protein_CV1081_pie_chart +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV1081) proteome'))
  )

venom_protein_CV1086_pie_chart3 <- venom_protein_CV1086_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. viridis'), ' (CV1086) proteome'))
  )

venom_protein_CV0985_pie_chart3 <- venom_protein_CV0985_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. concolor'), ' (CV0985) proteome'))
  )

venom_protein_CV0987_pie_chart3 <- venom_protein_CV0987_pie_chart +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 11)
  ) +
  labs(
    title = expression(paste(italic('C. lutosus'), ' (CV0987) proteome'))
  )

## Combine mRNA and Protein plot ----
# Create figure 1
figure_a <- plot_grid(
  # CV1086
  venom_mRNA_CV1086_pie_chart3,
  venom_protein_CV1086_pie_chart3,
  # CV1087
  venom_mRNA_CV1087_pie_chart3,
  venom_protein_CV1087_pie_chart3,
  # CV1081
  venom_mRNA_CV1081_pie_chart3,
  venom_protein_CV1081_pie_chart3,
  # CV0857
  venom_mRNA_CV0857_pie_chart3,
  venom_protein_CV0857_pie_chart3,
  # CV0987
  venom_mRNA_CV0987_pie_chart3,
  venom_protein_CV0987_pie_chart3,
  # CV0985
  venom_mRNA_CV0985_pie_chart3,
  venom_protein_CV0985_pie_chart3,
  ncol = 2,
  nrow = 6,
  align = 'v'
)
figure_a

# Add the tree to the plot
figure_a <- plot_grid(
  tree_plot,
  figure_a,
  align = 'v'
)
figure_a
# Add title
figure_a <- ggdraw() +
  draw_label("b. Transcriptome and proteome profiles of three species", fontface = 'bold', size = 14, x = 0.5, y = 0.95, hjust = 0.5) +
  draw_plot(figure_a, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
figure_a
ggsave("Figures/1_Main_Figures/Figure_a/Figure_a_2025.01.23.pdf", plot = figure_a, width = 10, height = 8, dpi = 2400, create.dir = T)




# Supplemental Table a ----

# Create df for transcriptome
mRNA_percentage_df <- left_join(
  # CV0857
  mRNA_CV0857_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV0857 = percentage),
  # CV0985
  mRNA_CV0985_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV0985 = percentage),
  by = 'venom.family'
) %>% 
  left_join(
  # CV1086
  mRNA_CV1086_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV1086 = percentage),
  by = 'venom.family'
  ) %>% 
  left_join(
  # CV1087
  mRNA_CV1087_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV1087 = percentage),
  by = 'venom.family'
  ) %>% 
  left_join(
  # CV1081
  mRNA_CV1081_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV1081 = percentage),
  by = 'venom.family'
  ) %>% 
  left_join(
  mRNA_CV0987_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV0987 = percentage),
  by = 'venom.family'
  ) %>% 
  pivot_longer(
    cols = contains('CV'),
    names_to = 'sample.id',
    values_to = 'Transcriptome %'
  )

# Create df for proteom
protein_percentage_df <- left_join(
  # CV0857
  prot_CV0857_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV0857 = percentage),
  # CV0985
  prot_CV0985_df %>% 
    select(venom.family, percentage) %>% 
    rename(CV0985 = percentage),
  by = 'venom.family'
) %>% 
  left_join(
    # CV1086
    prot_CV1086_df %>% 
      select(venom.family, percentage) %>% 
      rename(CV1086 = percentage),
    by = 'venom.family'
  ) %>% 
  left_join(
    # CV1087
    prot_CV1087_df %>% 
      select(venom.family, percentage) %>% 
      rename(CV1087 = percentage),
    by = 'venom.family'
  ) %>% 
  left_join(
    # CV1081
    prot_CV1081_df %>% 
      select(venom.family, percentage) %>% 
      rename(CV1081 = percentage),
    by = 'venom.family'
  ) %>% 
  left_join(
    prot_CV0987_df %>% 
      select(venom.family, percentage) %>% 
      rename(CV0987 = percentage),
    by = 'venom.family'
  ) %>% 
  pivot_longer(
    cols = contains('CV'),
    names_to = 'sample.id',
    values_to = 'Proteome %'
  )

# Fuse protein and mRNA together and create table a
table_a <- left_join(
  mRNA_percentage_df,
  protein_percentage_df,
  by = c('venom.family', 'sample.id')
) %>% 
  pivot_longer(
    cols = ends_with('%'),
    names_to = 'Type',
    values_to = 'percentage'
  ) %>% 
  pivot_wider(
    names_from = 'venom.family',
    values_from = 'percentage'
  )
write_csv3(table_a, file = 'Tables/2_Suplemental_Tables/Table_a/Table_a_2025.01.23.csv')



# mRNA Heat Maps ----

## All genes ----

# Create dataframe
mRNA_all_genes_df <- mi_df %>% 
  filter(in.library == 'Yes') %>% # Remove any proteins not in the library sent to Anthony
  dplyr::select(
    sample.id, genes, venom.family, mRNA.ntd
  ) %>% 
  distinct()
glimpse(mRNA_all_genes_df)
  

# Calculate variance in mRNA
mRNA_variance_df <- mRNA_all_genes_df %>%
  group_by(genes) %>%
  summarize(variance = var(mRNA.ntd)) %>% 
  ungroup()

# Create a heatmap of variances per gene
mRNA_variance_heatmap <- ggplot(mRNA_variance_df, aes(x = genes, y = 1, fill = log(variance + 1))) +
  geom_tile() +
  # scale_fill_viridis_c(option = 'viridis') +
  scale_fill_gradient(
    low = '#ededed',
    high = 'tomato3',
    breaks = seq(0, 4, by = 0.8)
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = 'Variance (Log)'
  ) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank()   # Remove y-axis ticks
  )
mRNA_variance_heatmap

# For these heat maps make sure to add a variance heat map over it 
# Create heat map
mRNA_all_heatmap <- ggplot(mRNA_all_genes_df, aes(y = sample.id, x = genes, fill = mRNA.ntd)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(
    y = 'Individuals', 
    x = 'Genes', 
    fill = 'mRNA expression' #, 
    # title = "mRNA expression Level for all Targeted genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margins
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 13)
  )
mRNA_all_heatmap
ggsave("Figures/Heat_Maps/mRNA/mRNA_all_genes_heatmap_2025.01.23.pdf", plot = mRNA_all_heatmap, width = 20, height = 4, dpi = 900, create.dir = T)

# Combined the main and variance heatmap
mRNA_combined_heatmap <- plot_grid(mRNA_variance_heatmap, mRNA_all_heatmap, ncol = 1, align = 'v', rel_heights = c(0.1, 1))
mRNA_combined_heatmap
ggsave("Figures/Heat_Maps/mRNA/mRNA_all_genes_with_variance_heatmap_2025.01.23.pdf", plot = mRNA_combined_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)


## Venom genes ----

# Create an order for the genes
venom_gene_order <- c(
  'BPP',
  'myotoxin',
  'ohanin',
  'PLA2B1', 'PLA2K', 'PLA2C1', 'PLA2A1',
  'SVSP1', 'SVSP2', 'SVSP3', 'SVSP4', 'SVSP10', 'SVSP5', 'SVSP6', 'SVSP11', 'SVSP7', 'SVSP8', 'SVSP9',
  'SVMP1', 'SVMP2', 'SVMP3', 'SVMP4', 'SVMP5', 'SVMP6', 'SVMP7', 'SVMP8', 'SVMP9', 'SVMP10', 'SVMP11', 'SVMP12',
  'CRISP1', 'CRISP2', 'CRISP3', 'CRISP4',
  'CTL1', 'CTL2', 'CTL3', 'CTL4', 'CTL5', 'CTL6',
  'EXO1', 'EXO2', 'EXO3',
  'LAAO1', 'LAAO2', 'LAAO3',
  'VEGF1', 'VEGF2',
  'vQC1', 'vQC2'
)

# Create an order for the samples
sample_order <- c(
  'CV0987_lutosus', 'CV0985_concolor', 'CV1086_viridis',
  'CV1081_viridis', 'CV0857_viridis', 'CV1087_viridis'
)

# Now for only venom
# Prepare data
all_venom_genes_mRNA_df <- mRNA_all_genes_df %>% 
  filter(
    str_starts(genes, 'Venom'),
    !str_detect(genes, 'ADAM')
  ) %>% 
  mutate(genes = str_remove(genes, '^Venom_')) %>% # Remove 'Venom_' from these genes
  distinct()

# Set the order of venom genes
all_venom_genes_mRNA_df$genes <- factor(all_venom_genes_mRNA_df$genes, levels = venom_gene_order)

# Set the order of the sample genes
all_venom_genes_mRNA_df$sample.id <- factor(all_venom_genes_mRNA_df$sample.id, levels = sample_order)

# Get the variance for the venom genes
venom_mRNA_variance_df <- all_venom_genes_mRNA_df %>% 
  group_by(genes) %>% 
  summarize(variance = var(mRNA.ntd)) %>% 
  ungroup()

# Create venom variance heatmap
venom_mRNA_variance_heatmap <- ggplot(venom_mRNA_variance_df, aes(x = genes, y = 1, fill = variance)) +
  geom_tile() +
  # scale_fill_viridis_c(option = 'viridis') +
  scale_fill_gradient(
    low = '#ededed',
    high = 'tomato3',
    breaks = seq(0, 20, by = 4)
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = 'Variance'
  ) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    legend.title = element_text(size = 18, face = 'italic')
    
  )
venom_mRNA_variance_heatmap

# Create venom expression heat map
all_venom_genes_mRNA_heatmap <- ggplot(all_venom_genes_mRNA_df, aes(y = sample.id, x = genes, fill = mRNA.ntd)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(
    y = 'Individuals', 
    x = 'Genes', 
    fill = 'mRNA expression' #, 
    # title = "mRNA expression Level for all Targeted Venom genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 15),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margins
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 18, face = 'italic'),
    legend.title = element_text(size = 18, face = 'italic')
  ) 
all_venom_genes_mRNA_heatmap
ggsave("Figures/Heat_Maps/mRNA/Venom_mRNA_heatmap_2025.01.23.pdf", plot = all_venom_genes_mRNA_heatmap, width = 14, height = 4, dpi = 900, create.dir = T)


# Create a combined heatmap of venom genes and their variance
venom_genes_with_variance_mRNA_heatmap <- plot_grid(venom_mRNA_variance_heatmap, all_venom_genes_mRNA_heatmap, ncol = 1, align ='v', rel_heights = c(0.1, 1))
venom_genes_with_variance_mRNA_heatmap
# Create title
venom_genes_with_variance_mRNA_heatmap <- ggdraw() +
  draw_label("a. Venom gene mRNA sample-wide expression heatmap with variances", fontface = 'bold', size = 14, x = 0.5, y = 0.99, hjust = 0.5) +
  draw_plot(venom_genes_with_variance_mRNA_heatmap, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
venom_genes_with_variance_mRNA_heatmap
ggsave("Figures/Heat_Maps/mRNA/Venom_mRNA_heatmap_with_variance_2025.01.23.pdf", plot = venom_genes_with_variance_mRNA_heatmap, width = 14, height = 8, dpi = 900, create.dir = T)


# Protein Heat Maps ----

## All protein expression ----

# Prepare data for heat-mapification
protein_all_df <- mi_df %>% 
  filter(in.library == 'Yes') %>% # Remove any proteins not in the library sent to Anthony
  dplyr::select(
    sample.id, genes, venom.family, intensity
  ) %>% 
  distinct()

# Calculate variance in protein expression
protein_variance_df <- protein_all_df %>%
  group_by(genes) %>%
  summarize(variance = var(intensity)) %>%
  ungroup()


# Create a heatmap of variances per gene
protein_variance_heatmap <- ggplot(protein_variance_df, aes(x = genes, y = 1, fill = log10(variance + 1))) +
  geom_tile() +
  # scale_fill_viridis_c(option = 'viridis') +
  scale_fill_gradient(
    low = '#ededed',
    high = 'tomato3',
    breaks = seq(0, 30, by = 4)
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = 'Variance (Log)'
  ) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank()   # Remove y-axis ticks
  )
protein_variance_heatmap


# Create heat map
protein_all_heatmap <- ggplot(protein_all_df, aes(y = sample.id, x = genes, fill = log10(intensity + 1))) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(
    y = 'Individuals', 
    x = 'Proteins', 
    fill = 'Protein expression' #, 
    # title = "Protein expression Level for all Targeted genes (log scaled)"
  ) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)  # Remove margins
  )
protein_all_heatmap
ggsave("Figures/Heat_Maps/Protein/Protein_all_heatmap_2025.01.23.pdf", plot = protein_all_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)

# Combined the main and variance heatmap
protein_combined_heatmap <- plot_grid(protein_variance_heatmap, protein_all_heatmap, ncol = 1, align = 'v', rel_heights = c(0.1, 1))
protein_combined_heatmap
ggsave("Figures/Heat_Maps/Protein/Protein_all_heatmap_with_variance_2025.01.23.pdf", plot = protein_combined_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)


## Venom protein expression ----

# Create an order for the genes
venom_gene_order <- c(
  'BPP',
  'myotoxin',
  'ohanin',
  'PLA2B1', 'PLA2K', 'PLA2C1', 'PLA2A1',
  'SVSP1', 'SVSP2', 'SVSP3', 'SVSP4', 'SVSP10', 'SVSP5', 'SVSP6', 'SVSP11', 'SVSP7', 'SVSP8', 'SVSP9',
  'SVMP1', 'SVMP2', 'SVMP3', 'SVMP4', 'SVMP5', 'SVMP6', 'SVMP7', 'SVMP8', 'SVMP9', 'SVMP10', 'SVMP11', 'SVMP12',
  'CRISP1', 'CRISP2', 'CRISP3', 'CRISP4',
  'CTL1', 'CTL2', 'CTL3', 'CTL4', 'CTL5', 'CTL6',
  'EXO1', 'EXO2', 'EXO3',
  'LAAO1', 'LAAO2', 'LAAO3',
  'VEGF1', 'VEGF2',
  'vQC1', 'vQC2'
)

# Create an order for the samples
sample_order <- c(
  'CV0987_lutosus', 'CV0985_concolor', 'CV1086_viridis',
  'CV1081_viridis', 'CV0857_viridis', 'CV1087_viridis'
)

# Now for only venom
# Prepare data
venom_protein_df <- protein_all_df %>% 
  filter(
    str_starts(genes, "Venom"),
    !str_detect(genes, 'ADAM')
  ) %>%
  mutate(genes = str_remove(genes, '^Venom_')) %>% 
  distinct()

# Set the order of venom genes
venom_protein_df$genes <- factor(venom_protein_df$genes, levels = venom_gene_order)

# Set the order of the sample genes
venom_protein_df$sample.id <- factor(venom_protein_df$sample.id, levels = sample_order)

# Get the variance for all of the venom genes
venom_protein_variance_df <- venom_protein_df %>% 
  group_by(genes) %>% 
  summarize(variance = var(intensity)) %>% 
  ungroup()

# Create venom variance heatmap
venom_protein_variance_heatmap <- ggplot(venom_protein_variance_df, aes(x = genes, y = 1, fill = log(variance + 1))) +
  geom_tile() +
  # scale_fill_viridis_c(option = 'viridis') +
  scale_fill_gradient(
    low = '#ededed',
    high = 'tomato3',
    breaks = seq(0, 55, by = 10)
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = 'Variance (Log)'
  ) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    legend.title = element_text(size = 18, face = 'italic')
  )
venom_protein_variance_heatmap

# Create heat map
venom_protein_heatmap <- ggplot(venom_protein_df, aes(y = sample.id, x = genes, fill = log(intensity + 1))) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(
    y = 'Individuals',
    x = 'Proteins', 
    fill = 'Protein expression' #, 
    # title = "Protein expression Level for all Targeted Venom genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 15),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margins
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 18, face = 'italic'),
    legend.title = element_text(size = 18, face = 'italic')
  )
venom_protein_heatmap
ggsave("Figures/Heat_Maps/Protein/Venom_protein_heatmap_2025.01.23.pdf", plot = venom_protein_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)

# Create a combined heatmap of venom genes and their variance
venom_protein_with_variance_heatmap <- plot_grid(venom_protein_variance_heatmap, venom_protein_heatmap, ncol = 1, align ='v', rel_heights = c(0.1, 1))
venom_protein_with_variance_heatmap
# Create title
venom_protein_with_variance_heatmap <- ggdraw() +
  draw_label("b. Venom gene protein sample-wide expression heatmap with variances", fontface = 'bold', size = 14, x = 0.5, y = 0.99, hjust = 0.5) +
  draw_plot(venom_protein_with_variance_heatmap, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
venom_protein_with_variance_heatmap
ggsave("Figures/Heat_Maps/Protein/Venom_protein_heatmap_with_variance_2025.01.23.pdf", plot = venom_protein_with_variance_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)



# miRNA Cluster Heat Maps ----

## All miRNAs ----

# Prepare data for heat-mapification
all_miRNA_df <- mi_df %>% 
  filter(
    in.library == 'Yes',
    !is.na(miRNA.cluster)
  ) %>% # Remove any proteins not in the library sent to Anthony
  dplyr::select(
    sample.id, genes, venom.family, miRNA.cluster, miRNA.ntd
  ) %>% 
  distinct() %>% 
  mutate(miRNA.cluster = str_remove(miRNA.cluster, '^cvi-')) # Remove the 'cvi-' in front of each miRNA name

# Use pheatmap to sort the figure
# Prepare the data in a format pheatmap can use
miRNA_wide_df <- all_miRNA_df %>% 
  dplyr::select(
    -genes, -venom.family
  ) %>% 
  distinct() %>% 
  pivot_wider(
    names_from = sample.id,
    values_from = miRNA.ntd
  )

# Convert miRNA expression to a matrix
miRNA_matrix <- as.matrix(
  miRNA_wide_df %>% 
    column_to_rownames(var = 'miRNA.cluster') %>% 
    distinct()
)

# Save as a high resolution PNG
pdf('Figures/Heat_Maps/miRNA/Pheatmap_miRNA_heatmap_2025.01.23.pdf', width = 10, height = 10)
# Create heatmap using pheatmap
pheatmap(
  miRNA_matrix,
  scale = 'none',
  cluster_cols = T,
  cluster_rows = T,
  show_rownames = T,
  show_colnames = T,
  fontsize_row = 3,
  fontsize_col = 10
)
dev.off()

# Save the pheatmap in memory to an object and extract row order
pheatmap_miRNA_res <- pheatmap(
  miRNA_matrix,
  scale = 'none',
  cluster_rows = T,
  cluster_cols = T
)

# Extract the row order (miRNA cluster order)
miRNA_cluster_order <- rownames(miRNA_matrix)[pheatmap_miRNA_res$tree_row$order]

# Reorder miRNA cluster in the df
all_miRNA_df$miRNA.cluster <- factor(all_miRNA_df$miRNA.cluster, levels = miRNA_cluster_order)

# Set the order of the sample genes
all_miRNA_df$sample.id <- factor(all_miRNA_df$sample.id, levels = sample_order)

# Calculate variance in miRNA
miRNA_variance_df <- all_miRNA_df %>%
  dplyr::select(-genes, -venom.family) %>% 
  distinct() %>% 
  group_by(miRNA.cluster) %>%
  summarize(variance = var(miRNA.ntd)) %>% 
  ungroup()

# Drop genes column so that it doesn't overwrite the clusters multiple times
# Apply this to the heat maps
all_miRNA_df2 <- all_miRNA_df %>% dplyr::select(-genes, -venom.family) %>% distinct()

# Create a heatmap of variances per miRNA locus
miRNA_variance_heatmap <- ggplot(miRNA_variance_df, aes(x = miRNA.cluster, y = 1, fill = variance)) +
  geom_tile() +
  # scale_fill_viridis_c(option = 'viridis') +
  scale_fill_gradient(
    low = '#ededed',
    high = 'tomato3',
    breaks = seq(0, 22, by = 5)
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = 'Variance (Log)'
  ) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    legend.title = element_text(size = 16, face = 'italic')
  )
miRNA_variance_heatmap

# Create heat map
miRNA_all_heatmap <- ggplot(all_miRNA_df2, aes(y = sample.id, x = miRNA.cluster, fill = miRNA.ntd)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(
    y = 'Individuals', 
    x = 'miRNA Cluster', 
    fill = 'miRNA expression' #, 
    # title = "miRNA expression Level for all miRNA Clusters (Log Scaled)"
  ) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margins
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 18, face = 'italic'),
    legend.title = element_text(size = 18, face = 'italic')
  )
miRNA_all_heatmap
ggsave("Figures/Heat_Maps/miRNA/miRNA_all_heatmap_2025.01.23.pdf", plot = miRNA_all_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)

# Create a combined heatmap of miRNAs expression and their variance
miRNA_combined_heatmap <- plot_grid(miRNA_variance_heatmap, miRNA_all_heatmap, ncol = 1, align = 'v', rel_heights = c(0.1, 1))
miRNA_combined_heatmap
ggsave("Figures/Heat_Maps/miRNA/miRNA_all_with_variance_heatmap_2025.01.23.pdf", plot = miRNA_combined_heatmap, width = 18, height = 6, dpi = 900, create.dir = T)



## Venom miRNAs ----

# Create a data frame with only venom targeting miRNAs
venom_miRNA_df <- all_miRNA_df %>% 
  dplyr::filter(
    str_detect(genes, 'Venom_')
  )

# Prepare the data in a format pheatmap can use
venom_miRNA_wide_df <- venom_miRNA_df %>% 
  dplyr::select(
    -genes, -venom.family
  ) %>% 
  distinct() %>% 
  pivot_wider(
    names_from = sample.id,
    values_from = miRNA.ntd
  )

# Convert miRNA expression to a matrix
venom_miRNA_matrix <- as.matrix(
  venom_miRNA_wide_df %>% 
    column_to_rownames(var = 'miRNA.cluster') %>% 
    distinct()
)

# Save as a high resolution PNG
pdf('Figures/Heat_Maps/miRNA/Venom_Pheatmap_miRNA_heatmap_2025.01.23.pdf', width = 10, height = 10)
# Create heatmap using pheatmap
pheatmap(
  venom_miRNA_matrix,
  scale = 'none',
  cluster_cols = T,
  cluster_rows = T,
  show_rownames = T,
  show_colnames = T,
  fontsize_row = 3,
  fontsize_col = 10
)
dev.off()

# Save the pheatmap in memory to an object and extract row order
venom_pheatmap_miRNA_res <- pheatmap(
  venom_miRNA_matrix,
  scale = 'none',
  cluster_rows = T,
  cluster_cols = T
)

# Extract the row order (miRNA cluster order)
venom_miRNA_cluster_order <- rownames(venom_miRNA_matrix)[venom_pheatmap_miRNA_res$tree_row$order]

# Reorder miRNA cluster in the df
venom_miRNA_df$miRNA.cluster <- factor(venom_miRNA_df$miRNA.cluster, levels = venom_miRNA_cluster_order)

# Set the order of the sample genes
venom_miRNA_df$sample.id <- factor(venom_miRNA_df$sample.id, levels = sample_order)

# Calculate variance in miRNA
venom_miRNA_variance_df <- venom_miRNA_df %>%
  dplyr::select(-genes, -venom.family) %>% 
  distinct() %>% 
  group_by(miRNA.cluster) %>%
  summarize(variance = var(miRNA.ntd)) %>% 
  ungroup()

# Drop genes column so that it doesn't overwrite the clusters multiple times
venom_miRNA_df2 <- venom_miRNA_df %>% dplyr::select(-genes, -venom.family) %>% distinct()


# Create venom variance heatmap
venom_miRNA_variance_heatmap <- ggplot(venom_miRNA_variance_df, aes(x = miRNA.cluster, y = 1, fill = variance)) +
  geom_tile() +
  # scale_fill_viridis_c(option = 'viridis') +
  scale_fill_gradient(
    low = '#ededed',
    high = 'tomato3',
    breaks = seq(0, 22, by = 5)
  ) +
  labs(
    y = NULL,
    x = NULL,
    fill = 'Variance'
  ) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),   # Remove y-axis ticks
    legend.title = element_text(size = 18, face = 'italic')
  )
venom_miRNA_variance_heatmap


# Create heat map
venom_miRNA_heatmap <- ggplot(venom_miRNA_df2, aes(y = sample.id, x = miRNA.cluster, fill = miRNA.ntd)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(
    y = 'Individuals', 
    x = 'miRNA Cluster', 
    fill = 'miRNA expression' #, 
    # title = "miRNA expression Level for all Venom targeting miRNA Clusters (Log Scale)"
  ) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Remove margins
    axis.text.y = element_text(size = 15),
    axis.title = element_text(size = 18, face = 'italic'),
    legend.title = element_text(size = 18, face = 'italic')
  )
venom_miRNA_heatmap
ggsave("Figures/Heat_Maps/miRNA/Venom_miRNA_heatmap_2025.01.23.pdf", plot = venom_miRNA_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)

# Create a combined heatmap of venom genes and their variance
venom_miRNAs_with_variance_heatmap <- plot_grid(venom_miRNA_variance_heatmap, venom_miRNA_heatmap, ncol = 1, align ='v', rel_heights = c(0.1, 1))
venom_miRNAs_with_variance_heatmap
# Create title
venom_miRNAs_with_variance_heatmap <- ggdraw() +
  draw_label("c. Venom gene targeting miRNAs sample-wide expression heatmap with variances", fontface = 'bold', size = 18, x = 0.5, y = 0.99, hjust = 0.5) +
  draw_plot(venom_miRNAs_with_variance_heatmap, y = 0, height = 0.9)  # Adjust height and y to fit the title and plots
venom_miRNAs_with_variance_heatmap
ggsave("Figures/Heat_Maps/miRNA/Venom_miRNA_heatmap_with_variance_2025.01.23.pdf", plot = venom_miRNAs_with_variance_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)



# Figure b ----

# Create figure b
figure_b <- plot_grid(
  venom_genes_with_variance_mRNA_heatmap,
  venom_protein_with_variance_heatmap,
  venom_miRNAs_with_variance_heatmap,
  ncol = 1,
  align = 'v'
)
figure_b
ggsave("Figures/1_Main_Figures/Figure_b/Figure_b_2025.03.04.pdf", plot = figure_b, width = 18, height = 16, dpi = 2400, create.dir = T)
