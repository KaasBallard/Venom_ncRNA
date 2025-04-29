# Last Edited: 2025/02/03

# Set up and Read Data ----

## Load in packages ----
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
library(arrow)

## Read data ----

# Create variable for the fused data set.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read both in as data frames
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)
glimpse(miRNA_mRNA_protein_df)

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    !str_detect(genes, 'maker-scaffold|augustus|XP_'),
    !is.na(feature.type)
  ) %>% 
  distinct()
rm(miRNA_mRNA_protein_df)


# Create color scheme for the venom genes
SVMP_color <- '#4A70B5'
ADAM_color <- '#9A70B5'
SVSP_color <- '#F0B830' 
PLA2_color <- '#7570B3'
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
miRNA_color <- 'grey'
three_prime_color <- 'black'
# five_prime_color <- '#0072b2'
five_prime_color <- '#1B9E77'
# cds_color <- '#d55e00'
cds_color <- '#4A70B5'



# Print Number of miRNA Relationships ----

# Create a data frame with only miRNA relationships
relationshps_df <- mi_df %>% 
  dplyr::select(
    miRNA.cluster, miRNA.start, miRNA.end, genes, feature.type, miRNA.target.start, miRNA.target.end, positions, total.energy, total.score
  ) %>% 
  distinct()

# Print number of putative miRNA relationships
print(nrow(relationshps_df))

# Create a data frame with only miRNA-venom relationships
venom_relationshps_df <- mi_df %>% 
  filter(str_detect(genes, 'Venom_')) %>%
  dplyr::select(
    miRNA.cluster, miRNA.start, miRNA.end, genes, feature.type, miRNA.target.start, miRNA.target.end, positions, total.energy, total.score
  ) %>% 
  distinct()

# Create a data frame with only miRNA-non-venom relationships
non_venom_relationshps_df <- mi_df %>% 
  filter(!str_detect(genes, 'Venom_')) %>%
  dplyr::select(
    miRNA.cluster, miRNA.start, miRNA.end, genes, feature.type, miRNA.target.start, miRNA.target.end, positions, total.energy, total.score
  ) %>% 
  distinct()



# Print number of putative miRNA relationships
print(nrow(venom_relationshps_df))

# The below is to find out why adding more than just the first three columns (miRNA.cluster, genes, and feature.type) increases row numbers
# Data frame with total.energy and total.score
with_energy_score_df <- mi_df %>%
  # filter(str_detect(genes, 'Venom_')) %>%
  dplyr::select(
    miRNA.cluster, genes, feature.type, positions, total.energy, total.score
  ) %>%
  distinct()

# Data frame without total.energy and total.score
without_energy_score_df <- mi_df %>%
  # filter(str_detect(genes, 'Venom_')) %>%
  dplyr::select(miRNA.cluster, genes, feature.type, positions) %>%
  distinct()

# Rows in 'with_energy_score_df' that are not in 'without_energy_score_df'
diff_in_with_energy_score <- anti_join(with_energy_score_df, without_energy_score_df, by = c("miRNA.cluster", "genes", "feature.type"))

# Rows in 'without_energy_score_df' that are not in 'with_energy_score_df'
diff_in_without_energy_score <- anti_join(without_energy_score_df, with_energy_score_df, by = c("miRNA.cluster", "genes", "feature.type"))

# Count the number of differences
nrow(diff_in_with_energy_score)
nrow(diff_in_without_energy_score)

# View the differences
print(diff_in_with_energy_score)
print(diff_in_without_energy_score)



# Count distinct rows with total.energy and total.score
count_with_energy_score <- mi_df %>%
  distinct(miRNA.cluster, genes, feature.type, total.energy, total.score) %>%
  nrow()

# Count distinct rows without total.energy and total.score
count_without_energy_score <- mi_df %>%
  distinct(miRNA.cluster, genes, feature.type) %>%
  nrow()

# Print the counts
count_with_energy_score
count_without_energy_score

# Find rows where the same miRNA.cluster, genes, and feature.type have multiple total.energy and total.score
multiple_values <- mi_df %>%
  group_by(miRNA.cluster, genes, feature.type) %>%
  summarise(
    count = n_distinct(total.energy, total.score)
  ) %>%
  filter(count > 1)

# View rows with multiple total.energy and total.score values
print(multiple_values)

# Inspect the rows where multiple total.energy and total.score exist
problematic_rows <- mi_df %>%
  filter(str_detect(genes, 'Venom_')) %>%
  semi_join(multiple_values, by = c("miRNA.cluster", "genes", "feature.type"))

# View the problematic rows
print(problematic_rows)


# Target Types Graph ----

## All targets ----
# Create a data frame with target type information by grouping by origin column
target_type_df <- relationshps_df %>% 
  group_by(feature.type) %>% 
  summarise(Counts = n()) %>%
  ungroup() %>% 
  mutate(
    color = case_when(
      grepl('three_prime_utr', feature.type) ~ three_prime_color,
      grepl('five_prime_utr', feature.type) ~ five_prime_color,
      grepl('CDS', feature.type) ~ cds_color
    )
  ) %>% 
  mutate(
    feature.type = case_when(
      grepl('three_prime_utr', feature.type) ~ "3' UTR",
      grepl('five_prime_utr', feature.type) ~ "5' UTR",
      grepl('CDS', feature.type) ~ "Coding Sequences",
      grepl('Total', feature.type) ~ "Total"
    )
  )

# Add a Total row using add_row, ensuring sum is calculated correctly
target_type_df <- target_type_df %>%
  add_row(feature.type = "Total targets", Counts = sum(target_type_df$Counts), color = 'grey')


# Create a bar plot that quantifies the relationships of each of the target types
target_types_bar_plot <- ggplot(target_type_df, aes(x = feature.type, y = Counts)) +
  geom_bar(stat = 'identity', aes(fill = color), width = 0.5) +
  geom_text(aes(label = Counts), vjust = -0.5, color = "black", size = 2.5) +  # Add count labels above bars
  labs(
    title = "d. miRNA - gene relationships", 
    x = "miRNA target type", 
    y = "Number of relationships"
  ) +
  theme_classic2() +
  scale_fill_identity() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.text.x = element_text(size = 7),  # Adjust x-axis text size
    axis.text.y = element_text(size = 7),   # Adjust y-axis text size
    axis.title = element_text(size = 8)
  )
target_types_bar_plot


## Venom targets ----
# Create a data frame with target type information by grouping by origin column
venom_target_type_df <- venom_relationshps_df %>% 
  group_by(feature.type) %>% 
  summarise(Counts = n()) %>% 
  ungroup() %>% 
  mutate(
    color = case_when(
      grepl('three_prime_utr', feature.type) ~ three_prime_color,
      grepl('five_prime_utr', feature.type) ~ five_prime_color,
      grepl('CDS', feature.type) ~ cds_color,
    )
  ) %>% 
  mutate(
    feature.type = case_when(
      grepl('three_prime_utr', feature.type) ~ "3' UTR",
      grepl('five_prime_utr', feature.type) ~ "5' UTR",
      grepl('CDS', feature.type) ~ "Coding Sequences",
    )
  )
# Add a Total row using add_row, ensuring sum is calculated correctly
venom_target_type_df <- venom_target_type_df %>%
  add_row(feature.type = "Total targets", Counts = sum(venom_target_type_df$Counts), color = 'orange')



# Create a bar plot that quantifies the relationships of each of the target types
venom_target_types_bar_plot <- ggplot(venom_target_type_df, aes(x = feature.type, y = Counts)) +
  geom_bar(stat = 'identity', aes(fill = color), width = 0.5) +
  geom_text(aes(label = Counts), vjust = -0.5, color = "black", size = 3) +  # Add count labels above bars
  labs(
    title = "b. miRNA - venom gene relationships", 
    x = "miRNA target type", 
    y = "Number of relationships"
  ) +
  theme_classic2() +
  scale_fill_identity() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.text.x = element_text(size = 7),  # Adjust x-axis text size
    axis.text.y = element_text(size = 7),   # Adjust y-axis text size
    axis.title = element_text(size = 8)
  )
venom_target_types_bar_plot


## Non-venom targets ----
# Create a data frame with target type information by grouping by origin column
non_venom_target_type_df <- non_venom_relationshps_df %>% 
  group_by(feature.type) %>% 
  summarise(Counts = n()) %>%
  ungroup() %>% 
  mutate(
    color = case_when(
      grepl('three_prime_utr', feature.type) ~ three_prime_color,
      grepl('five_prime_utr', feature.type) ~ five_prime_color,
      grepl('CDS', feature.type) ~ cds_color
    )
  ) %>% 
  mutate(
    feature.type = case_when(
      grepl('three_prime_utr', feature.type) ~ "3' UTR",
      grepl('five_prime_utr', feature.type) ~ "5' UTR",
      grepl('CDS', feature.type) ~ "Coding Sequences",
      grepl('Total', feature.type) ~ "Total"
    )
  )

# Add a Total row using add_row, ensuring sum is calculated correctly
non_venom_target_type_df <- non_venom_target_type_df %>%
  add_row(feature.type = "Total targets", Counts = sum(non_venom_target_type_df$Counts), color = 'yellow')


# Create a bar plot that quantifies the relationships of each of the target types
non_venom_target_types_bar_plot <- ggplot(non_venom_target_type_df, aes(x = feature.type, y = Counts)) +
  geom_bar(stat = 'identity', aes(fill = color), width = 0.5) +
  geom_text(aes(label = Counts), vjust = -0.5, color = "black", size = 2.5) +  # Add count labels above bars
  labs(
    title = "c. miRNA - non-venom gene relationships", 
    x = "miRNA target type", 
    y = "Number of relationships"
  ) +
  theme_classic2() +
  scale_fill_identity() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.text.x = element_text(size = 7),  # Adjust x-axis text size
    axis.text.y = element_text(size = 7),   # Adjust y-axis text size
    axis.title = element_text(size = 8)
  )
non_venom_target_types_bar_plot


# Number of miRNAs plot ----

# Create a data frame that only contains venom targeting miRNAs
venom_mirnas_df <- mi_df %>% 
  filter(
    str_detect(genes, 'Venom_')
  ) %>% 
  select(miRNA.cluster) %>% 
  distinct() %>% 
  summarise(
    Counts = n()
  ) %>% 
  mutate(
    Target.Gene = 'Venom'
  ) %>% 
  mutate(color = 'orange')

# # Reset index
# row.names(venom_mirnas_df) <- NULL # Don't actually need this

# Create a data frame that only contains non-venom targeting miRNAs
non_venom_mirnas_df <- mi_df %>% 
  filter(
    !str_detect(genes, 'Venom_')
  ) %>% 
  select(miRNA.cluster) %>% 
  distinct() %>% 
  summarise(
    Counts = n()
  ) %>% 
  mutate(
    Target.Gene = 'Non-Venom'
  ) %>% 
  mutate(color = 'yellow')

# Create a data frame that contains all miRNAs
all_mirnas_df <- mi_df %>% 
  select(miRNA.cluster) %>% 
  distinct() %>% 
  summarise(
    Counts = n()
  ) %>% 
  mutate(
    Target.Gene = 'Total putative miRNAs'
  ) %>% 
  mutate(color = 'grey')

# Create an order for the target types
order <- c('Venom', 'Non-Venom', 'Total putative miRNAs')

# Bind the rows
mirnas_df <- bind_rows(all_mirnas_df, venom_mirnas_df, non_venom_mirnas_df) %>% 
  mutate(Target.Gene = factor(Target.Gene, levels = order)) %>% 
  arrange(Target.Gene)

# Plot
mirnas_plot <- ggplot(mirnas_df, aes(x = Target.Gene, y = Counts)) +
  geom_bar(stat = 'identity', aes(fill = color), width = 0.5) +
  geom_text(aes(label = Counts), vjust = -0.5, color = "black", size = 3) +  # Add count labels above bars
  labs(
    title = "a. miRNAs by target gene type", 
    x = "miRNA targets", 
    y = "miRNAs"
  ) +
  theme_classic2() +
  scale_fill_identity() +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )
mirnas_plot

# Create a combined plot
three_graphs <- plot_grid(
  venom_target_types_bar_plot,
  non_venom_target_types_bar_plot,
  target_types_bar_plot,
  ncol = 1,
  nrow = 3,
  align = 'v',
  rel_widths = c(1, 1, 1),
  rel_heights = c(1, 1, 1)
)
three_graphs

# Create another combined plot
four_graphs <- plot_grid(
  mirnas_plot,
  three_graphs,
  align = 'v'
)
four_graphs
ggsave("Figures/miRNA_Numbers_and_Relationships/Unfiltered/miRNA_Counts_and_relationships_2025.02.03.pdf", plot = four_graphs, width = 9, height = 10, dpi = 900, create.dir = T)





# Find miRNA that doesn't target a venom gene ----

# Create a data frame for only miRNAs that target venom genes
venom_only_df <- mi_df %>% 
  filter(str_detect(genes, 'Venom_')) %>% 
  select(
    miRNA.cluster, miRNA.cluster.original
  ) %>% 
  distinct()
# Reset index
rownames(venom_only_df) <- NULL

# Create a data frame for only miRNAs that target non-venom genes
non_venom_only_df <- mi_df %>% 
  filter(!str_detect(genes, 'Venom_')) %>% 
  select(
    miRNA.cluster, miRNA.cluster.original
  ) %>% 
  distinct()
# Reset index
rownames(non_venom_only_df) <- NULL

# Find miRNAs that target only non-venom genes
miRNAs_non_venom_only <- anti_join(non_venom_only_df, venom_only_df, by = c("miRNA.cluster", "miRNA.cluster.original"))

print(miRNAs_non_venom_only)
# The nonvenom miRNA should be Cluster_1888



# miRNAs Chromosome Locations ----

# Create a data frame for miRNAS and their chromosomes
chromosome_df <- mi_df %>%
  # filter(Protein.Observed == 'Yes') %>% 
  # filter(feature.type == 'three_prime_utr') %>% 
  select(
    miRNA.cluster, miRNA.cluster.original, feature.type, miRNA.start, miRNA.end, miRNA.sequence.chrom, genes, miRNA.target.chrom, miRNA.target.start, miRNA.target.end, total.energy, total.score
  ) %>% 
  distinct() %>% 
  mutate(
    Chrom.Match = if_else(miRNA.sequence.chrom == miRNA.target.chrom, 'match', 'miss-match')
  ) %>% 
  group_by(Chrom.Match) %>% 
  summarise(Counts = n()) %>% 
  mutate(Percentage = round(Counts / sum(Counts) * 100, 1))
  
# Create plot based on this data
matching_chromosomes_pie_chart <- ggplot(chromosome_df, aes(x = "", y = Counts, fill = Chrom.Match)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  geom_text(
    aes(label = paste0(Percentage, "% (", Counts, ")")),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = 'black'
  ) +
  coord_polar(theta = 'y', start = 0) +
  # scale_fill_viridis(discrete = T, option = 'viridis') +
  theme_void() +
  # scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 12),
    # legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Match vs miss-match', 
    title = 'a. Chromosomal miRNA - gene matches vs miss-matches'
  )
matching_chromosomes_pie_chart


# Create a data frame for miRNAS and their chromosomes for non_venom genes only
non_venom_chrom_df <- mi_df %>%
  filter(!str_detect(genes, 'Venom_')) %>% 
  select(
    miRNA.cluster, miRNA.cluster.original, feature.type, miRNA.start, miRNA.end, miRNA.sequence.chrom, genes, miRNA.target.chrom, miRNA.target.start, miRNA.target.end, total.energy, total.score
  ) %>% 
  distinct() %>% 
  mutate(
    Chrom.Match = if_else(miRNA.sequence.chrom == miRNA.target.chrom, 'match', 'miss-match')
  ) %>% 
  group_by(Chrom.Match) %>% 
  summarise(Counts = n()) %>% 
  mutate(Percentage = round(Counts / sum(Counts) * 100, 1))

# Create plot based on this data
non_venom_matching_chrom_pie_chart <- ggplot(non_venom_chrom_df, aes(x = "", y = Counts, fill = Chrom.Match)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  geom_text(
    aes(label = paste0(Percentage, "% (", Counts, ")")),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = 'black'
  ) +
  coord_polar(theta = 'y', start = 0) +
  # scale_fill_viridis(discrete = T, option = 'viridis') +
  theme_void() +
  # scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 12),
    # legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Match vs miss-match', 
    title = 'b. Chromosomal miRNA - non-venom gene matches vs miss-matches'
  )
non_venom_matching_chrom_pie_chart


# Create a data frame for miRNAS and their chromosomes for non_venom genes only
venom_chrom_df <- mi_df %>%
  filter(str_detect(genes, 'Venom_')) %>% 
  # filter(
  #   total.score > 150,
  #   total.energy < -19
  # ) %>% # Use if I want to see if filtering by score changes anything
  select(
    miRNA.cluster, miRNA.cluster.original, feature.type, miRNA.start, miRNA.end, miRNA.sequence.chrom, genes, miRNA.target.chrom, miRNA.target.start, miRNA.target.end, total.energy, total.score
  ) %>% 
  distinct() %>% 
  mutate(
    Chrom.Match = if_else(miRNA.sequence.chrom == miRNA.target.chrom, 'match', 'miss-match')
  ) %>% 
  group_by(Chrom.Match) %>% 
  summarise(Counts = n()) %>% 
  mutate(Percentage = round(Counts / sum(Counts) * 100, 1))

# Create plot based on this data
venom_matching_chrom_pie_chart <- ggplot(venom_chrom_df, aes(x = "", y = Counts, fill = Chrom.Match)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  geom_text(
    aes(label = paste0(Percentage, "% (", Counts, ")")),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = 'black'
  ) +
  coord_polar(theta = 'y', start = 0) +
  # scale_fill_viridis(discrete = T, option = 'viridis') +
  theme_void() +
  # scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 14),
    # legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) +  # Adjust the margin and other title information
  labs(
    fill = 'Match vs miss-match', 
    title = 'c. Chromosomal miRNA - venom gene matches vs miss-matches'
  )
venom_matching_chrom_pie_chart

# Create a combined figure for this
chrom_pie_charts <- plot_grid(
  matching_chromosomes_pie_chart,
  non_venom_matching_chrom_pie_chart,
  venom_matching_chrom_pie_chart,
  ncol = 3,
  nrow = 1,
  align = 'v',
  rel_widths = c(1, 1, 1),
  rel_heights = c(1, 1, 1)
)
chrom_pie_charts
ggsave('Figures/Pie_Charts/Unfiltered/miRNA_Target_Chromosomal_Matches/miRNA-Gene_MatchesVsMissmatches_2025.02.03.pdf', plot = chrom_pie_charts, width = 20, height = 6, dpi = 900, create.dir = T)



# Create a combined figure of the venom relationships bar graph and the piechart
# First, I need to change the bar graph title
venom_target_types_bar_plot2 <- venom_target_types_bar_plot + 
  labs(title = 'a. miRNA - venom gene relationships') + 
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 14),
    axis.text.x = element_text(size = 10),  # Adjust x-axis text size
    axis.text.y = element_text(size = 10),   # Adjust y-axis text size
    axis.title = element_text(size = 11)
  )
venom_target_types_bar_plot2

# Second, change the pie chart title
venom_matching_chrom_pie_chart2 <- venom_matching_chrom_pie_chart + 
  labs(title = 'b. Chromosomal miRNA - venom gene matches vs miss-matches') +
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5, margin = margin(b = 5, t = 5), size = 14),
    legend.title = element_text(face = 'bold', size = 11),
    legend.text = element_text(size = 11)
  )
venom_matching_chrom_pie_chart2

# Combine
fig_s <- plot_grid(
  venom_target_types_bar_plot2,
  venom_matching_chrom_pie_chart2,
  align = 'v'
)
fig_s
ggsave('Figures/1_Main_Figures/Figure_s/Figure_s_2025.02.03.pdf', plot = fig_s, width = 14, height = 7, dpi = 1800, create.dir = T)
