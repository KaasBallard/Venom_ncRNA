#### Set up and Read Data ####

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

# Set working directory
# setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/')
setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")

# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/miRNA_mRNA_Protein_Combined_Data_IMPORTANT.2024.06.14.tsv'

# Read both in as data frames
miRNA_mRNA_protein_df <- read.table(file = miRNA_mRNA_protein_data, header = T)
# I am going to exclude the following genes because they are not well annotated
excluded_genes = c(
  'maker-scaffold-mi1-augustus-gene-59.13_crovir-transcript-12940',
  'maker-scaffold-mi1-augustus-gene-59.20_crovir-transcript-12947',
  'maker-scaffold-mi2-augustus-gene-22.17_crovir-transcript-739',
  'maker-scaffold-un11-augustus-gene-5.19',
  'XP_016876419',
  'XP_011528471'
)


# Let's keep the analysis limited to the 3UTR
# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  dplyr::rename(
    'miRNA.Cluster.Original' = 'miRNA.Cluster',
    'Genes' = 'Converted.Gene.IDs',
    'miRNA.Cluster' = 'Putative.miRNA.Name'
  ) %>%
  dplyr::select(miRNA.Cluster, everything()) %>% # Move the new miRNA.Clusters to the front
  filter(!(Genes %in% excluded_genes)) %>% 
  filter(Origin == 'three_prime_utr') %>% # Filter out everything not targeting the 3' UTR
  dplyr::select(-Origin)  # Remove to save memory
rm(miRNA_mRNA_protein_df)

# Create color scheme for the venom genes
SVMP_color <- '#4A70B5'
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

  
#### Check # of miRNAs ####

# Check all miRNAs 
miRNA_only_df <- mi_df %>% 
  dplyr::select(miRNA.Cluster) %>% distinct() %>% rowwise() %>% mutate(index = row_number()) %>% ungroup()
# Check venom targeting miRNAs
venom_miRNA_only_df <- mi_df %>% 
  filter(str_starts(Genes, 'Venom')) %>% 
  dplyr::select(miRNA.Cluster) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(index = row_number()) %>% 
  ungroup()


#### Gene Table ####

# Create table of genes regulated by miRNAs
gene_miRNA_table <- mi_df %>%
  distinct(Genes, miRNA.Cluster.Original, miRNA.Cluster) %>%
  group_by(Genes) %>%
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.Cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )

# Now save that table as a csv
write.csv(gene_miRNA_table, 'Tables/3UTR/Genes_regulated_by_miRNAs_in_3UTR_2024.06.23.csv', row.names = F)

# Now create the same thing but only for venom genes.
# Create table of genes regulated by miRNAs
venom_gene_miRNA_table <- mi_df %>%
  filter(str_starts(Genes, 'Venom_|PLA2G2E.1')) %>% 
  distinct(Genes, miRNA.Cluster.Original, miRNA.Cluster) %>%
  group_by(Genes) %>%
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.Cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
# Now save that table as a csv
write.csv(venom_gene_miRNA_table, 'Tables/3UTR/Venom_Genes_regulated_by_miRNAs_in_3UTR_2024.06.23.csv', row.names = F)

# Now create the same thing but only for non-venom genes.
# Create table of genes regulated by miRNAs
non_venom_gene_miRNA_table <- mi_df %>%
  filter(!str_starts(Genes, 'Venom')) %>% 
  distinct(Genes, miRNA.Cluster.Original, miRNA.Cluster) %>%
  group_by(Genes) %>%
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.Cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
# Now save that table as a csv
write.csv(non_venom_gene_miRNA_table, 'Tables/3UTR/Non-Venom_Genes_regulated_by_miRNAs_in_3UTR_2024.06.23.csv', row.names = F)

# Now I want to create the opposite thing, clusters with their genes
miRNA_gene_table <- mi_df %>% 
  distinct(Genes, miRNA.Cluster, miRNA.Cluster.Original) %>% 
  group_by(miRNA.Cluster) %>% 
  summarise(
    Genes = paste(Genes, collapse = ", "),
    Number.of.Genes.per.miRNA = n()
  ) %>% 
  dplyr::rename(
    'miRNA Locus' = 'miRNA.Cluster',
    'Number of Genes per miRNA' = 'Number.of.Genes.per.miRNA'
  )
# Now save that table as a csv
write.csv(miRNA_gene_table, 'Tables/3UTR/miRNAs_in_3UTR_and_regulated_genes_2024.06.23.csv', row.names = F)

# Create a table for the best blast hits for each miRNA
putative_miRNA_names <- mi_df %>% 
  distinct(miRNA.Cluster.Original, Best.miRNA.Blast.Hits, miRNA.Cluster) %>%
  group_by(miRNA.Cluster.Original) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.Cluster.Original',
    'Best Blast Hit' = 'Best.miRNA.Blast.Hits',
    'Putative miRNA Name' = 'miRNA.Cluster'
  )
write.csv(putative_miRNA_names, 'Tables/3UTR/Putative_miRNA_Names_2024.06.23.csv', row.names = F)

# Create a table for the best blast hits for each miRNA that targets a venom gene
venom_putative_miRNA_names <- mi_df %>% 
  filter(str_starts(Genes, 'Venom_|PLA2G2E.1')) %>% 
  distinct(miRNA.Cluster.Original, Best.miRNA.Blast.Hits, miRNA.Cluster) %>%
  group_by(miRNA.Cluster.Original) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.Cluster.Original',
    'Best Blast Hit' = 'Best.miRNA.Blast.Hits',
    'Putative miRNA Name' = 'miRNA.Cluster'
  )
write.csv(venom_putative_miRNA_names, 'Tables/3UTR/Venom_Putative_miRNA_Names_2024.06.23.csv', row.names = F)



#### Venom Transcriptome Profile Pie Charts ####

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


# Read mRNA data
mRNA_data <- 'Data/mRNA/RNAseq_VST_Formated_IMPORTANT_2024.06.14.tsv'
mRNA_df <- read.table(file = mRNA_data, header = T) %>% 
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove, collapse = "|"), '')) %>% 
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove2, collapse = "|"), '')) %>% 
  pivot_wider(
    names_from = Sample.ID,
    values_from = RNA.VST
  ) %>% 
  dplyr::rename('Genes' = 'Converted.Gene.IDs') %>% 
  filter(str_starts(Genes, 'Venom_')) %>% 
  filter(!str_starts(Genes, 'Venom_ADAM')) %>% 
  mutate(Venom.Family = case_when(grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
                                  # grepl('VEGF', Genes) ~ 'VEGF',
                                  # grepl('ohanin', Genes) ~ 'ohanin',
                                  # grepl('vQC', Genes) ~ 'vQC',
                                  grepl('SVSP', Genes) ~ 'SVSP',
                                  grepl('PLA2', Genes) ~ 'PLA2',
                                  grepl('CRISP', Genes) ~ 'CRISP',
                                  # grepl('CTL', Genes) ~ 'CTL',
                                  # grepl('EXO', Genes) ~ 'EXO',
                                  grepl('LAAO', Genes) ~ 'LAAO',
                                  grepl('myotoxin', Genes) ~ 'myotoxin',
                                  grepl('BPP', Genes) ~ 'BPP',
                                  TRUE ~ 'others')) %>% 
  mutate(Color = case_when(
    grepl('SVMP', Genes) ~ SVMP_color,
    # grepl('VEGF', Genes) ~ VEGF_color,
    # grepl('ohanin', Genes) ~ ohanin_color,
    # grepl('vQC', Genes) ~ vQC_color,
    grepl('SVSP', Genes) ~ SVSP_color,
    grepl('PLA2', Genes) ~ PLA2_color,
    grepl('CRISP', Genes) ~ CRISP_color,
    # grepl('CTL', Genes) ~ CTL_color,
    # grepl('EXO', Genes) ~ EXO_color,
    grepl('LAAO', Genes) ~ LAAO_color,
    grepl('myotoxin', Genes) ~ myotoxin_color,
    grepl('BPP', Genes) ~ BPP_color,
    TRUE ~ other_color
  ))



# Limit the data to CV0857
mRNA_CV0857_df <- mRNA_df %>% 
  dplyr::select(Venom.Family, CV0857, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Expression = mean(CV0857)) %>% 
  dplyr::select(-CV0857) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(mRNA_CV0857_df$Color, mRNA_CV0857_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- mRNA_CV0857_df$Venom.Family
value <- mRNA_CV0857_df$Expression

# Create pie chart for the mRNA data for CV0857
venom_mRNA_CV0857_pie_chart <- mRNA_CV0857_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
    ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Transcriptome in Northern ', italic('C. viridis'))))
venom_mRNA_CV0857_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Transcripome_viridis_CV0857_2024.06.23.pdf', venom_mRNA_CV0857_pie_chart, create.dir = T)



# Limit the data to CV1087 
mRNA_CV1087_df <- mRNA_df %>% 
  dplyr::select(Venom.Family, CV1087, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Expression = mean(CV1087)) %>% 
  dplyr::select(-CV1087) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(mRNA_CV1087_df$Color, mRNA_CV1087_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- mRNA_CV1087_df$Venom.Family
value <- mRNA_CV1087_df$Expression

# Create pie chart for the mRNA data for CV1087
venom_mRNA_CV1087_pie_chart <- mRNA_CV1087_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  # scale_fill_brewer(palette = 'Dark2') +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Transcriptome in Northern ', italic('C. viridis'))))
venom_mRNA_CV1087_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Transcripome_viridis_CV1087_2024.06.23.pdf', venom_mRNA_CV1087_pie_chart, create.dir = T)



# Limit the data to CV1081 
mRNA_CV1081_df <- mRNA_df %>% 
  dplyr::select(Venom.Family, CV1081, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Expression = mean(CV1081)) %>% 
  dplyr::select(-CV1081) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(mRNA_CV1081_df$Color, mRNA_CV1081_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- mRNA_CV1081_df$Venom.Family
value <- mRNA_CV1081_df$Expression

# Create pie chart for the mRNA data for CV1081
venom_mRNA_CV1081_pie_chart <- mRNA_CV1081_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Transcriptome in Central ', italic('C. viridis'))))
venom_mRNA_CV1081_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Transcripome_viridis_CV1081_2024.06.23.pdf', venom_mRNA_CV1081_pie_chart, create.dir = T)


# Limit the data to CV1086
mRNA_CV1086_df <- mRNA_df %>% 
  dplyr::select(Venom.Family, CV1086, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Expression = mean(CV1086)) %>% 
  dplyr::select(-CV1086) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(mRNA_CV1086_df$Color, mRNA_CV1086_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- mRNA_CV1086_df$Venom.Family
value <- mRNA_CV1086_df$Expression

# Create pie chart for the mRNA data for CV1086
venom_mRNA_CV1086_pie_chart <- mRNA_CV1086_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Transcriptome in Southern ', italic('C. viridis'))))
venom_mRNA_CV1086_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Transcripome_viridis_CV1086_2024.06.23.pdf', venom_mRNA_CV1086_pie_chart, create.dir = T)


# Limit the data to CV0987 
mRNA_CV0987_df <- mRNA_df %>% 
  dplyr::select(Venom.Family, CV0987, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Expression = mean(CV0987)) %>% 
  dplyr::select(-CV0987) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(mRNA_CV0987_df$Color, mRNA_CV0987_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- mRNA_CV0987_df$Venom.Family
value <- mRNA_CV0987_df$Expression

# Create pie chart for the mRNA data for CV0987
venom_mRNA_CV0987_pie_chart <- mRNA_CV0987_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Transcriptome in ', italic('C. lutosus'))))
venom_mRNA_CV0987_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Transcripome_lutosus_CV0987_2024.06.23.pdf', venom_mRNA_CV0987_pie_chart, create.dir = T)


# Limit the data to CV0985
mRNA_CV0985_df <- mRNA_df %>% 
  dplyr::select(Venom.Family, CV0985, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Expression = mean(CV0985)) %>% 
  dplyr::select(-CV0985) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(mRNA_CV0985_df$Color, mRNA_CV0985_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- mRNA_CV0985_df$Venom.Family
value <- mRNA_CV0985_df$Expression

# Create pie chart for the mRNA data for CV0985
venom_mRNA_CV0985_pie_chart <- mRNA_CV0985_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Transcriptome in ', italic('C. concolor'))))
venom_mRNA_CV0985_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Transcripome_concolor_CV0985_2024.06.23.pdf', venom_mRNA_CV0985_pie_chart, create.dir = T)




#### Venom Proteome Profile Pie Charts ####

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


# Set and read the data
protein_data <- 'Data/Protein/Formated_Protein_Data_IMPORTANT_2024.06.14.tsv'
protein_df <- read.table(file = protein_data, header = T) %>% 
  dplyr::rename('Genes' = 'Converted.Gene.IDs') %>%
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove, collapse = "|"), '')) %>% 
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove2, collapse = "|"), '')) %>% 
  filter(Feature == 'Intensity') %>% 
  filter(str_starts(Genes, 'Venom_')) %>% 
  dplyr::select(Intensity, Sample.ID, Genes) %>% 
  pivot_wider(
    names_from = Sample.ID,
    values_from = Intensity
  ) %>% 
  mutate(Venom.Family = case_when(grepl('SVMP', Genes) ~ 'SVMP', # Add a Venom.Families Column
                                  # grepl('VEGF', Genes) ~ 'VEGF',
                                  # grepl('ohanin', Genes) ~ 'ohanin',
                                  # grepl('vQC', Genes) ~ 'vQC',
                                  grepl('SVSP', Genes) ~ 'SVSP',
                                  grepl('PLA2', Genes) ~ 'PLA2',
                                  grepl('CRISP', Genes) ~ 'CRISP',
                                  # grepl('CTL', Genes) ~ 'CTL',
                                  # grepl('EXO', Genes) ~ 'EXO',
                                  grepl('LAAO', Genes) ~ 'LAAO',
                                  grepl('myotoxin', Genes) ~ 'myotoxin',
                                  grepl('BPP', Genes) ~ 'BPP',
                                  TRUE ~ 'others')) %>% 
  mutate(Color = case_when(
    grepl('SVMP', Genes) ~ SVMP_color,
    # grepl('VEGF', Genes) ~ VEGF_color,
    # grepl('ohanin', Genes) ~ ohanin_color,
    # grepl('vQC', Genes) ~ vQC_color,
    grepl('SVSP', Genes) ~ SVSP_color,
    grepl('PLA2', Genes) ~ PLA2_color,
    grepl('CRISP', Genes) ~ CRISP_color,
    # grepl('CTL', Genes) ~ CTL_color,
    # grepl('EXO', Genes) ~ EXO_color,
    grepl('LAAO', Genes) ~ LAAO_color,
    grepl('myotoxin', Genes) ~ myotoxin_color,
    grepl('BPP', Genes) ~ BPP_color,
    TRUE ~ other_color
  ))

# Create a new data frame for only CV0987
prot_CV0987_df <- protein_df %>% 
  dplyr::select(Venom.Family, CV0987, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Total.Intensity = mean(CV0987)) %>% 
  dplyr::select(-CV0987) %>% 
  # filter(Total.Intensity > 0) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(prot_CV0987_df$Color, prot_CV0987_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- prot_CV0987_df$Venom.Family
value <- prot_CV0987_df$Total.Intensity

# Create plot
venom_protein_CV0987_pie_chart <- prot_CV0987_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Proteome in ', italic('C. lutosus'))))
venom_protein_CV0987_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Proteome_lutosus_CV0987_2024.06.23.pdf', venom_protein_CV0987_pie_chart, create.dir = T)


# Create a new data frame for only CV0985
prot_CV0985_df <- protein_df %>% 
  dplyr::select(Venom.Family, CV0985, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Total.Intensity = mean(CV0985)) %>% 
  dplyr::select(-CV0985) %>% 
  # filter(Total.Intensity > 0) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(prot_CV0985_df$Color, prot_CV0985_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- prot_CV0985_df$Venom.Family
value <- prot_CV0985_df$Total.Intensity

# Create plot
venom_protein_CV0985_pie_chart <- prot_CV0985_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Proteome in ', italic('C. concolor'))))
venom_protein_CV0985_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Proteome_concolor_CV0985_2024.06.23.pdf', venom_protein_CV0985_pie_chart, create.dir = T)



# Create a new data frame for only CV1087
prot_CV1087_df <- protein_df %>% 
  dplyr::select(Venom.Family, CV1087, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Total.Intensity = mean(CV1087)) %>% 
  dplyr::select(-CV1087) %>% 
  # filter(Total.Intensity > 0) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(prot_CV1087_df$Color, prot_CV1087_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- prot_CV1087_df$Venom.Family
value <- prot_CV1087_df$Total.Intensity

# Create plot
venom_protein_CV1087_pie_chart <- prot_CV1087_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Proteome in Northern ', italic('C. viridis'))))
venom_protein_CV1087_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Proteome_viridis_CV1087_2024.06.23.pdf', venom_protein_CV1087_pie_chart, create.dir = T)



# Create a new data frame for only CV1086
prot_CV1086_df <- protein_df %>% 
  dplyr::select(Venom.Family, CV1086, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Total.Intensity = mean(CV1086)) %>% 
  dplyr::select(-CV1086) %>% 
  # filter(Total.Intensity > 0) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(prot_CV1086_df$Color, prot_CV1086_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- prot_CV1086_df$Venom.Family
value <- prot_CV1086_df$Total.Intensity

# Create plot
venom_protein_CV1086_pie_chart <- prot_CV1086_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Proteome in Southern ', italic('C. viridis'))))
venom_protein_CV1086_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Proteome_viridis_CV1086_2024.06.23.pdf', venom_protein_CV1086_pie_chart, create.dir = T)


# Create a new data frame for only CV1081
prot_CV1081_df <- protein_df %>% 
  dplyr::select(Venom.Family, CV1081, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Total.Intensity = mean(CV1081)) %>% 
  dplyr::select(-CV1081) %>% 
  # filter(Total.Intensity > 0) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(prot_CV1081_df$Color, prot_CV1081_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- prot_CV1081_df$Venom.Family
value <- prot_CV1081_df$Total.Intensity

# Create plot
venom_protein_CV1081_pie_chart <- prot_CV1081_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Proteome in Central ', italic('C. viridis'))))
venom_protein_CV1081_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Proteome_viridis_CV1081_2024.06.23.pdf', venom_protein_CV1081_pie_chart, create.dir = T)


# Create a new data frame for only CV0857
prot_CV0857_df <- protein_df %>% 
  dplyr::select(Venom.Family, CV0857, Color) %>% 
  group_by(Venom.Family) %>% 
  mutate(Total.Intensity = mean(CV0857)) %>% 
  dplyr::select(-CV0857) %>% 
  # filter(Total.Intensity > 0) %>% 
  distinct()

# Create a vector for colors
colors <- setNames(prot_CV0857_df$Color, prot_CV0857_df$Venom.Family)

# Set y variable and fill for the pie chart
venom_family <- prot_CV0857_df$Venom.Family
value <- prot_CV0857_df$Total.Intensity

# Create plot
venom_protein_CV0857_pie_chart <- prot_CV0857_df %>% 
  ggplot(aes(x = "", y = value, fill = venom_family)) +
  geom_bar(stat = 'identity', width = 1, color = 'black', linewidth = 0.25) +
  coord_polar(theta = 'y') +
  theme_void() +
  scale_fill_manual(values = colors) +
  ggtitle(expression(paste('Venom Proteome in Northern ', italic('C. viridis')))) +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +  # Adjust the margin and other title information
  labs(fill = 'Venom Family', title = expression(paste('Venom Proteome in Northern ', italic('C. viridis'))))
venom_protein_CV0857_pie_chart
ggsave('Figures/Pie_Charts/3UTR/Proteome_viridis_CV0857_2024.06.23.pdf', venom_protein_CV0857_pie_chart, create.dir = T)



#### miRNA Stacked Bar Graphs ####


# Create a data frame that can be used for stacked bars
miRNA_rpm_df <- mi_df %>% 
  dplyr::select(
    miRNA.Cluster, Genes, contains('miRNA.RPM.')
  ) %>% 
  pivot_longer(
    cols = contains('miRNA.RPM'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.RPM'
  ) %>% 
  mutate(Sample.ID = str_extract(Sample.ID, 'CV[0-9]+')) %>%
  mutate(miRNA.Cluster = gsub('cvi-', '', miRNA.Cluster)) %>% 
  # mutate(Sample.ID = gsub('miRNA.RPM.', '', Sample.ID)) %>% 
  distinct()

# Create graph for all targeted genes
miRNA_relative_frequencies <- ggplot(
  data = miRNA_rpm_df, aes(x = miRNA.Cluster, y = miRNA.RPM, fill = Sample.ID)
) +
  geom_bar(stat = 'identity', position = 'fill') + 
  labs(y = 'Relative Frequency of miRNAs per Sample', x = 'miRNA Sequence', fill = 'Sample ID') +
  scale_fill_brewer(palette = 'Spectral') +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 70, hjust = 1))
  coord_flip()
miRNA_relative_frequencies
# Save
ggsave('Figures/Bar_Graphs/3UTR/miRNA_relative_frequencies_all_genes_2024.06.23.pdf', plot = miRNA_relative_frequencies, width = 10, height = 15, create.dir = T)

# Create a graph for only venom gene targeting miRNAs
miRNA_relative_frequencies_venom <- miRNA_rpm_df %>% 
  filter(str_detect(Genes, 'Venom_|PLA2G2E.1')) %>% 
  ggplot(aes(x = miRNA.Cluster, y = miRNA.RPM, fill = Sample.ID)) +
    geom_bar(stat = 'identity', position = 'fill') + 
    labs(y = 'Relative Frequency of miRNAs per Sample', x = 'miRNA Sequence', fill = 'Sample ID') +
    scale_fill_brewer(palette = 'Spectral') +
    theme_classic() +
    # theme(axis.text.x = element_text(angle = 70, hjust = 1))
    coord_flip()
miRNA_relative_frequencies_venom
# Save
ggsave('Figures/Bar_Graphs/3UTR/miRNA_relative_frequencies_venom_genes_2024.06.23.pdf', plot = miRNA_relative_frequencies, width = 10, height = 15, create.dir = T)

# Create a bar graph of total miRNAs for each sample
miRNA_totals_plot <- miRNA_rpm_df %>%
  dplyr::select(-Genes) %>% 
  distinct() %>% # Use distinct so that the same miRNA isn't counted twice for the same gene  
  group_by(Sample.ID) %>%
  summarise(total_miRNA_RPM = sum(miRNA.RPM, na.rm = TRUE)) %>%
  mutate(
    Color = case_when(
      grepl('CV0857', Sample.ID) ~ 'CV0857',
      grepl('CV0985', Sample.ID) ~ 'CV0985',
      grepl('CV0987', Sample.ID) ~ 'CV0987',
      grepl('CV1081', Sample.ID) ~ 'CV1081',
      grepl('CV1086', Sample.ID) ~ 'CV1086',
      grepl('CV1087', Sample.ID) ~ 'CV1087',
    )
  ) %>%
  ggplot(aes(x = Sample.ID, y = total_miRNA_RPM, fill = Color)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Spectral') +
  labs(y = 'Total miRNA per sample (RPM)', x = 'Sample', title = 'Total miRNA per Sample') +
  theme_classic() +
  theme(legend.position = 'none')
miRNA_totals_plot
ggsave('Figures/Bar_Graphs/3UTR/miRNA_totals_per_sample_2024.06.23.pdf', plot = miRNA_totals_plot, create.dir = T)



#### Protein Levels Bar Graph ####

# Define path for data
protein_data <- 'Data/Protein/Formated_Protein_Data_IMPORTANT_2024.06.14.tsv'

# Read in the Protein expression data
protein_only_df <- read.table(file = protein_data, header = T) %>% 
  filter(Feature == 'Intensity') %>% 
  dplyr::select(
    Sample.ID, Converted.Gene.IDs, Intensity
  ) %>% 
  mutate(Sample.ID = str_extract(Sample.ID, 'CV[0-9]+')) %>% 
  group_by(Sample.ID) 

# Create a plot for all proteins
all_prot_plot <- protein_only_df %>% 
  summarise(Total.Protein = sum(Intensity, na.rm = T)) %>%
  mutate(
    Color = case_when(
      grepl('CV0857', Sample.ID) ~ 'CV0857',
      grepl('CV0985', Sample.ID) ~ 'CV0985',
      grepl('CV0987', Sample.ID) ~ 'CV0987',
      grepl('CV1081', Sample.ID) ~ 'CV1081',
      grepl('CV1086', Sample.ID) ~ 'CV1086',
      grepl('CV1087', Sample.ID) ~ 'CV1087',
    )
  ) %>%
  ggplot(aes(x = Sample.ID, y = Total.Protein, fill = Color)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Spectral') +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  labs(y = 'Total spectral intensity', x = 'Sample', title = 'Total Spectral Intensity per Sample') +
  theme_classic() + 
  theme(legend.position = 'none')
all_prot_plot  
ggsave('Figures/Bar_Graphs/3UTR/Protein_totals_per_sample_2024.06.23.pdf', plot = all_prot_plot, create.dir = T)


# Create plot for AGO protein levels only
ago_prot_plot <- protein_only_df %>% 
  filter(Converted.Gene.IDs == grepl('AGO', Converted.Gene.IDs))
# Doesn't seem that the AGO proteins were detected


#### mRNA Expression Levels Bar Graph ####

# Define path for the data
mrna_data <- 'Data/mRNA/RNAseq_Raw_Formated_IMPORTANT_2024.06.14.tsv'

# Read in the mRNA data
mrna_only_df <- read.table(file = mrna_data, header = T) %>% 
  dplyr::select(Sample.ID, Converted.Gene.IDs, RNA.Raw.Counts) %>% 
  mutate(Sample.ID = str_extract(Sample.ID, 'CV[0-9]+')) %>% 
  filter(!Sample.ID == 'CV1082') %>%  # Filter out the one without protein data
  group_by(Sample.ID)

# Create a plot for all mRNAs
all_mrna_plot <- mrna_only_df %>% 
  summarise(Total.mRNA = sum(RNA.Raw.Counts, na.rm = T)) %>% 
  mutate(
    Color = case_when(
      grepl('CV0857', Sample.ID) ~ 'CV0857',
      grepl('CV0985', Sample.ID) ~ 'CV0985',
      grepl('CV0987', Sample.ID) ~ 'CV0987',
      grepl('CV1081', Sample.ID) ~ 'CV1081',
      grepl('CV1086', Sample.ID) ~ 'CV1086',
      grepl('CV1087', Sample.ID) ~ 'CV1087',
    )
  ) %>%
  ggplot(aes(x = Sample.ID, y = Total.mRNA, fill = Color)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Spectral') +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  labs(y = 'Total mRNA Expression', x = 'Sample', title = 'Total mRNA Expression per Sample') +
  theme_classic() + 
  theme(legend.position = 'none')
all_mrna_plot
ggsave('Figures/Bar_Graphs/3UTR/mRNA_totals_per_sample_2024.06.23.pdf', plot = all_mrna_plot, create.dir = T)



#### mRNA Totals vs Protein Totals vs miRNA Correlations ####

# Create a summ of the protein data
protein_sum_df <- protein_only_df %>% 
  summarise(Total.Protein = sum(Intensity, na.rm = T))

# Create a summed mRNA data frame
mrna_sum_df <- mrna_only_df %>% 
  summarise(Total.mRNA = sum(RNA.Raw.Counts, na.rm = T))

# Create a summed miRNA data frame
mirna_sum_df <-  miRNA_rpm_df %>% dplyr::select(-Genes) %>% 
  distinct() %>% # Use distinct so that the same miRNA isn't counted twice for the same gene
  group_by(Sample.ID) %>%
  summarise(Total.miRNA.RPM = sum(miRNA.RPM, na.rm = TRUE))

# Protein and mRNA joined data frame
mRNAvsProtein_df <- left_join(
  mrna_sum_df,
  protein_sum_df,
  by = c('Sample.ID')
) 

# Protein, mRNA, and miRNA joined data frame
mRNA_vs_Protein_vs_miRNA_df <- left_join(
  mRNAvsProtein_df,
  mirna_sum_df,
  by = c('Sample.ID')
) %>% 
  mutate(
    Color = case_when(
      grepl('CV0857', Sample.ID) ~ 'CV0857',
      grepl('CV0985', Sample.ID) ~ 'CV0985',
      grepl('CV0987', Sample.ID) ~ 'CV0987',
      grepl('CV1081', Sample.ID) ~ 'CV1081',
      grepl('CV1086', Sample.ID) ~ 'CV1086',
      grepl('CV1087', Sample.ID) ~ 'CV1087',
    )
  )

# Create plot to correlate mRNA and Protein expression levels
mRNA_vs_protein_sum_plot <- ggplot(data = mRNA_vs_Protein_vs_miRNA_df, aes(x = Total.mRNA, y = Total.Protein)) +
  geom_point(aes(color = Sample.ID), size = 3, alpha = 0.8) +
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
  scale_fill_brewer(palette = 'Spectral') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  labs(y = 'Total Protein Expression per Sample', x = 'Total mRNA Expression per Sample', title = 'Total mRNA vs Protein Expression')
mRNA_vs_protein_sum_plot
ggsave('Figures/Expression_Plots/3UTR/Total_Sample_Regressions/mRNA_vs_Protein_Expression_Totals.2024.06.23.pdf', plot = mRNA_vs_protein_sum_plot, create.dir = T)


# Create plot to correlate miRNA and Protein expression levels
miRNA_vs_protein_sum_plot <- ggplot(data = mRNA_vs_Protein_vs_miRNA_df, aes(x = Total.miRNA.RPM, y = Total.Protein)) +
  geom_point(aes(color = Sample.ID), size = 3, alpha = 0.8) +
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
  scale_fill_brewer(palette = 'Spectral') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  labs(y = 'Total Protein Expression per Sample', x = 'Total miRNA Expression per Sample', title = 'Total miRNA vs Protein Expression')
miRNA_vs_protein_sum_plot
ggsave('Figures/Expression_Plots/3UTR/Total_Sample_Regressions/miRNA_vs_Protein_Expression_Totals.2024.06.23.pdf', plot = mRNA_vs_protein_sum_plot, create.dir = T)



#### Linear Model For mRNA-Protein-miRNA Totals and Residuals against miRNA graph ####

# Set variables for protein-mRNA linear model
total_mRNA_expression <- mRNA_vs_Protein_vs_miRNA_df$Total.mRNA
total_protein_expression <- mRNA_vs_Protein_vs_miRNA_df$Total.Protein
total_miRNA_expression <- mRNA_vs_Protein_vs_miRNA_df$Total.miRNA.RPM
sample_color <- mRNA_vs_Protein_vs_miRNA_df$Color
  
  
# Linear model
total_expression_lm <- lm(Total.Protein ~ Total.mRNA, data = mRNA_vs_Protein_vs_miRNA_df)
total_expression_slope <- coef(total_expression_lm)[2] 
total_expression_intercept <- coef(total_expression_lm)[1]
summary(total_expression_lm)

# Get sample residuals
totals_residuals_df <- mRNA_vs_Protein_vs_miRNA_df %>% mutate(Residuals = resid(total_expression_lm))

# Graph of residuals of total mRNA vs Total protein vs miRNA total expression per sample
totals_residuals_plot <- ggplot(data = totals_residuals_df, mapping = aes(x = total_miRNA_expression, y = Residuals)) +
  geom_point(aes(color = Sample.ID), size = 3, alpha = 0.8) +
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
  scale_fill_brewer(palette = 'Spectral') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  ) +
  labs(y = 'Residuals', x = 'Total miRNA Expression per Sample', title = 'Residuals of Linear Regression Against Total miRNA Expression Levels per Sample')
totals_residuals_plot
ggsave('Figures/Expression_Plots/3UTR/Total_Sample_Regressions/miRNA_vs_Residuals_of_mRNAvsProtein_Linear_Regression.2024.06.23.pdf', plot = totals_residuals_plot, create.dir = T)



#### miRNA Binding Energy histogram Plot ####

# Create data frame for a histogram of miRNA origins
mi_histogram_df <- read.table(file = miRNA_mRNA_protein_data, header = T) %>% 
  dplyr::select(miRNA.Cluster, Max.Energy, Origin) %>% 
  distinct()

# Create plot
mirna_histogram_plot <- ggplot(data = mi_histogram_df, mapping = aes(x = Max.Energy, fill = Origin)) +
  geom_histogram(alpha = 0.8)
mirna_histogram_plot



#### Gene family plots dataframe for Venom Gene Dot Plots (VST) ####

# Create a data frame for venom genes only
venom_genes_vst_df <- mi_df %>%
  filter(str_starts(Genes, 'Venom_')) %>%
  dplyr::select(
    miRNA.Cluster,
    Genes,
    contains('miRNA.Counts.'),
    contains('RNA.VST.'),
    contains('Intensity.'),
    contains('miRNA.RPM.')
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


#### mRNA (VST) vs Protein for all samples ####


# Pivot data frame so that all samples are in a single column
all_samples_plot_protein_df <- venom_genes_vst_df %>%
  dplyr::select(-contains('miRNA.'), -contains('Intensity')) %>%
  distinct() %>%
  pivot_longer(
    cols = contains('RNA.VST.'),
    names_to = 'Sample.ID',
    values_to = 'RNA.VST'
  ) %>%
  mutate(Sample.ID = gsub('RNA.VST.', '', Sample.ID))

# Pivot this one too
all_samples_plot_mRNA_df <- venom_genes_vst_df %>%
  dplyr::select(-contains('miRNA.'), -contains('RNA.')) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('Intensity.'),
    names_to = 'Sample.ID',
    values_to = 'Intensity'
  ) %>% 
  mutate(Sample.ID = gsub('Intensity.', '', Sample.ID))

# Fuse back together
all_samples_plot_df <- left_join(
  all_samples_plot_protein_df,
  all_samples_plot_mRNA_df,
  by = c('Sample.ID', 'Genes', 'Venom.Family')
) 

# write.table(all_samples_plot_df, file = '/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Fxn_NSF/_Manuscripts/3_Venom_ncRNA/Data/Normalized_RNA_Protein.tsv', sep = '\t')


# Calculate CLR
all_samples_clr_df <- all_samples_plot_df %>%
  group_by(Sample.ID) %>%
  mutate(Scaled.Protein.Expression = Intensity / sum(Intensity)) %>% # Sum protein expression to 1
  mutate(Scaled.mRNA.Expression = RNA.VST / sum(RNA.VST)) %>% # Sum mRNA expression to 1
  ungroup() %>%
  mutate(CLR.Scaled.Protein.Expression = as.numeric(compositions::clr(Scaled.Protein.Expression))) %>%
  mutate(CLR.Scaled.mRNA.Expression = as.numeric(compositions::clr(Scaled.mRNA.Expression)))


# Set variables for protein and RNA expression levels
all_mRNA_expression <- all_samples_clr_df$CLR.Scaled.mRNA.Expression
all_protein_expression <- all_samples_clr_df$CLR.Scaled.Protein.Expression
venom_family <- all_samples_clr_df$Venom.Family
# gene_label <- all_samples_clr_df$Genes


# Create color scheme for the venom genes
venom_colors <- c(
  SVMP = '#4A70B5',
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


# Find linear regression equation
venom_lm <- lm(all_protein_expression ~ all_mRNA_expression - 1, data = all_samples_clr_df)
venom_slope <- coef(venom_lm)[2]
venom_intercept <- coef(venom_lm)[1]

# Create figure for mRNA and Protein expression
all_gene_expression_plot <- all_samples_clr_df %>% 
  ggplot(aes(x = all_mRNA_expression, y = all_protein_expression)) +
  # ggplot(aes(x = all_mRNA_expression, y = all_protein_expression, label = gene_label)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x - 1
  ) + 
  # geom_text(check_overlap = T, size = 1, hjust = -0.2) +
  #xlim(0, NA) +  # Set the desired x-axis limits
  #ylim(0, NA) +  # Set y-axis limits to start at 0 and extend to the maximum value
  geom_abline(slope = 1, color = 'black', linetype = 'solid') +
  scale_x_continuous() +  # Set x-axis to continuous scale
  scale_y_continuous() +  # Set y-axis to continuous scale
  xlab('clr(gene expression)') +
  ylab('clr(peak intensity)') +
  scale_color_manual(values = venom_colors) +  # Apply color scheme  
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.title = element_blank())
all_gene_expression_plot
# It seems the data, at least if you plot it for all of the individuals is Heteroskedastic and violates the assumptions of the linear model

# Save 
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/All_Individuals_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = all_gene_expression_plot, create.dir = T)



#### Dataframe for residuals and miRNA (VST) ####

# Find linear regression equation
venom_lm <- lm(all_protein_expression ~ all_mRNA_expression - 1, data = all_samples_clr_df)
venom_slope <- coef(venom_lm)[2]
venom_intercept <- coef(venom_lm)[1]

# Find the predicted values from the model and put them into the dataframe
all_samples_residuals_df <- all_samples_clr_df %>% 
  mutate(Residuals = resid(venom_lm))



# Create a dataframe that has miRNAs and all the protein information in it that I can fuse later
miRNA_protein_pivot_df <- mi_df %>%
  filter(str_starts(Genes, 'Venom_')) %>%
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
  dplyr::select(-contains('miRNA.Counts.'), -contains('RNA.VST')) %>% 
  distinct() %>% 
  pivot_longer(
    cols = contains('Intensity.'),
    names_to = 'Sample.ID',
    values_to = 'Intensity'
  ) %>% 
  mutate(Sample.ID = gsub('Intensity.', '', Sample.ID)) %>% 
  distinct()

dup_rows1 <- miRNA_protein_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, Sample.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1)


# Create dataframe that has miRNAs in it that I can fuse to this one.
miRNA_mRNA_pivot_df <- mi_df %>% 
  filter(str_starts(Genes, 'Venom_')) %>%
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
  dplyr::select(-contains('miRNA.Counts.'), -contains('Intensity')) %>% 
  pivot_longer(
    cols = contains('RNA.VST'),
    names_to = 'Sample.ID',
    values_to = 'RNA.VST'
  ) %>% 
  mutate(Sample.ID = gsub('RNA.VST.', '', Sample.ID)) %>% 
  distinct()

# Identify duplicates
dup_rows2 <- miRNA_mRNA_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Start, miRNA.End, miRNA.Target.Start, Sample.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1)

# Create a variable to share common names between columns in the protein and mRNA data sets
shared_columns = intersect(names(miRNA_protein_pivot_df), names(miRNA_mRNA_pivot_df))

# Create a dataframe by fusing the protein and mRNA together
mRNA_protein_pivot_df <- left_join(
  miRNA_protein_pivot_df,
  miRNA_mRNA_pivot_df,
  by = shared_columns
) %>% 
  distinct()

# Identify duplicates
dup_rows3 <- mRNA_protein_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, Sample.ID) %>%
  summarise(n = n()) %>%
  filter(n > 1)



# Create dataframe of only miRNAs and genes with their samples id
miRNA_pivot_df <- mi_df %>% 
  filter(str_starts(Genes, 'Venom_')) %>%
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
  dplyr::select(-contains('miRNA.Counts.'), -contains('Intensity'), -contains('RNA.VST')) %>% 
  pivot_longer(
    cols = contains('miRNA.RPM'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.RPM'
  ) %>%
  mutate(Sample.ID = gsub('miRNA.RPM.', '', Sample.ID)) %>% 
  distinct()

# Check for duplicated rows  
dup_rows4 <- miRNA_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, miRNA.Target.End, Venom.Family) %>%
  summarise(n = n()) %>%
  filter(n > 1)

# Get common columns between the two dataframes so they can be joined efficiently.
common_columns <- intersect(names(mRNA_protein_pivot_df), names(miRNA_pivot_df))

# Fuse dataframes back together
miRNA_mRNA_protein_pivot_df <- left_join(
  mRNA_protein_pivot_df,
  miRNA_pivot_df,
  by = common_columns
) %>% 
  dplyr::select(-contains('miRNA.RPM.')) %>% 
  distinct()

# Check for duplicated rows  
dup_rows5 <- miRNA_mRNA_protein_pivot_df %>%
  group_by(Genes, miRNA.Cluster, miRNA.Cluster.Original, miRNA.Start, miRNA.End, miRNA.Target.Start, miRNA.Target.End, Venom.Family, Intensity, RNA.VST) %>%
  summarise(n = n()) %>%
  filter(n > 1)


# Create a set of shared columns
columns_shared <- intersect(names(miRNA_mRNA_protein_pivot_df), names(all_samples_residuals_df))

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

# IMPORTANT: DON'T DELETE
# Fuse to the residuals dataframe
miRNA_with_residuals_df <- left_join(
  all_samples_residuals_df,
  miRNA_mRNA_protein_pivot_df,
  by = columns_shared
) %>% 
  distinct() %>% 
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove, collapse = "|"), '')) %>% 
  mutate(Sample.ID = str_replace(Sample.ID, paste(characters_to_remove2, collapse = "|"), ''))

# IMPORTANT: DON'T DELETE
# Create a wider version of the above
wider_miRNA_with_residuals_df <- miRNA_with_residuals_df %>% 
  pivot_wider(
    names_from = 'miRNA.Cluster',
    values_from = 'miRNA.RPM'
  ) %>% 
  dplyr::select(
    -contains('miRNA'), -contains('Blast'), -contains('Score'), -contains('Energy'), -Positions, -Strand, -E.value
  )




#### Residuals Graphs (based on VST) ####
# 
# # Create a dataframe just for ohanin and filter at any missing data
# ohanin_df <- wider_miRNA_with_residuals_df %>% 
#   filter(Genes == 'Venom_ohanin') %>% 
#   dplyr::select(where(~any(!is.na(.))))
# 
# # Create plot for cvi-miR-145-5p
# # Create variable to control x and y of the ggplot
# miRNA_expression <- as.numeric(ohanin_df$Cluster_2576)
# exp_residuals <- abs(as.numeric(ohanin_df$Residuals))
# gene <- ohanin_df$Genes
# 
# # Create graph of values
# miR_145_5p_plot <- ohanin_df %>% 
#   ggplot(aes(x = miRNA_expression, y = exp_residuals)) +
#   geom_point(aes(color = gene), size = 2.5, alpha = 0.8) +
#   ggtitle('cvi-miR-145-5p miRNA Abundance compared to residuals') +
#   geom_smooth(
#     method = 'lm',
#     se = T,
#     color = 'black',
#     linetype = 'dashed',
#     formula = y ~ x - 1
#   ) +
#   stat_poly_eq(
#     aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
#     formula = y ~ x - 1
#   ) + 
#   xlab('miRNA Abundance (RPM)') +
#   ylab('| Residuals |') +
#   scale_color_manual(values = venom_colors) +
#   theme_classic() +
#   theme(legend.position = 'bottom',
#         legend.title = element_blank())
# miR_145_5p_plot
# 
# 

source('Scripts/R/Protein-miRNA_RNA_Analysis/Residuals_Correlation_Function.R')

# Create a list of venom genes for a function to iterate through
venom_genes <- wider_miRNA_with_residuals_df$Genes
# Create path for the plots to go in
path <- 'Figures/Expression_Plots/3UTR/Linear_Regression_VST/Individual_miRNA_Gene_Relationships'

# Loop the function 
for (gene in venom_genes) {
  # Run the function in a loop
  plots <- residuals_vs_miRNA_plot2(wider_miRNA_with_residuals_df, gene, path, '2024.07.30', filter_r_squared = F)
}


# Add only the plots with R squared higher than 
# Create a second path
path2 <- 'Figures/Expression_Plots/3UTR/Linear_Regression_VST/Individual_miRNA_Gene_Relationships/Good_R_squared'

# Loop through the function
for (gene in venom_genes) {
  # Run the function in a loop
  plots <- residuals_vs_miRNA_plot2(wider_miRNA_with_residuals_df, gene, path2, '2024.07.30', filter_r_squared = T)
}



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
  ggplot(aes(x = viridis_mRNA_expression, y = viridis_protein_expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x -1
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
viridis_only_plot

# Save 
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/Viridis_Only_Individuals_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = viridis_only_plot, create.dir = T)



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
  ggplot(aes(x = mRNA_expression, y = protein_expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x -1
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
CV1081_plot

# Save 
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/CV1081_viridis_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = CV1081_plot, create.dir = T)


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
  ggplot(aes(x = mRNA_expression, y = protein_expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x -1
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
CV0857_plot

# Save 
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/CV0857_viridis_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = CV0857_plot, create.dir = T)



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
  ggplot(aes(x = mRNA_expression, y = protein_expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x -1
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
CV1086_plot

# Save 
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/CV1086_viridis_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = CV1086_plot, create.dir = T)



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
  ggplot(aes(x = mRNA_expression, y = protein_expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x -1
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
CV1087_plot

# Save 
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/CV1087_viridis_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = CV1087_plot, create.dir = T)


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
  ggplot(aes(x = mRNA_expression, y = protein_expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x -1
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
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/CV0987_lutosus_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = CV0987_plot, create.dir = T)



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
  ggplot(aes(x = mRNA_expression, y = protein_expression)) +
  geom_point(aes(color = venom_family), size = 2.5, alpha = 0.8) +
  geom_smooth(
    method = 'lm', 
    se = T, 
    color = 'black', 
    linetype = 'dashed', 
    formula = y ~ x - 1
  ) +
  stat_poly_eq(
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
    formula = y ~ x -1
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
CV0985_plot

# Save 
ggsave('Figures/Expression_Plots/3UTR/Linear_Regression_VST/CV0985_concolor_mRNAvsProtein_Expression_Plot_2024.06.23.pdf', plot = CV0985_plot, create.dir = T)


### mRNA Heat Maps ###

#### Prepare general data frame for Heat maps ####

# # I am going to exclude they don't have very high expression
# low_transcripts = c(
#   'AHSG',
#   'AHSG.1',
#   'AHSG.2',
#   'ALB',
#   'PINLYP.1',
#   'PINLYP.10'
# )

# Create dataframe
mRNA_all_genes_df <- mi_df %>% 
  dplyr::rename(
    'C. viridis - CV1081' = 'RNA.VST.LVG.2.CV1081.viridis.Mid.M',
    'C. viridis - CV0857' = 'RNA.VST.LVG.4.CV0857.viridis.North.M',
    'C. viridis - CV1086' = 'RNA.VST.LVG.9.CV1086.viridis.South.M',
    'C. viridis - CV1087' = 'RNA.VST.RVG.5S.CV1087.viridis.North.F',
    'C. o. lutosus - CV0987' = 'RNA.VST.RVG.6S.CV0987.lutosus.Other.M',
    'C. o. concolor - CV0985' = 'RNA.VST.RVG.7S.CV0985.concolor.Other.F'
  ) %>% 
  #filter(!(Genes %in% low_transcripts)) %>% 
  dplyr::select('Genes', 'C. viridis - CV1081', 'C. viridis - CV0857', 'C. viridis - CV1086', 'C. viridis - CV1087', 'C. o. lutosus - CV0987', 'C. o. concolor - CV0985') %>% 
  distinct() %>% 
  pivot_longer(cols = -Genes, names_to = 'Individuals', values_to = 'mRNA.Expression') %>%  # Pivot the information for heatmap
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
  group_by(Venom.Family)

# # Calculate variance in mRNA
# mRNA_variance_df <- mRNA_all_genes_df %>% 
#   group_by(Genes) %>% 
#   summarize(variance = var(mRNA.Expression))

# # Reorder Genes based on variance
# var_sorted_genes <- mRNA_variance_df %>% 
#   arrange((desc(variance))) %>% # Rearder Genes based on variance
#   pull(Genes)

# # Reorder Genes in the heat_df
# mRNA_all_genes_df$Genes <- factor(mRNA_all_genes_df$Genes, levels = var_sorted_genes)


#### mRNA Heat Map Everything ####


# For these heat maps make sure to add a variance heat map over it 


# Create heat map
mRNA_all_heatmap <- ggplot(mRNA_all_genes_df, aes(y = Individuals, x = Genes, fill = mRNA.Expression)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(
    y = 'Individuals', 
    x = 'Genes', 
    fill = 'mRNA Expression', 
    title = "mRNA Expression Level for all Targeted Genes"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
  
# Print
mRNA_all_heatmap
# Save
ggsave("Figures/Heat_Maps/3UTR/mRNA_all_genes_3UTR_heatmap_2024.06.25.pdf", plot = mRNA_all_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)


# Now for only venom
# Prepare data
all_venom_genes_mRNA_df <- mRNA_all_genes_df %>% 
  filter(str_starts(Genes, "Venom")) %>% 
  distinct()

# Create heat map
all_venom_genes_mRNA_heatmap <- ggplot(all_venom_genes_mRNA_df, aes(y = Individuals, x = Genes, fill = mRNA.Expression)) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(y = 'Individuals', x = 'Genes', fill = 'mRNA Expression', title = "mRNA Expression Level for all Targeted Venom Genes") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
# Print
all_venom_genes_mRNA_heatmap
# Save
ggsave("Figures/Heat_Maps/3UTR/venom_mRNA_3UTR_heatmap_2024.06.25.pdf", plot = all_venom_genes_mRNA_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)


#### miRNA Cluster Heat Map Everything ####

# Prepare data for heat-mapification
all_miRNA_df <- mi_df %>% 
  dplyr::rename(
    'C. viridis - CV1081' = 'miRNA.RPM.LVG.2.CV1081.viridis.Mid.M',
    'C. viridis - CV0857' = 'miRNA.RPM.LVG.4.CV0857.viridis.North.M',
    'C. viridis - CV1086' = 'miRNA.RPM.LVG.9.CV1086.viridis.South.M',
    'C. viridis - CV1087' = 'miRNA.RPM.RVG.5S.CV1087.viridis.North.F',
    'C. o. lutosus - CV0987' = 'miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M',
    'C. o. concolor - CV0985' = 'miRNA.RPM.RVG.7S.CV0985.concolor.Other.F'
  ) %>% 
  #filter(!(Genes %in% low_miRNAs)) %>% 
  dplyr::select('Genes', 'miRNA.Cluster', 'C. viridis - CV1081', 'C. viridis - CV0857', 'C. viridis - CV1086', 'C. viridis - CV1087', 'C. o. lutosus - CV0987', 'C. o. concolor - CV0985') %>% 
  pivot_longer(cols = c(-miRNA.Cluster, -Genes), names_to = 'Individuals', values_to = 'miRNA.Expression') %>%  # Pivot the information for heatmap
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
  group_by(Venom.Family)

# # Calculate total expression in miRNA
# miRNA_total_expression_df <- all_miRNA_df %>%
#   group_by(miRNA.Cluster) %>%
#   summarize(total_expression = sum(miRNA.Expression))
# 
# # Reorder miRNA cluster based on total expression
# total_sorted_genes <- miRNA_total_expression_df %>%
#   arrange(dplyr::desc(total_expression)) %>%  
#   pull(miRNA.Cluster)
# 
# # Reorder miRNA cluster in the df
# all_miRNA_df$miRNA.Cluster <- factor(all_miRNA_df$miRNA.Cluster, levels = total_sorted_genes)


# Calculate variance in miRNA
miRNA_variance_df <- all_miRNA_df %>%
  group_by(miRNA.Cluster) %>%
  summarize(variance = var(miRNA.Expression))

# Reorder miRNA cluster based on variance
var_sorted_genes <- miRNA_variance_df %>%
  arrange((dplyr::desc(variance))) %>% # Rearder miRNA cluster based on variance
  pull(miRNA.Cluster)

# Reorder miRNA cluster in the df
all_miRNA_df$miRNA.Cluster <- factor(all_miRNA_df$miRNA.Cluster, levels = var_sorted_genes)

# Drop Genes column so that it doesn't overwrite the clusters multiple times
# Apply this to the heat maps
all_miRNA_df <- all_miRNA_df %>% dplyr::select(-Genes) %>% distinct()


# Create heat map
miRNA_all_heatmap <- ggplot(all_miRNA_df, aes(y = Individuals, x = miRNA.Cluster, fill = log(miRNA.Expression + 1))) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(y = 'Individuals', x = 'miRNA Cluster', fill = 'miRNA Expression', title = "miRNA Expression Level for all miRNA Clusters (Log Scaled)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
miRNA_all_heatmap
# Save
ggsave("Figures/Heat_Maps/3UTR/miRNA_all_3UTR_heatmap_2024.06.25.pdf", plot = miRNA_all_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)

# Now for only venom
# Prepare data
venom_miRNA_df <- all_miRNA_df %>% 
  filter(!Venom.Family == 'others')

# Create heat map
venom_miRNA_heatmap <- ggplot(venom_miRNA_df, aes(y = Individuals, x = miRNA.Cluster, fill = log(miRNA.Expression + 1))) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(y = 'Individuals', x = 'miRNA Cluster', fill = 'miRNA Expression', title = "miRNA Expression Level for all Venom targeting miRNA Clusters (Log Scale)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
venom_miRNA_heatmap
# Save
ggsave("Figures/Heat_Maps/3UTR/venom_miRNA_3UTR_heatmap_2024.06.25.pdf", plot = venom_miRNA_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)



#### Protein Heat Map Everything ####

# Proteins to be removed because their expression is low accross the board
# low_proteins <- c(
#   'CCN1.1',
#   'KRT17.1',
#   'OS9'
# )


# Prepare data for heat-mapification
protein_all_df <- mi_df %>% 
  dplyr::rename(
    'C. viridis - CV1081' = 'Intensity.LVG.2.CV1081.viridis.Mid.M',
    'C. viridis - CV0857' = 'Intensity.LVG.4.CV0857.viridis.North.M',
    'C. viridis - CV1086' = 'Intensity.LVG.9.CV1086.viridis.South.M',
    'C. viridis - CV1087' = 'Intensity.RVG.5S.CV1087.viridis.North.F',
    'C. o. lutosus - CV0987' = 'Intensity.RVG.6S.CV0987.lutosus.Other.M',
    'C. o. concolor - CV0985' = 'Intensity.RVG.7S.CV0985.concolor.Other.F'
  ) %>% 
  #filter(!(Genes %in% low_proteins)) %>% 
  dplyr::select('Genes', 'C. viridis - CV1081', 'C. viridis - CV0857', 'C. viridis - CV1086', 'C. viridis - CV1087', 'C. o. lutosus - CV0987', 'C. o. concolor - CV0985') %>% 
  distinct() %>% 
  pivot_longer(cols = c(-Genes), names_to = 'Individuals', values_to = 'Intensity') %>%  # Pivot the information for heatmap
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
  group_by(Venom.Family)

# # Calculate variance in miRNA
# protein_variance_df <- protein_all_df %>% 
#   group_by(Genes) %>% 
#   summarize(variance = var(Intensity))

# # Reorder miRNA cluster based on variance
# var_sorted_genes <- protein_variance_df %>% 
#   arrange((desc(variance))) %>% # Rearder miRNA cluster based on variance
#   pull(Genes)

# # Reorder miRNA cluster in the df
# protein_all_df$Genes <- factor(protein_all_df$Genes, levels = var_sorted_genes)


# Create heat map
protein_all_heatmap <- ggplot(protein_all_df, aes(y = Individuals, x = Genes, fill = log(Intensity + 1))) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(y = 'Individuals', x = 'Proteins', fill = 'Protein Expression', title = "Protein Expression Level for all Targeted Genes (log scaled)") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
protein_all_heatmap
# Save
ggsave("Figures//Heat_Maps/3UTR/protein_all_3UTR_heatmap_2024.06.25.pdf", plot = protein_all_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)

# Now for only venom
# Prepare data
venom_protein_df <- protein_all_df %>% 
  filter(str_starts(Genes, "Venom")) %>%
  distinct()

# Create heat map
venom_protein_heatmap <- ggplot(venom_protein_df, aes(y = Individuals, x = Genes, fill = log(Intensity + 1))) +
  geom_tile() +
  scale_fill_viridis_c(option = 'magma') +
  labs(y = 'Individuals', x = 'Proteins', fill = 'Protein Expression', title = "Protein Expression Level for all Targeted Venom Genes") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
venom_protein_heatmap
# Save
ggsave("Figures//Heat_Maps/3UTR/venom_protein_3UTR_heatmap_2024.06.25.pdf", plot = venom_protein_heatmap, width = 18, height = 4, dpi = 900, create.dir = T)






#### Create General Dataframe for CIRCOS ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
source('Scripts/R/Protein-miRNA_RNA_Analysis/CIRCOS_Functions.R')

# Create smaller data frame for circos plot circos.link function to draw lines
venom_circos_df <- mi_df %>% 
  dplyr::select(
    miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, Positions,
    contains('miRNA.RPM.'),
    contains('RNA.'),
    contains('Intensity.')
  ) %>% 
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  filter(!miRNA.Target.Chrom %in% c('PE-reconstructed-10x-myo')) %>% 
  filter(!str_detect(miRNA.Sequence.Chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  filter(!str_detect(miRNA.Target.Chrom, 'scaffold-un')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  mutate(Positions = strsplit(as.character(Positions), ' ')) %>% # Separate the Positions column into multiple rows
  unnest(Positions) %>% 
  mutate(Positions = as.numeric(Positions)) %>%  # Convert Positions to numeric
  filter(!is.na(Positions)) %>% 
  # Create the new column miRNA.Target.Position
  mutate(miRNA.Target.Position = miRNA.Target.Start + Positions) %>% 
  mutate(Color = case_when(
    str_starts(Genes, 'Venom_PLA2') ~ PLA2_color,
    str_starts(Genes, 'Venom_SVSP') ~ SVSP_color,
    str_starts(Genes, 'Venom_SVMP') ~ SVMP_color,
    str_starts(Genes, 'Venom_vQC') ~ vQC_color,
    str_starts(Genes, 'Venom_ohanin') ~ ohanin_color,
    str_starts(Genes, 'Venom_LAAO') ~ LAAO_color,
    str_starts(Genes, 'Venom_CTL') ~ CTL_color,
    str_starts(Genes, 'Venom_CRISP') ~ CRISP_color,
    str_starts(Genes, 'Venom_BPP') ~ BPP_color,
    str_starts(Genes, 'Venom_VEGF') ~ VEGF_color,
    str_starts(Genes, 'Venom_myotoxin') ~ myotoxin_color,
    str_starts(Genes, 'Venom_EXO') ~ EXO_color,
    T ~ other_color
  ))

# Set track height
trackheight <- 0.2

#### CIRCOS Plot for RVG_5S colored by venom gene family ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)
source('Scripts/R/Protein-miRNA_RNA_Analysis/CIRCOS_Functions.R')


# Create smaller data frame for circos plot circos.link function to draw lines
# This data frame also only has RVG_5S information in it and is color coded for venom gene families
RVG5_venom_circos_df <- venom_circos_df %>% 
  dplyr::select(
    -contains('.RVG.6S.'),
    -contains('.RVG.7S.'),
    -contains('.LVG.2.'),
    -contains('.LVG.4.'),
    -contains('.LVG.9.')
  ) %>% 
  filter(!miRNA.RPM.RVG.5S.CV1087.viridis.North.F == 0) %>% # filter out any miRNA with zero expression, I don't think there are any, but just to be safe
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  distinct()

# Initialize max count values
max_protein_value <- range(log(RVG5_venom_circos_df$Intensity.RVG.5S.CV1087.viridis.North.F + 1), na.rm = T)
max_miRNA_count <- range(log(RVG5_venom_circos_df$miRNA.RPM.RVG.5S.CV1087.viridis.North.F + 1), na.rm = T)
max_mRNA_expression <- range(log(RVG5_venom_circos_df$RNA.VST.RVG.5S.CV1087.viridis.North.F + 1), na.rm = T)

# Define breaks for color scales
protein_breaks <- quantile(log(RVG5_venom_circos_df$Intensity.RVG.5S.CV1087.viridis.North.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
miRNA_breaks <- quantile(log(RVG5_venom_circos_df$miRNA.RPM.RVG.5S.CV1087.viridis.North.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
mRNA_breaks <- quantile(log(RVG5_venom_circos_df$RNA.VST.RVG.5S.CV1087.viridis.North.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)

# Create color ramps using magma color map
protein_color_ramp <- colorRamp2(protein_breaks, magma(3))
miRNA_color_ramp <- colorRamp2(miRNA_breaks, magma(3))
mRNA_color_ramp <- colorRamp2(mRNA_breaks, magma(3))



# Load scaffold sizes
# This will be used to initiate the sizes of each sector of the
# circos plot
scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
  dplyr::select(-V4, -V5)

# rename columns in scaffold_size
names(scaffold_size) = c('Chrom','size','genome_position')

# create a column of zeros, to indicate the starting
# value of each sector
scaffold_size = scaffold_size %>% mutate(start = 0)

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/3UTR/Venom_Genes_Colored_CIRCOS_CV1087_viridis_2024.06.24.pdf', width = 10, height = 10)

circos.clear()

# Set basic parameters for the circos plot
circos.par('track.height' = trackheight,
           cell.padding = c(0.01, 0, 0.01, 0),
           start.degree = 90,
           gap.degree = c(rep(1, 17), 20))

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))

# # Draw circos tracks for protein
# draw_circos_tracks_points2(
#   RVG5_venom_circos_df,
#   max_value = max_protein_value,
#   track_height = 0.2,
#   chromosome_col = 'miRNA.Target.Chrom',
#   loci_location_col =  'miRNA.Target.Start',
#   value_col = 'Intensity.RVG.5S.CV1087.viridis.North.F',
#   color = 'purple'
#   )

# Draw circos tracks for mRNA
draw_circos_tracks_points2(
  RVG5_venom_circos_df,
  max_value = max_mRNA_expression,
  track_height = trackheight,
  chromosome_col = 'miRNA.Target.Chrom',
  loci_location_col = 'miRNA.Target.Start',
  value_col = 'RNA.VST.RVG.5S.CV1087.viridis.North.F',
  color_ramp = mRNA_color_ramp,
  chrom_labels = T
  )

# Draw circos tracks for miRNA
draw_circos_tracks_points2(
  RVG5_venom_circos_df,
  max_value = max_miRNA_count,
  track_height = trackheight,
  chromosome_col = 'miRNA.Sequence.Chrom',
  loci_location_col = 'miRNA.Start',
  value_col = 'miRNA.RPM.RVG.5S.CV1087.viridis.North.F',
  color_ramp = miRNA_color_ramp,
  chrom_labels = F
)

# Call the function to draw lines
draw_target_arrows(
  RVG5_venom_circos_df,
  sector1_col = 'miRNA.Sequence.Chrom',
  position1_col = 'miRNA.Start',
  sector2_col = 'miRNA.Target.Chrom',
  position2_col = 'miRNA.Target.Position',
  line_color_col = 'Color'
  )

# Finish pdf
dev.off()


#### CIRCOS Plot for RVG_6S colored by venom gene family ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)
source('Scripts/R/Protein-miRNA_RNA_Analysis/CIRCOS_Functions.R')


# Create smaller data frame for circos plot circos.link function to draw lines
# This data frame also only has RVG_6S information in it and is color coded for venom gene families
RVG6_venom_circos_df <- venom_circos_df %>% 
  dplyr::select(
    -contains('.RVG.5S.'),
    -contains('.RVG.7S.'),
    -contains('.LVG.2.'),
    -contains('.LVG.4.'),
    -contains('.LVG.9.')
  ) %>% 
  filter(!miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M == 0) %>% # filter out any miRNA with zero expression, I don't think there are any, but just to be safe
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  distinct()

# Initialize max count values
max_protein_value <- range(log(RVG6_venom_circos_df$Intensity.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)
max_miRNA_count <- range(log(RVG6_venom_circos_df$miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)
max_mRNA_expression <- range(log(RVG6_venom_circos_df$RNA.VST.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)

# Define breaks for color scales
protein_breaks <- quantile(log(RVG6_venom_circos_df$Intensity.RVG.6S.CV0987.lutosus.Other.M + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
miRNA_breaks <- quantile(log(RVG6_venom_circos_df$miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
mRNA_breaks <- quantile(log(RVG6_venom_circos_df$RNA.VST.RVG.6S.CV0987.lutosus.Other.M + 1), probs = c(0, 0.5, 1), na.rm = TRUE)

# Create color ramps using magma color map
protein_color_ramp <- colorRamp2(protein_breaks, magma(3))
miRNA_color_ramp <- colorRamp2(miRNA_breaks, magma(3))
mRNA_color_ramp <- colorRamp2(mRNA_breaks, magma(3))


# Load scaffold sizes
# This will be used to initiate the sizes of each sector of the
# circos plot
scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
  dplyr::select(-V4, -V5)

# rename columns in scaffold_size
names(scaffold_size) = c('Chrom','size','genome_position')

# create a column of zeros, to indicate the starting
# value of each sector
scaffold_size = scaffold_size %>% mutate(start = 0)

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/3UTR/Venom_Genes_Colored_CIRCOS_CV0987_lutosus_2024.06.24.pdf', width = 10, height = 10)

circos.clear()

# Set basic parameters for the circos plot
circos.par('track.height' = trackheight,
           cell.padding = c(0.01, 0, 0.01, 0),
           start.degree = 90,
           gap.degree = c(rep(1, 17), 20))

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))

# # Draw circos tracks for protein
# draw_circos_tracks_points2(
#   RVG6_venom_circos_df,
#   max_value = max_protein_value,
#   track_height = trackheight,
#   chromosome_col = 'miRNA.Target.Chrom',
#   loci_location_col =  'miRNA.Target.Start',
#   value_col = 'Intensity.RVG.6S.CV0987.lutosus.Other.M',
#   color_ramp = protein_color_ramp,
#   chrom_labels = T
#   )

# Draw circos tracks for mRNA
draw_circos_tracks_points2(
  RVG6_venom_circos_df,
  max_value = max_mRNA_expression,
  track_height = trackheight,
  chromosome_col = 'miRNA.Target.Chrom',
  loci_location_col = 'miRNA.Target.Start',
  value_col = 'RNA.VST.RVG.6S.CV0987.lutosus.Other.M',
  color_ramp = mRNA_color_ramp,
  chrom_labels = T
)

# Draw circos tracks for miRNA
draw_circos_tracks_points2(
  RVG6_venom_circos_df,
  max_value = max_miRNA_count,
  track_height = trackheight,
  chromosome_col = 'miRNA.Sequence.Chrom',
  loci_location_col = 'miRNA.Start',
  value_col = 'miRNA.RPM.RVG.6S.CV0987.lutosus.Other.M',
  color_ramp = miRNA_color_ramp,
  chrom_labels = F
  )


# Call the function to draw lines
draw_target_arrows(RVG6_venom_circos_df,
                   sector1_col = 'miRNA.Sequence.Chrom',
                   position1_col = 'miRNA.Start',
                   sector2_col = 'miRNA.Target.Chrom',
                   position2_col = 'miRNA.Target.Position',
                   line_color_col = 'Color')

# Finish pdf
dev.off()




#### Mini CICROS plot for PLA2 genes only ####


# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
source('Scripts/R/Protein-miRNA_RNA_Analysis/CIRCOS_Functions.R')

# Create smaller data frame for circos plot circos.link function to draw lines
PLA2_circos_df <- mi_df %>% 
  dplyr::select(miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, Positions,
         miRNA.Counts.RVG.5S.CV1087.viridis.North.F, miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M, miRNA.Counts.RVG.7S.CV0985.concolor.Other.F,
         RNA.VST.RVG.5S.CV1087.viridis.North.F, RNA.VST.RVG.6S.CV0987.lutosus.Other.M, RNA.VST.RVG.7S.CV0985.concolor.Other.F,
         Intensity.RVG.5S.CV1087.viridis.North.F, Intensity.RVG.6S.CV0987.lutosus.Other.M, Intensity.RVG.7S.CV0985.concolor.Other.F) %>%
  filter(str_detect(Genes, 'PLA2')) %>% # Filter out everything but PLA2
  #arrange(desc(Intensity.RVG.5S.CV1087.viridis.North.F)) %>%
  filter(!miRNA.Target.Chrom %in% c('PE-reconstructed-10x-myo', 'scaffold-un187')) %>%
  filter(!miRNA.Sequence.Chrom %in% c('scaffold-un11', 'scaffold-un619', 'scaffold-un31', 'scaffold-un147')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
  mutate(Color = case_when(
    str_starts(Genes, 'Venom_PLA2') ~ 'purple',
    str_starts(Genes, 'Venom_SVSP') ~ 'blue',
    str_starts(Genes, 'Venom_SVMP') ~ 'red',
    str_starts(Genes, 'Venom_vQC') ~ 'violet',
    str_starts(Genes, 'Venom_ohanin') ~ 'green',
    str_starts(Genes, 'Venom_LAAO') ~ 'orange',
    str_starts(Genes, 'Venom_CTL') ~ 'magenta',
    str_starts(Genes, 'Venom_CRISP') ~ 'coral',
    str_starts(Genes, 'Venom_BPP') ~ 'cyan',
    str_starts(Genes, 'Venom_VEGF') ~ 'navy',
    str_starts(Genes, 'Venom_myotoxin') ~ 'maroon',
    str_starts(Genes, 'Venom_BPP') ~ 'gray',
    str_starts(Genes, 'Venom_EXO') ~ 'pink',
    str_starts(Genes, 'Venom_ADAM') ~ 'aquamarine',
    T ~ 'black'
  ))


# Create individualized circos data frames for each sample
RVG5_PLA2_circos_df <- PLA2_circos_df %>% 
  dplyr::select(-miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M, -miRNA.Counts.RVG.7S.CV0985.concolor.Other.F,
         -RNA.VST.RVG.6S.CV0987.lutosus.Other.M, -RNA.VST.RVG.7S.CV0985.concolor.Other.F,
         -Intensity.RVG.6S.CV0987.lutosus.Other.M, -Intensity.RVG.7S.CV0985.concolor.Other.F
  ) %>% 
  dplyr::rename(
    miRNA = miRNA.Counts.RVG.5S.CV1087.viridis.North.F,
    Protein = Intensity.RVG.5S.CV1087.viridis.North.F,
    RNA = RNA.VST.RVG.5S.CV1087.viridis.North.F
  ) %>%
  mutate(miRNA.Target.End = miRNA.Target.End - miRNA.Target.Start) %>%
  mutate(miRNA.Target.Start = 0) %>% 
  mutate(miRNA.End = miRNA.End - miRNA.Start) %>% 
  mutate(miRNA.Start = 0) %>% 
  #dplyr::select(-miRNA.Cluster, -miRNA.Sequence.Chrom, -miRNA) %>% 
  distinct()

# Reset the row names of the above dataframe
rownames(RVG5_PLA2_circos_df) <- NULL

# Remove and alter only the miRNA loci from the dataframe
miRNA_scaffold_df <- RVG5_PLA2_circos_df %>% 
  dplyr::select(miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End) %>% 
  dplyr::rename(
    Loci = miRNA.Cluster,
    Chrom = miRNA.Sequence.Chrom,
    Start = miRNA.Start,
    End = miRNA.End
  )

# Now take only the target loci information
target_scaffold_df <- RVG5_PLA2_circos_df %>% 
  dplyr::select(Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End) %>% 
  dplyr::rename(
    Loci = Genes,
    Chrom = miRNA.Target.Chrom,
    Start = miRNA.Target.Start,
    End = miRNA.Target.End
  )

# Fuse the bottom of the two above data frames so that I now have a scaffold file I can use to draw sectors
scaffold_pla2_df <- bind_rows(miRNA_scaffold_df, target_scaffold_df) %>% 
  dplyr::rename(Size = End) %>% 
  distinct()

# Reset indexing because it all weird and now:
rownames(scaffold_pla2_df) <- NULL

# Initialize max count values
max_protein_value <- range(log(RVG5_PLA2_circos_df$Protein + 1), na.rm = T)
max_miRNA_count <- range(log(RVG5_PLA2_circos_df$miRNA + 1), na.rm = T)
max_mRNA_expression <- range(log(RVG5_PLA2_circos_df$RNA + 1), na.rm = T)

# Create filtere

# Set sector height:
sector_height <- 0.2

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/3UTR/PLA_mini_CIRCOS.2024.06.23.pdf', width = 10, height = 10)

circos.clear()

# Set basic parameters for the circos plot
circos.par(track.height = sector_height,
           cell.padding = c(0.01, 0, 0.01, 0),
           start.degree = 90,
           gap.degree = 15)
#gap.degree = c(rep(1, 15), 20))

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_pla2_df$Loci, xlim = cbind(scaffold_pla2_df$Start, scaffold_pla2_df$Size))

# Create data frame with only miRNA information
mi_expression_df <- RVG5_PLA2_circos_df %>% 
  dplyr::select(miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, miRNA) %>% 
  distinct() %>% 
  dplyr::rename(
    Loci = miRNA.Cluster,
    Chrom = miRNA.Sequence.Chrom,
    Start = miRNA.Start,
    End = miRNA.End,
    Expression = miRNA
  ) %>% 
  mutate(
    Normalized.Expression = scale(Expression)
  )

# Now for a data frame for only geneic information
gene_expression_df <- RVG5_PLA2_circos_df %>% 
  dplyr::select(Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, RNA) %>% 
  distinct() %>% 
  dplyr::rename(
    Loci = Genes,
    Chrom = miRNA.Target.Chrom,
    Start = miRNA.Target.Start,
    End = miRNA.Target.End,
    Expression = RNA
  ) %>% 
  mutate(
    Normalized.Expression = scale(Expression)
  )


# Define new circos data frame that combines all of the expression into a single column
combined_expression_df <- bind_rows(mi_expression_df, gene_expression_df) %>% 
  distinct()

circos.track(
  ylim = c(0, sector_height),
  track.height = sector_height,
  panel.fun = function(x, y) {
    
    # Filter data for current chromosome
    # Create variable that holds the sector data using the genes as sectors
    sector <- combined_expression_df$Loci
    
    # Create variable to hold the information for chromomes
    chrom <- combined_expression_df$Chrom
    
    # Variable for the current sector
    current_sector <- combined_expression_df %>% 
      filter(sector == CELL_META$sector.index)
    
    # Iterate over the rows of current_sector
    for (loci in seq_len(nrow(current_sector))) {
      
      # Variables for x and y
      start <- current_sector$Start[loci]
      end <- current_sector$End[loci]
      value <- log(current_sector$Expression[loci] + 1)
      norm_value <- current_sector$Normalized.Expression[loci]
      
      # # Assign color based on value
      # color <- heat.colors(10)[findInterval(value, seq(min(max_mRNA_expression), max(max_mRNA_expression), length.out = 10))]
      # Assign color based on normalized value
      color <- heat.colors(10)[findInterval(norm_value, seq(-3, 3, length.out = 10))]
      
      # Debugging
      print(start)
      print(end)
      print(value)
      print(norm_value)
      
      # Draw rectangles
      circos.rect(xleft = start,
                  xright = end,
                  ybottom = 0,
                  ytop = sector_height,
                  col = color)
      
      # Add text to the CIRCOS that represents the chromosome and the gene
      circos.text(
        CELL_META$xcenter,
        CELL_META$cell.ylim[2] + mm_y(7),
        gsub('Venom_', '', CELL_META$sector.index),
        facing = 'downward',
        cex = 0.6
      )
    }
  }
)

# Again create a new data frame for the data so that I can draw arrows
arrows_PLA2_df <- RVG5_PLA2_circos_df %>% 
  mutate(miRNA.Midpoint = (miRNA.Start + miRNA.End) / 2) %>% 
  mutate(Positions = as.numeric(Positions))

# Plot each miRNA target line one by one with a for loop
for(cluster in 1:nrow(arrows_PLA2_df)) {
  
  # Need to specify sector1, position1 and sector2, position2 to plot the miRNA targeting line
  sector1 <- arrows_PLA2_df$miRNA.Cluster[cluster]
  position1 <- arrows_PLA2_df$miRNA.Midpoint[cluster]
  sector2 <- arrows_PLA2_df$Genes[cluster]
  position2 <- arrows_PLA2_df$Positions[cluster]
  line_and_arrow_color <- arrows_PLA2_df$Color[cluster]
  
  # Print for debugging
  print(arrows_PLA2_df[cluster,])
  
  # Draw links between sectors based on sector and position.
  circos.link(
    sector1,
    position1,
    sector2,
    position2,
    col = line_and_arrow_color,
    #col = add_transparency(line_and_arrow_color, 0.95), # Adjust the line color as needed
    directional = 1,
    lwd = 1,
    arr.width = 0.1,
    arr.length = 0.1
  )
}

# Save
dev.off()

#### CIRCOS Plot for RVG_7S colored by venom gene family ####

# Load CIRCOS plot packages
library(circlize)
library(rlist)
library(officer)
library(rvg)
library(foreach)
library(viridis)
source('Scripts/R/Protein-miRNA_RNA_Analysis/CIRCOS_Functions.R')


# Create smaller data frame for circos plot circos.link function to draw lines
# This data frame also only has RVG_7S information in it and is color coded for venom gene families
RVG7_venom_circos_df <- venom_circos_df %>% 
  dplyr::select(
    -contains('.RVG.5S.'),
    -contains('.RVG.6S.'),
    -contains('.LVG.2.'),
    -contains('.LVG.4.'),
    -contains('.LVG.9.')
  ) %>% 
  filter(!miRNA.RPM.RVG.7S.CV0985.concolor.Other.F == 0) %>% # filter out any miRNA with zero expression, I don't think there are any, but just to be safe
  filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
  distinct()

# Initialize max count values
max_protein_value <- range(log(RVG7_venom_circos_df$Intensity.RVG.7S.CV0985.concolor.Other.F + 1), na.rm = T)
max_miRNA_count <- range(log(RVG7_venom_circos_df$miRNA.RPM.RVG.7S.CV0985.concolor.Other.F + 1), na.rm = T)
max_mRNA_expression <- range(log(RVG7_venom_circos_df$RNA.VST.RVG.7S.CV0985.concolor.Other.F + 1), na.rm = T)

# Define breaks for color scales
protein_breaks <- quantile(log(RVG7_venom_circos_df$Intensity.RVG.7S.CV0985.concolor.Other.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
miRNA_breaks <- quantile(log(RVG7_venom_circos_df$miRNA.RPM.RVG.7S.CV0985.concolor.Other.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)
mRNA_breaks <- quantile(log(RVG7_venom_circos_df$RNA.VST.RVG.7S.CV0985.concolor.Other.F + 1), probs = c(0, 0.5, 1), na.rm = TRUE)

# Create color ramps using magma color map
protein_color_ramp <- colorRamp2(protein_breaks, magma(3))
miRNA_color_ramp <- colorRamp2(miRNA_breaks, magma(3))
mRNA_color_ramp <- colorRamp2(mRNA_breaks, magma(3))

# Load scaffold sizes
# This will be used to initiate the sizes of each sector of the
# circos plot
scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
  dplyr::select(-V4, -V5)

# rename columns in scaffold_size
names(scaffold_size) = c('Chrom','size','genome_position')

# create a column of zeros, to indicate the starting
# value of each sector
scaffold_size = scaffold_size %>% mutate(start = 0)

# Begin saving the plot as a pdf
pdf(file = 'Figures/CIRCOS/3UTR/Venom_Genes_Colored_CIRCOS_CV0985_concolor_2024.06.24.pdf', width = 10, height = 10)

circos.clear()

# Set basic parameters for the circos plot
circos.par('track.height' = trackheight,
           cell.padding = c(0.01, 0, 0.01, 0),
           start.degree = 90,
           gap.degree = c(rep(1, 17), 20))

# Initiate the circos plot
# xlim is used to specify the start and the end value of each chromosome (sectors)
circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))

# # Draw circos tracks for protein
# draw_circos_tracks_points2(
#   RVG7_venom_circos_df,
#   max_value = max_protein_value,
#   track_height = trackheight,
#   chromosome_col = 'miRNA.Target.Chrom',
#   loci_location_col =  'miRNA.Target.Start',
#   value_col = 'Intensity.RVG.7S.CV0985.concolor.Other.F',
#   color_ramp = protein_color_ramp,
#   chrom_labels = T
#   )

# Draw circos tracks for mRNA
draw_circos_tracks_points2(
  RVG7_venom_circos_df,
  max_value = max_mRNA_expression,
  track_height = trackheight,
  chromosome_col = 'miRNA.Target.Chrom',
  loci_location_col = 'miRNA.Target.Start',
  value_col = 'RNA.VST.RVG.7S.CV0985.concolor.Other.F',
  color_ramp = mRNA_color_ramp,
  chrom_labels = T
)

# Draw circos tracks for miRNA
draw_circos_tracks_points2(
  RVG7_venom_circos_df,
  max_value = max_miRNA_count,
  track_height = trackheight,
  chromosome_col = 'miRNA.Sequence.Chrom',
  loci_location_col = 'miRNA.Start',
  value_col = 'miRNA.RPM.RVG.7S.CV0985.concolor.Other.F',
  color_ramp = miRNA_color_ramp,
  chrom_labels = F
  )


# Call the function to draw lines
draw_target_arrows(RVG7_venom_circos_df,
                   sector1_col = 'miRNA.Sequence.Chrom',
                   position1_col = 'miRNA.Start',
                   sector2_col = 'miRNA.Target.Chrom',
                   position2_col = 'miRNA.Target.Position',
                   line_color_col = 'Color')

# Finish pdf
dev.off()



# #### PLA2 Array plot ####
# 
# # Load packages
# library(tidyverse)
# library(rtracklayer)
# library(ggcoverage)
# library(scales)
# library(gggenes)
# library(viridis)
# library(ggforce)
# 
# 
# # CHANGE PATH UP TO /Dropbox
# setwd('~/Dropbox/CastoeLabFolder/projects/z_Done_or_Abandoned_Projects-MSs/_VenomGeneRegulation_NEW_Aug2021')
# 
# # Read priority venom genes data
# pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)
# 
# # Read gene expression data and filter for priority venom genes
# exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
#   filter(txid %in% pri_venom_genes$X6) %>%
#   left_join(pri_venom_genes,by=c('txid'='X6')) %>%
#   mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))
# 
# # Read all gene information and filter for priority venom genes
# all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>%
#   filter(str_detect(X9,'trnascan',negate = T)) %>%
#   mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>%
#   mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>%
#   filter(tx_id %in% pri_venom_genes$X6) %>%
#   left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>%
#   dplyr::select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>%
#   mutate(strand = ifelse(strand == '+','forward','reverse')) %>%
#   mutate(direction = ifelse(strand == 'forward',1,-1)) %>%
#   mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>%
#   left_join(exp) %>%
#   mutate(gene = ifelse(str_detect(gene,'ADAM28'),paste('NVP: ',gene,sep = ''),gene)) %>%
#   mutate(prom_start = ifelse(strand=='forward',start,end))
# 
# # Filter for PLA2 genes
# PLA2_info <- all_info %>%
#   filter(str_detect(gene,'PLA2'))
# 
# 
# # ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.
# 
# # Define region start and end points for PLA2 plot
# PLA2_reg_start = floor(3019890 / 1000 ) * 1000
# PLA2_reg_end = ceiling(3043778 / 1000 ) * 1000
# PLA2_reg_length <- paste(c(round((PLA2_reg_end-PLA2_reg_start)/1000,digits = 2),'kb'),collapse = ' ')
# 
# 
# # Get venom gene expression 
# VST_count_mat <- read.table('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Data/mRNA/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)
# 
# # Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
# VST_count_mat <- VST_count_mat %>%
#   dplyr::select(matches('LVG|RVG')) %>%
#   filter(grepl('Venom', row.names(.)))
# 
# 
# column_numbers <- 1:13
# target_names <- paste0('VG_', column_numbers)
# 
# for (i in column_numbers) {
#   target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))
#   
#   if (length(target_columns) > 0) {
#     if (length(target_columns) == 1) {
#       # If only one column is found, copy it to become VG_i
#       VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
#     } else {
#       # Calculate row mean for multiple columns
#       VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
#     }
#   }
# }
# 
# VST_count_mat_filtered <- VST_count_mat %>%
#   dplyr::select(-contains('LVG'),-contains('RVG')) %>%
#   dplyr::select(contains('VG')) %>%
#   #dplyr::select('RVG_5', 'RVG_6', 'RVG_7', 'RVG_12') %>%
#   rownames_to_column(var = 'Gene') %>%
#   filter(grepl('PLA2', Gene)) %>%
#   mutate(AverageExp = rowMeans(dplyr::select(., starts_with("VG_")), na.rm = T)) %>%
#   dplyr::select(Gene, AverageExp)
# 
# VST_count_mat_filtered <- VST_count_mat_filtered %>%
#   mutate(Gene = gsub('Venom_PLA2', 'PLA2 ', Gene)) %>%
#   mutate(Gene = gsub('Venom_', 'NVP: ', Gene))
# 
# # ADD TO PLA2 INFO
# PLA2_info <- PLA2_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
# PLA2_info <- PLA2_info %>% mutate(gene = gsub(' ','',gene))
# all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
# all_info <- all_info %>% mutate(gene = gsub(' ','',gene))
# 
# # Add venom gene expression data to PLA2_info
# PLA2_info <- PLA2_info %>%
#   mutate(
#     # Swap start and end if direction is -1
#     start2 = case_when(direction == -1 ~ end,
#                        T ~ start),
#     end2 = case_when(direction == -1 ~ start,
#                      T ~ end)
#   )
# 
# # Set color for PLA2_info
# PLA2_info$color = "blue"
# 
# # Plotting 
# # Create plot using ggplot
# PLA2_figure <- ggplot(PLA2_info, aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp+1))) +
#   # Add text labels for genes
#   ggrepel::geom_text_repel(data = all_info %>% 
#                              filter(str_detect(gene,'PLA2')) %>% 
#                              mutate(start = (start + end)/2), 
#                            aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2,size=2.5) +
#   # Add gray segment indicating the region limits                         
#   geom_segment(aes(x=PLA2_reg_start,xend=PLA2_reg_end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
#   # Add gene arrows
#   geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
#   # Set labels and scales
#   ylab('') +
#   xlab('') +
#   scale_fill_viridis_c(option = 'B') +
#   scale_x_continuous(labels = comma,limits=c(PLA2_reg_start,PLA2_reg_end),expand=c(0,0)) +
#   # Set theme
#   theme_classic(base_size = 14) +
#   theme(axis.line.y = element_blank(),
#         plot.title.position = 'plot',
#         plot.title = element_text(color='black',face='bold',size = 14),
#         axis.title.x=element_blank())
# PLA2_figure
# # Save plot to a PDF file
# ggsave("/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Figures/Gene_Arrays/PLA2_array_plot.2024.5.22.pdf", PLA2_figure, width = 10, height = 2, units = "in", create.dir = T)
# 
# 
# 
# #### SVSP Array plot ####
# 
# # Load packages
# library(tidyverse)
# library(rtracklayer)
# library(ggcoverage)
# library(scales)
# library(gggenes)
# library(viridis)
# library(ggforce)
# 
# 
# # CHANGE PATH UP TO /Dropbox
# setwd('~/Dropbox/CastoeLabFolder/projects/z_Done_or_Abandoned_Projects-MSs/_VenomGeneRegulation_NEW_Aug2021')
# 
# pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)
# exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
#   filter(txid %in% pri_venom_genes$X6) %>%
#   left_join(pri_venom_genes,by=c('txid'='X6')) %>%
#   mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))
# 
# 
# all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>%
#   filter(str_detect(X9,'trnascan',negate = T)) %>%
#   mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>%
#   mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>%
#   filter(tx_id %in% pri_venom_genes$X6) %>%
#   left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>%
#   dplyr::select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>%
#   mutate(strand = ifelse(strand == '+','forward','reverse')) %>%
#   mutate(direction = ifelse(strand == 'forward',1,-1)) %>%
#   mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>%
#   left_join(exp) %>%
#   mutate(gene = ifelse(str_detect(gene,'ADAM28'),paste('NVP: ',gene,sep = ''),gene)) %>%
#   mutate(prom_start = ifelse(strand=='forward',start,end))
# 
# SVSP_info <- all_info %>%
#   filter(str_detect(gene,'SVSP'))
# 
# # ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.
# 
# SVSP_reg_start = floor(8568727 / 1000 ) * 1000
# 
# SVSP_reg_end = ceiling(8981362 / 1000 ) * 1000
# 
# SVSP_reg_length <- paste(c(round((SVSP_reg_end-SVSP_reg_start)/1000,digits = 2),'kb'),collapse = ' ')
# 
# 
# # Get venom gene expression 
# VST_count_mat <- read.table('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Data/mRNA/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)
# 
# # Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
# VST_count_mat <- VST_count_mat %>%
#   dplyr::select(matches('LVG|RVG')) %>%
#   filter(grepl('Venom', row.names(.)))
# 
# 
# column_numbers <- 1:13
# target_names <- paste0('VG_', column_numbers)
# 
# for (i in column_numbers) {
#   target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))
#   
#   if (length(target_columns) > 0) {
#     if (length(target_columns) == 1) {
#       # If only one column is found, copy it to become VG_i
#       VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
#     } else {
#       # Calculate row mean for multiple columns
#       VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
#     }
#   }
# }
# 
# VST_count_mat_filtered <- VST_count_mat %>%
#   dplyr::select(-contains('LVG'),-contains('RVG')) %>%
#   dplyr::select(contains('VG')) %>%
#   #dplyr::select('RVG_5', 'RVG_6', 'RVG_7', 'RVG_12') %>%
#   rownames_to_column(var = 'Gene') %>%
#   filter(grepl('SVSP', Gene)) %>%
#   mutate(AverageExp = rowMeans(dplyr::select(., starts_with("VG_")), na.rm = T)) %>%
#   dplyr::select(Gene, AverageExp)
# 
# VST_count_mat_filtered <- VST_count_mat_filtered %>%
#   mutate(Gene = gsub('Venom_SVSP', 'SVSP ', Gene)) %>%
#   mutate(Gene = gsub('Venom_', 'NVP: ', Gene))
# 
# # ADD TO SVSP INFO
# SVSP_info <- SVSP_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
# SVSP_info <- SVSP_info %>% mutate(gene = gsub(' ','',gene))
# all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
# all_info <- all_info %>% mutate(gene = gsub(' ','',gene))
# SVSP_info <- SVSP_info %>%
#   mutate(
#     # Swap start and end if direction is -1
#     start2 = case_when(direction == -1 ~ end,
#                        T ~ start),
#     end2 = case_when(direction == -1 ~ start,
#                      T ~ end)
#   )
# SVSP_info$color = "blue"
# 
# # Plotting 
# SVSP_figure <- ggplot(SVSP_info,aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp+1))) +
#   ggrepel::geom_text_repel(data = all_info %>% 
#                              filter(str_detect(gene,'SVSP')) %>% # Change for SVSP-SVSP
#                              mutate(start = (start + end)/2), 
#                            aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2,size=2.5) +
#   geom_segment(aes(x=SVSP_reg_start,xend=SVSP_reg_end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
#   geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
#   ylab('') +
#   xlab('') +
#   scale_fill_viridis_c(option = 'B') +
#   scale_x_continuous(labels = comma,limits=c(SVSP_reg_start,SVSP_reg_end),expand=c(0,0)) +
#   theme_classic(base_size = 14) +
#   theme(axis.line.y = element_blank(),
#         plot.title.position = 'plot',
#         plot.title = element_text(color='black',face='bold',size = 14),
#         axis.title.x=element_blank())
# SVSP_figure
# # Save SVSP_figure to a .pdf
# ggsave("/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Figures/Gene_Arrays/SVSP_array_plot.2024.5.22.pdf", SVSP_figure, width = 10, height = 2, units = "in", create.dir = T)
# 
# 
# 
# 
# #### SVMP Array plot ####
# 
# # Load packages
# library(tidyverse)
# library(rtracklayer)
# library(ggcoverage)
# library(scales)
# library(gggenes)
# library(viridis)
# library(ggforce)
# 
# # CHANGE PATH UP TO /Dropbox
# setwd('~/Dropbox/CastoeLabFolder/projects/z_Done_or_Abandoned_Projects-MSs/_VenomGeneRegulation_NEW_Aug2021')
# 
# pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)
# exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
#   filter(txid %in% pri_venom_genes$X6) %>%
#   left_join(pri_venom_genes,by=c('txid'='X6')) %>%
#   mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))
# 
# 
# all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>%
#   filter(str_detect(X9,'trnascan',negate = T)) %>%
#   mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>%
#   mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>%
#   filter(tx_id %in% pri_venom_genes$X6) %>%
#   left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>%
#   dplyr::select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>%
#   mutate(strand = ifelse(strand == '+','forward','reverse')) %>%
#   mutate(direction = ifelse(strand == 'forward',1,-1)) %>%
#   mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>%
#   left_join(exp) %>%
#   mutate(gene = ifelse(str_detect(gene,'ADAM28'),paste('NVP: ',gene,sep = ''),gene)) %>%
#   mutate(prom_start = ifelse(strand=='forward',start,end))
# 
# SVMP.info <- all_info %>%
#   filter(str_detect(gene,'SVMP|ADAM'))
# 
# # ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.
# 
# SVMP_reg_start = floor(13901005 / 1000 ) * 1000
# 
# SVMP_reg_end = ceiling(14424729 / 1000 ) * 1000
# 
# SVSP_reg_length <- paste(c(round((SVMP_reg_end-SVMP_reg_start)/1000,digits = 2),'kb'),collapse = ' ')
# 
# 
# # Get venom gene expression
# VST_count_mat <- read.table('Data/mRNA/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)
# 
# # Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
# VST_count_mat <- VST_count_mat %>%
#   dplyr::select(matches('LVG|RVG')) %>%
#   filter(grepl('Venom', row.names(.)))
# 
# 
# column_numbers <- 1:13
# target_names <- paste0('VG_', column_numbers)
# 
# for (i in column_numbers) {
#   target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))
#   
#   if (length(target_columns) > 0) {
#     if (length(target_columns) == 1) {
#       # If only one column is found, copy it to become VG_i
#       VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
#     } else {
#       # Calculate row mean for multiple columns
#       VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
#     }
#   }
# }
# 
# VST_count_mat_filtered <- VST_count_mat %>%
#   dplyr::select(-contains('LVG'),-contains('RVG')) %>%
#   dplyr::select(contains('VG')) %>%
#   #dplyr::select('RVG_5', 'RVG_6', 'RVG_7', 'RVG_12') %>%
#   rownames_to_column(var = 'Gene') %>%
#   filter(grepl('SVMP|ADAM', Gene)) %>%
#   mutate(AverageExp = rowMeans(dplyr::select(., starts_with("VG_")), na.rm = T)) %>%
#   dplyr::select(Gene, AverageExp)
# 
# VST_count_mat_filtered <- VST_count_mat_filtered %>%
#   mutate(Gene = gsub('Venom_SVMP', 'SVMP ', Gene)) %>%
#   mutate(Gene = gsub('Venom_', 'NVP: ', Gene))
# 
# # ADD TO SVMP INFO
# SVMP.info <- SVMP.info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
# SVMP.info <- SVMP.info %>% mutate(gene = gsub(' ','',gene))
# all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
# all_info <- all_info %>% mutate(gene = gsub(' ','',gene))
# SVMP.info <- SVMP.info %>%
#   mutate(
#     # Swap start and end if direction is -1
#     start2 = case_when(direction == -1 ~ end,
#                        T ~ start),
#     end2 = case_when(direction == -1 ~ start,
#                      T ~ end)
#   )
# SVMP.info$color = "blue"
# 
# # Plotting 
# SVMP_figure <- ggplot(SVMP.info,aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp+1))) +
#   ggrepel::geom_text_repel(data = all_info %>% 
#                              filter(str_detect(gene,'SVMP')) %>% # Change for SVSP-SVMP
#                              mutate(start = (start + end)/2), 
#                            aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2,size=2.5) +
#   geom_segment(aes(x=SVMP_reg_start,xend=SVMP_reg_end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
#   geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
#   ylab('') +
#   xlab('') +
#   scale_fill_viridis_c(option = 'B') +
#   scale_x_continuous(labels = comma,limits=c(SVMP_reg_start,SVMP_reg_end),expand=c(0,0)) +
#   theme_classic(base_size = 14) +
#   theme(axis.line.y = element_blank(),
#         plot.title.position = 'plot',
#         plot.title = element_text(color='black',face='bold',size = 14),
#         axis.title.x=element_blank())
# SVMP_figure
# # Save SVMP_figure to a .pdf
# ggsave("/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA/Figures/Gene_Arrays/SVMP_array_plot.2024.5.22.pdf", SVMP_figure, width = 10, height = 2, units = "in", create.dir = T)



#### DESeq2 Analysis ####

# Load in DESeq2
library(DESeq2)
# Load regionReport
library(regionReport)
library("IHW")
library('ashr')
library('vsn')
library('pheatmap')


# Set up

# Set the working directory
setwd('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Documents/Lab/Projects/Venom_grant/ncRNA')
# setwd("C:/Users/kaasb/OneDrive - UT Arlington (1)/Documents/Lab/Projects/Venom_grant/ncRNA/")


# # I am going to exclude the following genes because they are not well annotated
# excluded_genes = c(
#   'maker-scaffold-mi1-augustus-gene-59.13_crovir-transcript-12940',
#   'maker-scaffold-mi1-augustus-gene-59.20_crovir-transcript-12947',
#   'maker-scaffold-mi2-augustus-gene-22.17_crovir-transcript-739',
#   'maker-scaffold-un11-augustus-gene-5.19',
#   'XP_016876419',
#   'XP_011528471'
# )

# Establish path for miRNA-mRNA data
rna_mirna_data <- 'Data/Merged/miRNA_mRNA_Combined_Data_IMPORTANT.2024.06.14.tsv'
# Establish the path for the sample metadata
metadata <- 'Data/metadata/DESeq_Sample_metadata_2024.6.4.csv'

# Read in other miRNA-mRNA data
rna_mirna_df <- read.table(file = rna_mirna_data, header = T)  %>% 
  dplyr::rename(
    'miRNA.Cluster.Original' = 'miRNA.Cluster',
    'Genes' = 'Converted.Gene.IDs',
    'miRNA.Cluster' = 'Putative.miRNA.Name'
  ) %>%
  dplyr::select(miRNA.Cluster, everything()) %>% # Move the new miRNA.Clusters to the front
  # filter(!(Genes %in% excluded_genes)) %>% 
  # filter(str_starts(Genes, 'Venom_|PLA2G2E.1')) %>%
  filter(Origin == 'three_prime_utr') %>% # Filter out everything not targeting the 3' UTR
  dplyr::select(-Origin)  # Remove to save memory


# Create miRNA data frame with only miRNAs
mirna_df <- rna_mirna_df %>% 
  # filter(str_starts(Genes, 'Venom')) %>% 
  dplyr::select(miRNA.Cluster, contains('miRNA.Counts')) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(
    Genes = miRNA.Cluster,
    CV1081.viridis = miRNA.Counts.LVG.2.CV1081.viridis.Mid.M, CV0857.viridis = miRNA.Counts.LVG.4.CV0857.viridis.North.M, CV1086.viridis = miRNA.Counts.LVG.9.CV1086.viridis.South.M, CV1082.viridis = miRNA.Counts.RVG.12S.CV1082.viridis.South.M,
    CV1087.viridis = miRNA.Counts.RVG.5S.CV1087.viridis.North.F, CV0987.lutosus = miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M, CV0985.concolor = miRNA.Counts.RVG.7S.CV0985.concolor.Other.F
  ) 
rownames(mirna_df) <- NULL # Reset row names


# Create mRNA data frame with only mRNA
mrna_df <- rna_mirna_df %>% 
  # filter(str_starts(Genes, 'Venom')) %>% 
  dplyr::select(Genes, contains('RNA.Raw.Counts.')) %>%
  dplyr::distinct() %>% 
  dplyr::rename(
    CV1081.viridis = RNA.Raw.Counts.LVG.2.CV1081.viridis.Mid.M, CV0857.viridis = RNA.Raw.Counts.LVG.4.CV0857.viridis.North.M, CV1086.viridis = RNA.Raw.Counts.LVG.9.CV1086.viridis.South.M, CV1082.viridis = RNA.Raw.Counts.RVG.12S.CV1082.viridis.South.M,
    CV1087.viridis = RNA.Raw.Counts.RVG.5S.CV1087.viridis.North.F, CV0987.lutosus = RNA.Raw.Counts.RVG.6S.CV0987.lutosus.Other.M, CV0985.concolor = RNA.Raw.Counts.RVG.7S.CV0985.concolor.Other.F
  )
rownames(mrna_df) <- NULL

# Concatenate the data frames
counts_matrix <- bind_rows(mrna_df, mirna_df)

# Turn genes column into rownames
rownames(counts_matrix) <- counts_matrix$Genes
# Remove the genes column
counts_matrix <- counts_matrix %>% dplyr::select(-Genes)

# Read the sample data
sample_data <- read.csv(metadata, row.names = 1)
sample_data$species <- factor(sample_data$species) # Turn the sample data into factors
sample_data$sub_species <- factor(sample_data$sub_species) # Turn the sample data into factors
sample_data$state <- factor(sample_data$state)
sample_data$locality <- factor(sample_data$locality)
sample_data$sex <- factor(sample_data$sex)
sample_data$tissue <- factor(sample_data$tissue)


# Create deseq data set for miRNAs
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_data,
  design = ~ species
)
dds


# Run DESeq for the miRNA data
dds <- DESeq(dds)
res <- results(dds)

# Trying an MA-plot on both
plotMA(res, ylim = c(-2,2))
# Nothing really stands out

# Let's create reports for both
# Create report to do some initial analysis
report1 <- DESeq2Report(
  dds = dds,
  res = res,
  intgroup = c('species'),
  project = 'miRNA-mRNA Analysis',
  author = 'Kaas Ballard',
  outdir = 'Figures/DESeq2/3UTR/DESeq2_miRNA-mRNA-Region_Report',
  output = 'Species_vs_Species'
)
browseURL(report1)


# Extract results for miRNAs and mRNAs separately
res_combined_df <- as.data.frame(res)
res_combined_df$Genes <- rownames(res_combined_df) # Create new column based off of rownames


# Identify miRNA and mRNA results based on prefixes
# Filter for any miRNAs
res_mirna <- res_combined_df %>% filter(grepl("^cvi-|^Cluster_", Genes)) %>% dplyr::rename(miRNA.Cluster = Genes)
rownames(res_mirna) <- NULL
# Filter out any miRNAs
res_mrna <- res_combined_df %>% filter(!grepl("^cvi-|^Cluster_", Genes))
rownames(res_mrna) <- NULL

# Create a data frame for interaction data
interactions_df <- rna_mirna_df %>% dplyr::select(miRNA.Cluster, Genes) %>% distinct()

# Join results into the same data frame again
merged_diff_df <- interactions_df %>% 
  left_join(res_mirna, by = 'miRNA.Cluster') %>% 
  left_join(res_mrna, by = 'Genes') %>% 
  dplyr::rename(
    baseMean.mi = baseMean.x, log2FoldChange.mi = log2FoldChange.x, lfcSE.mi = lfcSE.x, stat.mi = stat.x, pvalue.mi = pvalue.x, padj.mi = padj.x,
    baseMean.mr = baseMean.y, log2FoldChange.mr = log2FoldChange.y, lfcSE.mr = lfcSE.y, stat.mr = stat.y, pvalue.mr = pvalue.y, padj.mr = padj.y
  )

# Filter down the above to only significant pairs
significant_diff_paired_df <- merged_diff_df %>% 
  filter(padj.mi < 0.05 & padj.mr < 0.05) %>% 
  filter(!is.na(padj.mi) & !is.na(padj.mr))
# Filter for each individually as well
significant_mi_diff_df <- merged_diff_df %>% filter(padj.mi < 0.05) %>% filter(!is.na(padj.mi))
significant_mr_diff_df <- merged_diff_df %>% filter(padj.mr < 0.05) %>% filter(!is.na(padj.mr))

# Extract normalized counts for significant pairs
significant_genes <- unique(c(significant_diff_paired_df$miRNA.Cluster, significant_diff_paired_df$Genes))
normalized_counts <- counts(dds, normalized = TRUE)
sig_normalized_counts <- normalized_counts[significant_genes, ]

# Create a heatmap
pheatmap(sig_normalized_counts, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = TRUE, show_colnames = TRUE,
         main = "Heatmap of Co-differentially Expressed miRNA-mRNA Pairs")
# Doesn't look promising

# # Scatter plot of log2 fold changes
# ggplot(significant_diff_paired_df, aes(x = log2FoldChange.mi, y = log2FoldChange.mr)) +
#   geom_point(alpha = 0.6) +
#   labs(x = "log2 Fold Change of miRNA", y = "log2 Fold Change of mRNA", 
#        title = "Scatter Plot of log2 Fold Changes for Co-differentially Expressed miRNA-mRNA Pairs") +
#   theme_minimal()


# Create PCA for this data
# Transform and extracting transformed values
vsd <- vst(dds, blind = F)

# Plot different interactions
plotPCA(vsd, intgroup = c("sub_species", "sex"))
plotPCA(vsd, intgroup = c('species'))

# Extract data for ggplot PCA
pca_df <- plotPCA(vsd, intgroup = c('species'), returnData = T)

# Create percent variation numbers
percentVar <- round(100 * attr(pca_df, 'percentVar'))

# Plot using ggplot
pca1 <- ggplot(
  data = pca_df, aes(PC1, PC2, color = species)
) +
  geom_point(size = 3) +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) +
  coord_fixed() +
  labs(title = 'PCA for differential expression between species') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  stat_ellipse(level = 0.95)
pca1
ggsave('Figures/DESeq2/3UTR/PCA/Species_DESeq_PCA.2024.06.23.pdf', pca1, create.dir = T)



# Gene family plots

# SVMP
plotCounts(dds, gene = 'Venom_SVMP1', intgroup = 'species')
# Extract data
svmp1_df <- plotCounts(dds, gene = 'Venom_SVMP1', intgroup = 'species', returnData = T)
svmp1_plot <- ggplot(svmp1_df, aes(x = species, y = count)) +
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  scale_y_log10(breaks = c(25, 100, 400)) +
  labs(y = 'Normalized Counts', x = 'Species', title = 'SVMP1 Differential Expression') +
  theme(plot.title = element_text(hjust = 0.5))
svmp1_plot

# Import counts plotting function
source('Scripts/R/Protein-miRNA_RNA_Analysis/Plot_Gene_Counts_Function.R')


#### DESeq Count Plots For Venom Genes and Targeting miRNAs ####

# Create a merged diff for only venom genes
venom_merged_diff_df <- merged_diff_df %>% filter(str_detect(Genes, 'Venom_|PLA2G2E.1'))
# Create a new df for the genes only
venom_genes_only_df <- venom_merged_diff_df %>% dplyr::select(-miRNA.Cluster, -contains('.mi')) %>% distinct()
# Create a new df for the miRNAs only
venom_mirna_only_df <- venom_merged_diff_df %>% dplyr::select(-Genes, -contains('.mr')) %>% distinct()

# Create a list of genes to be run the function on based on the merged data frame
genes <- venom_genes_only_df$Genes

# Run the function in a loop
for (gene in genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Venom/Genes/Significant_and_Nonsignificant', date = '2024.06.23')
}

# Create a list of miRNAs
miRNAs <- venom_mirna_only_df$miRNA.Cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Venom/miRNAs/Significant_and_Nonsignificant', date = '2024.06.23')
}


#### DESeq Count Plots For Significant Venom Genes and Targeting miRNA pairs ####

# Create a data frame of significant venom genes
venom_significant_diff_paired_df <- significant_diff_paired_df %>% filter(str_detect(Genes, 'Venom_|PLA2GE.1')) %>% distinct()
# Create a new df for the genes only
venom_sig_genes_paired_df <- venom_significant_diff_paired_df %>% dplyr::select(-miRNA.Cluster, -contains('.mi')) %>% distinct()
# Create a new df for the miRNAs only
venom_sig_mi_paired_df <- venom_significant_diff_paired_df %>% dplyr::select(-Genes, -contains('.mr')) %>% distinct()


# Create a list of genes that are significantly differentially expressed
sig_genes <- venom_sig_genes_paired_df$Genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Venom/Genes/Significant_Paired', date = '2024.06.23')
}

# Create a list of miRNAs
sig_miRNAs <- venom_sig_mi_paired_df$miRNA.Cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Venom/miRNAs/Significant_Paired', date = '2024.06.23')
}


#### DESeq Count Plots For Significant Venom Genes and Targeting miRNAs that are not paired ####

# Create a data frame of significant venom genes (not paired)
venom_sig_mr_diff_df <- significant_mr_diff_df %>% filter(str_detect(Genes, 'Venom_|PLA2GE.1')) %>% dplyr::select(-miRNA.Cluster, -contains('.mi')) %>% distinct()
# Create a new df for the miRNAs only
venom_sig_mi_diff_df <- significant_mi_diff_df %>% filter(str_detect(Genes, 'Venom_|PLA2GE.1')) %>% dplyr::select(-Genes, -contains('.mr')) %>% distinct()


# Create a list of genes that are significantly differentially expressed
sig_genes <- venom_sig_mr_diff_df$Genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Venom/Genes/Significant', date = '2024.06.23')
}

# Create a list of miRNAs
sig_miRNAs <- venom_sig_mi_diff_df$miRNA.Cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Venom/miRNAs/Significant', date = '2024.06.23')
}



#### DESeq Count Plots For Non-Venom Genes and Targeting miRNAs ####

# Create a merged diff for only nonvenom genes
non_venom_merged_diff_df <- merged_diff_df %>% filter(!str_detect(Genes, 'Venom_|PLA2G2E.1')) %>% distinct()
# # Create a new df for the genes only
# non_venom_genes_df <- non_venom_merged_diff_df %>% dplyr::select(-miRNA.Cluster) %>% distinct()
# Create a new df for the miRNAs only
non_venom_mirna_df <- non_venom_merged_diff_df %>% dplyr::select(-Genes, -contains('.mr')) %>% distinct()

# I don't actually care about this
# # Create a list of genes to be run the function on based on the merged data frame
# genes <- non_venom_genes_df$Genes
# 
# # Run the function in a loop
# for (gene in genes) {
#   plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Non_Venom/Genes/Significant_and_Nonsignificant', date = '2024.06.23')
# }

# Create a list of miRNAs
miRNAs <- non_venom_mirna_df$miRNA.Cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Non_Venom/miRNAs/Significant_and_Nonsignificant', date = '2024.06.23')
}



# Create a data frame of significant non-venom genes
non_venom_significant_diff_paired_df <- significant_diff_paired_df %>% filter(!str_detect(Genes, 'Venom_|PLA2GE.1')) %>% distinct()
# Create a new df for the genes only
non_venom_sig_genes_paired_df <- non_venom_significant_diff_paired_df %>% dplyr::select(-miRNA.Cluster, -contains('.mi')) %>% distinct()
# Create a new df for the miRNAs only
non_venom_sig_mi_paired_df <- non_venom_significant_diff_paired_df %>% dplyr::select(-Genes, -contains('.mr')) %>% distinct()


# Create a list of genes that are significantly differentially expressed
sig_genes <- non_venom_sig_genes_paired_df$Genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Non_Venom/Genes/Significant_Paired', date = '2024.06.23')
}

# Create a list of miRNAs
sig_miRNAs <- non_venom_sig_mi_paired_df$miRNA.Cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Non_Venom/miRNAs/Significant_Paired', date = '2024.06.23')
}


#### DESeq Count Plots For Significant Non-Venom Genes and Targeting miRNAs that are not paired ####

# Create a data frame of significant venom genes (not paired)
non_venom_sig_mr_diff_df <- significant_mr_diff_df %>% filter(!str_detect(Genes, 'Venom_|PLA2GE.1')) %>% dplyr::select(-miRNA.Cluster, -contains('.mi')) %>% distinct()
# Create a new df for the miRNAs only
non_venom_sig_mi_diff_df <- significant_mi_diff_df %>% filter(!str_detect(Genes, 'Venom_|PLA2GE.1')) %>% dplyr::select(-Genes, -contains('.mr')) %>% distinct()


# Create a list of genes that are significantly differentially expressed
sig_genes <- non_venom_sig_mr_diff_df$Genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Non_Venom/Genes/Significant', date = '2024.06.23')
}

# Create a list of miRNAs
sig_miRNAs <- non_venom_sig_mi_diff_df$miRNA.Cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/3UTR/DESeq_Counts/Species/Non_Venom/miRNAs/Significant', date = '2024.06.23')
}




#### Partial Correlation Analysis ####

# Load package
library(ppcor)

# Get normalized vst counts
vsd <- vst(dds, blind = F)

# Extract the vst transformed values
vst <- as.data.frame(assay(vsd)) %>% dplyr::rename(
  LVG.2.CV1081.viridis.Mid.M = CV1081.viridis, 
  LVG.4.CV0857.viridis.North.M = CV0857.viridis, 
  LVG.9.CV1086.viridis.South.M = CV1086.viridis,
  RVG.12S.CV1082.viridis.South.M = CV1082.viridis, 
  RVG.5S.CV1087.viridis.North.F = CV1087.viridis,
  RVG.6S.CV0987.lutosus.Other.M = CV0987.lutosus, 
  RVG.7S.CV0985.concolor.Other.F = CV0985.concolor
) 
# Turn rows into columns
vst <- rownames_to_column(vst, var = "Genes")

# Create a data frame only for the miRNAs
mi_vst <- vst %>% 
  filter(str_detect(Genes, 'cvi-|Cluster_')) %>% 
  pivot_longer(
    cols = contains('CV'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.VST'
  ) %>% 
  dplyr::rename(miRNA.Cluster = Genes)


# Create a data frame only for the regular genes
genes_vst <- vst %>% 
  filter(!str_detect(Genes, 'cvi-|Cluster_')) %>% 
  pivot_longer(
    cols = contains('CV'),
    names_to = 'Sample.ID',
    values_to = 'mRNA.VST'
  )

# Create a data frame for interaction data
interactions_df <- rna_mirna_df %>% dplyr::select(miRNA.Cluster, Genes, Max.Energy, Max.Score) %>% distinct() %>% 
  group_by(Genes) %>% 
  filter(Max.Energy == min(Max.Energy, na.rm = T)) %>% 
  ungroup() %>% 
  dplyr::select(-contains('Max.'))

# Find shared columns between the above and the normalized miRNA count data
shared_columns <- intersect(names(mi_vst), names(interactions_df))

# Fuse miRNA data with it's interactions
mi_vst <- mi_vst %>% left_join(interactions_df, by = shared_columns, relationship = 'many-to-many')

# Find shared columns between mi_vst and genes_vst
shared_columns <- intersect(names(mi_vst), names(genes_vst))

# Fuse data frames
mi_mrna_vst_df <- left_join(mi_vst, genes_vst, by = shared_columns) %>% dplyr::select(Sample.ID, miRNA.Cluster, miRNA.VST, Genes, mRNA.VST)

# Read protein data in:
# Set and read the data
protein_data <- 'Data/Protein/Formated_Protein_Data_IMPORTANT_2024.06.14.tsv'
protein_df <- read.table(file = protein_data, header = T) %>% 
  dplyr::rename('Genes' = 'Converted.Gene.IDs') %>%
  filter(Feature == 'Intensity') %>% 
  dplyr::select(Intensity, Sample.ID, Genes) %>% 
  dplyr::rename(Protein.Intensity = Intensity)

# Find shared columns between mi_vst and genes_vst
shared_columns <- intersect(names(mi_mrna_vst_df), names(protein_df))

# Fuse data frames
partial_correlation_df <- left_join(mi_mrna_vst_df, protein_df, by = shared_columns) %>% dplyr::select(Sample.ID, miRNA.Cluster, miRNA.VST, Genes, mRNA.VST, Protein.Intensity) %>% 
  filter(!is.na(Protein.Intensity))

# Cut down data size and rename stuff

partial_corr_df <- partial_correlation_df %>% 
  dplyr::select(miRNA.VST, mRNA.VST, Protein.Intensity) %>% 
  dplyr::rename(
    Z = miRNA.VST,
    X = mRNA.VST,
    Y = Protein.Intensity
  ) %>% 
  dplyr::mutate(across(everything(), as.numeric))

# Actually do the partial correlation
# Calculate partial correlation
partial_correlation_results <- pcor(partial_corr_df)$estimate
# Doesn't work because each X has multiple Zs



#### Multiple Regression ####

# Load package
library(effects)


# Get normalized vst counts
vsd <- vst(dds, blind = F)

# Extract the vst transformed values
vst <- as.data.frame(assay(vsd)) %>% dplyr::rename(
  LVG.2.CV1081.viridis.Mid.M = CV1081.viridis, 
  LVG.4.CV0857.viridis.North.M = CV0857.viridis, 
  LVG.9.CV1086.viridis.South.M = CV1086.viridis,
  RVG.12S.CV1082.viridis.South.M = CV1082.viridis, 
  RVG.5S.CV1087.viridis.North.F = CV1087.viridis,
  RVG.6S.CV0987.lutosus.Other.M = CV0987.lutosus, 
  RVG.7S.CV0985.concolor.Other.F = CV0985.concolor
) 
# Turn rows into columns
vst <- rownames_to_column(vst, var = "Genes")

# Create a data frame only for the miRNAs
mi_vst <- vst %>% 
  filter(str_detect(Genes, 'cvi-|Cluster_')) %>% 
  pivot_longer(
    cols = contains('CV'),
    names_to = 'Sample.ID',
    values_to = 'miRNA.VST'
  ) %>% 
  dplyr::rename(miRNA.Cluster = Genes)


# Create a data frame only for the regular genes
genes_vst <- vst %>% 
  filter(!str_detect(Genes, 'cvi-|Cluster_')) %>% 
  pivot_longer(
    cols = contains('CV'),
    names_to = 'Sample.ID',
    values_to = 'mRNA.VST'
  )

# Create a data frame for interaction data
interactions_df <- rna_mirna_df %>% dplyr::select(miRNA.Cluster, Genes, Max.Energy, Max.Score) %>% distinct() %>% 
  group_by(Genes) %>%
  filter(Max.Energy == min(Max.Energy, na.rm = T)) %>%
  ungroup() %>%
  dplyr::select(-contains('Max.'))

# Find shared columns between the above and the normalized miRNA count data
shared_columns <- intersect(names(mi_vst), names(interactions_df))

# Fuse miRNA data with it's interactions
mi_vst <- mi_vst %>% left_join(interactions_df, by = shared_columns, relationship = 'many-to-many')

# Find shared columns between mi_vst and genes_vst
shared_columns <- intersect(names(mi_vst), names(genes_vst))

# Fuse data frames
mi_mrna_vst_df <- left_join(mi_vst, genes_vst, by = shared_columns) %>% dplyr::select(Sample.ID, miRNA.Cluster, miRNA.VST, Genes, mRNA.VST)

# Read protein data in:
# Set and read the data
protein_data <- 'Data/Protein/Formated_Protein_Data_IMPORTANT_2024.06.14.tsv'
protein_df <- read.table(file = protein_data, header = T) %>% 
  dplyr::rename('Genes' = 'Converted.Gene.IDs') %>%
  filter(Feature == 'Intensity') %>% 
  dplyr::select(Intensity, Sample.ID, Genes) %>% 
  dplyr::rename(Protein.Intensity = Intensity)

# Find shared columns between mi_vst and genes_vst
shared_columns <- intersect(names(mi_mrna_vst_df), names(protein_df))

# Fuse data frames
multiple_correlation_df <- left_join(mi_mrna_vst_df, protein_df, by = shared_columns) %>% dplyr::select(Sample.ID, miRNA.Cluster, miRNA.VST, Genes, mRNA.VST, Protein.Intensity) %>% 
  filter(!is.na(Protein.Intensity)) %>% 
  filter(str_detect(Genes, 'Venom_|PLA2G2E.1'))

# Fit a multiple linear regression model to the mRNA and miRNA data
mlr_model <- lm(Protein.Intensity ~ mRNA.VST - miRNA.VST, data = multiple_correlation_df)

# Show summary
summary(mlr_model)


# Visualize model for mRNAs effect on protein
mlr_mrna_plot <- ggplot(data = multiple_correlation_df, aes(x = mRNA.VST, y = Protein.Intensity)) +
  geom_point() +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x - 1
  ) +
  labs(title = 'Protein Intensity vs. mRNA VST')
mlr_mrna_plot

# Visualize model for miRNAs effect on protein
mlr_mirna_plot <- ggplot(multiple_correlation_df, aes(x = miRNA.VST, y = Protein.Intensity)) +
  geom_point() +
  geom_smooth(
    method = 'lm',
    se = T,
    color = 'black',
    linetype = 'dashed',
    formula = y ~ x - 1
  ) +
  labs(title = "Protein Intensity vs. miRNA VST")
mlr_mirna_plot

# Create effect plots
effect_plot_miRNA <- allEffects(mlr_model)[["miRNA.VST"]]
effect_plot_mRNA <- allEffects(mlr_model)[["mRNA.VST"]]
# Didn't work either
