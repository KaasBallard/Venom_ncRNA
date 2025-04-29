# Last Edited: 2024/9/14


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
library(scales)

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


# Let's keep the analysis limited to the 3UTR
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
  dplyr::select(-'gtf.gene', -'crovir.transcript', -'Protein.Probability', -'Top.Peptide.Probability', -'Blast.Alignment.Length')  # Remove columns to save memory
rm(miRNA_mRNA_protein_df)

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



#### Binding Data for PLA2, SVSP, SVMP ####

# Create a data frame to hold venom gene info
venom_genes_table <- mi_df %>% 
  filter(
    str_detect(Genes, 'Venom_|PLA2G2E.1|PLA2_A1')
  ) %>% 
  distinct(Genes, miRNA.Cluster, Origin)
rownames(venom_genes_table) <- NULL

# Create a table for miRNAs that target PLA2
pla2_table_3utr <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'three_prime_utr'),
    str_detect(Genes, 'PLA2')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(pla2_table_3utr, 'Figures/Gene_Arrays/PLA2_three_utr_targets_2024.09.14.csv', row.names = F)

pla2_table_5utr <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'five_prime_utr'),
    str_detect(Genes, 'PLA2')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(pla2_table_5utr, 'Figures/Gene_Arrays/PLA2_five_utr_targets_2024.09.14.csv', row.names = F)


pla2_table_cds <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'CDS'),
    str_detect(Genes, 'PLA2')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(pla2_table_cds, 'Figures/Gene_Arrays/PLA2_cds_targets_2024.09.14.csv', row.names = F)



# Create a table for miRNAs that target SVSP
SVSP_table_3utr <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'three_prime_utr'),
    str_detect(Genes, 'SVSP')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(SVSP_table_3utr, 'Figures/Gene_Arrays/SVSP_three_utr_targets_2024.09.14.csv', row.names = F)

SVSP_table_5utr <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'five_prime_utr'),
    str_detect(Genes, 'SVSP')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(SVSP_table_5utr, 'Figures/Gene_Arrays/SVSP_five_utr_targets_2024.09.14.csv', row.names = F)


SVSP_table_cds <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'CDS'),
    str_detect(Genes, 'SVSP')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(SVSP_table_cds, 'Figures/Gene_Arrays/SVSP_cds_targets_2024.09.14.csv', row.names = F)


# Create a table for miRNAs that target SVMP
SVMP_table_3utr <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'three_prime_utr'),
    str_detect(Genes, 'SVMP')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(SVMP_table_3utr, 'Figures/Gene_Arrays/SVMP_three_utr_targets_2024.09.14.csv', row.names = F)

SVMP_table_5utr <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'five_prime_utr'),
    str_detect(Genes, 'SVMP')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(SVMP_table_5utr, 'Figures/Gene_Arrays/SVMP_five_utr_targets_2024.09.14.csv', row.names = F)


SVMP_table_cds <- venom_genes_table %>% 
  filter(
    str_detect(Origin, 'CDS'),
    str_detect(Genes, 'SVMP')
  ) %>% 
  distinct() %>% 
  group_by(Genes) %>% 
  summarize(
    miRNA.Cluster = paste(miRNA.Cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  )
write.csv(SVMP_table_cds, 'Figures/Gene_Arrays/SVMP_cds_targets_2024.09.14.csv', row.names = F)




#### PLA2 Array plot ####

# Load packages
library(tidyverse)
library(rtracklayer)
library(scales)
library(gggenes)
library(viridis)
library(ggforce)
library(cowplot)


# CHANGE PATH UP TO /Dropbox
setwd('~/Dropbox/CastoeLabFolder/projects/z_Done_or_Abandoned_Projects-MSs/_VenomGeneRegulation_NEW_Aug2021')

# Read priority venom genes data
pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt', col_names = F)

# Read gene expression data and filter for priority venom genes
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
  select(-contains('RVG')) %>% 
  filter(txid %in% pri_venom_genes$X6) %>% # Filter to keep only transcript IDs present in the priority venom genes list
  left_join(pri_venom_genes, by = c('txid'='X6')) %>%
  mutate(gene = ifelse(str_detect(X7, 'ADAM28', negate = T),str_replace_all(X7, '_', ' '), X7))

# Read all gene information and filter for priority venom genes
all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff', col_names = F) %>%
  filter(str_detect(X9, 'trnascan', negate = T)) %>%
  mutate(txid = str_split_fixed(X9, ';', 4)[ , 3]) %>%
  mutate(txid = str_remove_all(txid, 'Crovir_Transcript_ID=')) %>%
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes, by = c('txid' = 'X6')) %>%
  dplyr::select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7, txid) %>%
  mutate(strand = ifelse(strand == '+','forward', 'reverse')) %>%
  mutate(direction = ifelse(strand == 'forward', 1, -1)) %>%
  mutate(gene = ifelse(str_detect(gene, 'ADAM28', negate = T), str_replace_all(gene, '_', ' '), gene)) %>%
  left_join(exp) %>%
  mutate(gene = ifelse(str_detect(gene, 'ADAM28'), paste('NVP: ', gene, sep = ''), gene)) %>%
  mutate(prom_start = ifelse(strand == 'forward', start, end))

# Filter for PLA2 genes
PLA2_info <- all_info %>%
  filter(str_detect(gene, 'PLA2'))


# ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.

# Define region start and end points for PLA2 plot
PLA2_reg_start = floor(3019890 / 1000 ) * 1000
PLA2_reg_end = ceiling(3043778 / 1000 ) * 1000
PLA2_reg_length <- paste(c(round((PLA2_reg_end-PLA2_reg_start)/1000, digits = 2), 'kb'), collapse = ' ')


# Get venom gene expression
VST_count_mat <- read.table('/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Data/mRNA/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)

# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>%
  dplyr::select(matches('LVG|RVG')) %>%
  filter(grepl('Venom|PLA2G2E.1', row.names(.)))

# Create a list of column numbers to operate over
column_numbers <- 1:13
# Use the above list to create a list of VG_i columns to operate over
target_names <- paste0('VG_', column_numbers)

# Calculate a pairwise mean for each sample
for (i in column_numbers) {
  target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))

  if (length(target_columns) > 0) {
    if (length(target_columns) == 1) {
      # If only one column is found, copy it to become VG_i
      VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
    } else {
      # Calculate row mean for multiple columns
      VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
    }
  }
}

VST_count_mat_filtered <- VST_count_mat %>%
  dplyr::select(-contains('LVG'), -contains('RVG')) %>%
  dplyr::select(contains('VG')) %>%
  dplyr::select('VG_5', 'VG_6', 'VG_7', 'VG_9', 'VG_2', 'VG_4') %>% # Filter for samples being used in study
  rownames_to_column(var = 'Gene') %>%
  filter(grepl('PLA2', Gene)) %>%
  mutate(AverageExp = rowMeans(dplyr::select(., starts_with("VG_")), na.rm = T)) %>%
  dplyr::select(Gene, AverageExp)

VST_count_mat_filtered <- VST_count_mat_filtered %>%
  mutate(Gene = gsub('Venom_PLA2', 'PLA2', Gene)) %>%
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene)) %>% 
  mutate(Gene = replace(Gene, Gene == 'PLA2G2E.1', 'PLA2gIIE'))

# ADD TO PLA2 INFO
PLA2_info <- PLA2_info %>% mutate(gene = gsub(' ', '', gene))
PLA2_info <- PLA2_info %>% left_join(VST_count_mat_filtered, by = c('gene' = 'Gene')) %>% mutate(gene = gsub(' ', '', gene))
all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene' = 'Gene'))
all_info <- all_info %>% mutate(gene = gsub(' ', '', gene))

# Add venom gene expression data to PLA2_info
PLA2_info <- PLA2_info %>%
  mutate(
    # Swap start and end if direction is -1
    start2 = case_when(direction == -1 ~ end,
                       T ~ start),
    end2 = case_when(direction == -1 ~ start,
                     T ~ end)
  )

# Set color for PLA2_info
PLA2_info$color = "blue"

# Plotting
# Create plot using ggplot
PLA2_figure <- ggplot(PLA2_info, aes(xmin = start2, xmax = end2, y = str_remove(molecule, 'scaffold\\-'), fill = log10(AverageExp + 1))) +
  # Add text labels for genes
  ggrepel::geom_text_repel(data = all_info %>%
                             filter(str_detect(gene, 'PLA2')) %>%
                             mutate(start = (start + end)/2),
                           aes(x = start, y = str_remove(molecule, 'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2, size = 2.5) +
  # Add gray segment indicating the region limits
  geom_segment(aes(x = PLA2_reg_start, xend = PLA2_reg_end, y = str_remove(molecule, 'scaffold\\-'), yend = str_remove(molecule, 'scaffold\\-')), lwd = 1, color = 'grey70') +
  # Add gene arrows
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"), show.legend = F) +
  # Set labels and scales
  ylab('') +
  xlab('') +
  labs(
    title = 'a. PLA2s'
  ) +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma, limits = c(PLA2_reg_start, PLA2_reg_end), expand = c(0, 0)) +
  # Set theme
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    plot.title.position = 'plot',
    plot.title = element_text(color= 'black', hjust = 0.5, face= 'bold', size = 14),
    axis.title.x=element_blank()
  )
PLA2_figure
# Save plot to a PDF file
ggsave("/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Figures/Gene_Arrays/PLA2_array_plot.2024.09.17.pdf", PLA2_figure, width = 10, height = 2, units = "in", create.dir = T)



#### SVSP Array plot ####

# Load packages
library(tidyverse)
library(rtracklayer)
library(scales)
library(gggenes)
library(viridis)
library(ggforce)


# CHANGE PATH UP TO /Dropbox
setwd('~/Dropbox/CastoeLabFolder/projects/z_Done_or_Abandoned_Projects-MSs/_VenomGeneRegulation_NEW_Aug2021')

# Read priority venom genes data
pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)

# Read gene expression data and filter for priority venom genes
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
  select(-contains('RVG')) %>%
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes, by = c('txid' = 'X6')) %>%
  mutate(gene = ifelse(str_detect(X7, 'ADAM28', negate = T), str_replace_all(X7, '_', ' '), X7))

# Read all gene information and filter for priority venom genes
all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff', col_names = F) %>%
  filter(str_detect(X9, 'trnascan', negate = T)) %>%
  mutate(tx_id = str_split_fixed(X9, ';', 4)[, 3]) %>%
  mutate(tx_id = str_remove_all(tx_id, 'Crovir_Transcript_ID=')) %>%
  filter(tx_id %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes, by = c('tx_id' = 'X6')) %>%
  dplyr::select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7, tx_id) %>%
  mutate(strand = ifelse(strand == '+', 'forward', 'reverse')) %>%
  mutate(direction = ifelse(strand == 'forward', 1, -1)) %>%
  mutate(gene = ifelse(str_detect(gene, 'ADAM28', negate = T), str_replace_all(gene, '_', ' '), gene)) %>%
  left_join(exp) %>%
  mutate(gene = ifelse(str_detect(gene, 'ADAM28'), paste('NVP: ', gene, sep = ''), gene)) %>%
  mutate(prom_start = ifelse(strand == 'forward', start, end))

# Filter for SVSP genes
SVSP_info <- all_info %>%
  filter(str_detect(gene, 'SVSP'))

# ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.

SVSP_reg_start = floor(8568727 / 1000 ) * 1000

SVSP_reg_end = ceiling(8981362 / 1000 ) * 1000

SVSP_reg_length <- paste(c(round((SVSP_reg_end-SVSP_reg_start)/1000, digits = 2), 'kb'), collapse = ' ')


# Get venom gene expression
VST_count_mat <- read.table('/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Data/mRNA/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)

# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>%
  dplyr::select(matches('LVG|RVG')) %>%
  filter(grepl('Venom|PLA2G2E.1', row.names(.)))

# Create a list of column numbers to operate over
column_numbers <- 1:13
# Use the above list to create a list of VG_i columns to operate over
target_names <- paste0('VG_', column_numbers)

# Calculate a pairwise mean for each sample
for (i in column_numbers) {
  target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))

  if (length(target_columns) > 0) {
    if (length(target_columns) == 1) {
      # If only one column is found, copy it to become VG_i
      VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
    } else {
      # Calculate row mean for multiple columns
      VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
    }
  }
}

VST_count_mat_filtered <- VST_count_mat %>%
  dplyr::select(-contains('LVG'), -contains('RVG')) %>%
  dplyr::select(contains('VG')) %>%
  dplyr::select('VG_5', 'VG_6', 'VG_7', 'VG_9', 'VG_2', 'VG_4') %>% # Filter for samples being used in study
  rownames_to_column(var = 'Gene') %>%
  filter(grepl('SVSP', Gene)) %>%
  mutate(AverageExp = rowMeans(dplyr::select(., starts_with("VG_")), na.rm = T)) %>%
  dplyr::select(Gene, AverageExp)

VST_count_mat_filtered <- VST_count_mat_filtered %>%
  mutate(Gene = gsub('Venom_SVSP', 'SVSP ', Gene)) %>%
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene))

# ADD TO SVSP INFO
SVSP_info <- SVSP_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
SVSP_info <- SVSP_info %>% mutate(gene = gsub(' ', '', gene))
all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
all_info <- all_info %>% mutate(gene = gsub(' ', '', gene))
SVSP_info <- SVSP_info %>%
  mutate(
    # Swap start and end if direction is -1
    start2 = case_when(direction == -1 ~ end,
                       T ~ start),
    end2 = case_when(direction == -1 ~ start,
                     T ~ end)
  )
SVSP_info$color = "blue"

# Plotting
SVSP_figure <- ggplot(SVSP_info, aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp + 1))) +
  ggrepel::geom_text_repel(data = all_info %>%
                             filter(str_detect(gene,'SVSP')) %>% # Change for SVSP-SVSP
                             mutate(start = (start + end)/2),
                           aes(x = start, y = str_remove(molecule, 'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2, size = 2.5) +
  geom_segment(aes(x = SVSP_reg_start, xend=SVSP_reg_end, y=str_remove(molecule,'scaffold\\-'), yend = str_remove(molecule, 'scaffold\\-')), lwd = 1, color = 'grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  labs(
    title = 'b. SVSPs'
  ) +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma, limits = c(SVSP_reg_start, SVSP_reg_end), expand = c(0,0)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    plot.title.position = 'plot',
    plot.title = element_text(color= 'black', hjust = 0.5, face= 'bold', size = 14),
    axis.title.x=element_blank()
  )
SVSP_figure
# Save SVSP_figure to a .pdf
ggsave("/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Figures/Gene_Arrays/SVSP_array_plot.2024.09.17.pdf", SVSP_figure, width = 10, height = 2, units = "in", create.dir = T)




#### SVMP Array plot ####

# Load packages
library(tidyverse)
library(rtracklayer)
library(scales)
library(gggenes)
library(viridis)
library(ggforce)

# CHANGE PATH UP TO /Dropbox
setwd('~/Dropbox/CastoeLabFolder/projects/z_Done_or_Abandoned_Projects-MSs/_VenomGeneRegulation_NEW_Aug2021')

# Read priority venom genes data
pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)

# Read gene expression data and filter for priority venom genes
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>%
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes, by=c('txid'='X6')) %>%
  mutate(gene = ifelse(str_detect(X7, 'ADAM28', negate = T), str_replace_all(X7, '_', ' '), X7))

# Read all gene information and filter for priority venom genes
all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff', col_names = F) %>%
  filter(str_detect(X9, 'trnascan', negate = T)) %>%
  mutate(tx_id = str_split_fixed(X9, ';', 4)[, 3]) %>%
  mutate(tx_id = str_remove_all(tx_id, 'Crovir_Transcript_ID=')) %>%
  filter(tx_id %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes, by = c('tx_id' = 'X6')) %>%
  dplyr::select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7, tx_id) %>%
  mutate(strand = ifelse(strand == '+', 'forward', 'reverse')) %>%
  mutate(direction = ifelse(strand == 'forward', 1, -1)) %>%
  mutate(gene = ifelse(str_detect(gene, 'ADAM28', negate = T), str_replace_all(gene, '_', ' '), gene)) %>%
  left_join(exp) %>%
  mutate(gene = ifelse(str_detect(gene, 'ADAM28'), paste('NVP: ', gene, sep = ''), gene)) %>%
  mutate(prom_start = ifelse(strand == 'forward', start, end))

# Filter for SVMP genes
SVMP.info <- all_info %>%
  filter(str_detect(gene, 'SVMP|ADAM'))

# ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions.

SVMP_reg_start = floor(13901005 / 1000 ) * 1000

SVMP_reg_end = ceiling(14424729 / 1000 ) * 1000

SVSP_reg_length <- paste(c(round((SVMP_reg_end-SVMP_reg_start)/1000,digits = 2),'kb'),collapse = ' ')


# Get venom gene expression
VST_count_mat <- read.table('/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Data/mRNA/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)

# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>%
  dplyr::select(matches('LVG|RVG')) %>%
  filter(grepl('Venom|PLA2G2E.1', row.names(.)))


# Create a list of column numbers to operate over
column_numbers <- 1:13
# Use the above list to create a list of VG_i columns to operate over
target_names <- paste0('VG_', column_numbers)

# Calculate a pairwise mean for each sample
for (i in column_numbers) {
  target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))
  
  if (length(target_columns) > 0) {
    if (length(target_columns) == 1) {
      # If only one column is found, copy it to become VG_i
      VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
    } else {
      # Calculate row mean for multiple columns
      VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = T)
    }
  }
}

VST_count_mat_filtered <- VST_count_mat %>%
  dplyr::select(-contains('LVG'), -contains('RVG')) %>%
  dplyr::select(contains('VG')) %>%
  dplyr::select('VG_5', 'VG_6', 'VG_7', 'VG_9', 'VG_2', 'VG_4') %>% # Filter for samples being used in study
  rownames_to_column(var = 'Gene') %>%
  filter(grepl('SVMP|ADAM', Gene)) %>%
  mutate(AverageExp = rowMeans(dplyr::select(., starts_with("VG_")), na.rm = T)) %>%
  dplyr::select(Gene, AverageExp)

VST_count_mat_filtered <- VST_count_mat_filtered %>%
  mutate(Gene = gsub('Venom_SVMP', 'SVMP ', Gene)) %>%
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene))

# ADD TO SVMP INFO
SVMP.info <- SVMP.info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
SVMP.info <- SVMP.info %>% mutate(gene = gsub(' ', '', gene))
all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
all_info <- all_info %>% mutate(gene = gsub(' ', '', gene))
SVMP.info <- SVMP.info %>%
  mutate(
    # Swap start and end if direction is -1
    start2 = case_when(direction == -1 ~ end,
                       T ~ start),
    end2 = case_when(direction == -1 ~ start,
                     T ~ end)
  )
SVMP.info$color = "blue"

# Plotting
SVMP_figure <- ggplot(SVMP.info,aes(xmin = start2, xmax = end2, y = str_remove(molecule,'scaffold\\-'), fill = log10(AverageExp + 1))) +
  ggrepel::geom_text_repel(data = all_info %>%
                             filter(str_detect(gene, 'SVMP')) %>% # Change for SVSP-SVMP
                             mutate(start = (start + end)/2),
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.2, size = 2.5) +
  geom_segment(aes(x = SVMP_reg_start, xend = SVMP_reg_end, y = str_remove(molecule, 'scaffold\\-'), yend = str_remove(molecule, 'scaffold\\-')), lwd = 1, color = 'grey70') +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  ylab('') +
  xlab('') +
  labs(
    title = 'c. SVMPs'
  ) +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma, limits = c(SVMP_reg_start, SVMP_reg_end), expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    plot.title.position = 'plot',
    plot.title = element_text(color= 'black', hjust = 0.5, face= 'bold', size = 14),
    axis.title.x=element_blank()
  )
SVMP_figure
# Save SVMP_figure to a .pdf
ggsave("/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Figures/Gene_Arrays/SVMP_array_plot.2024.09.17.pdf", SVMP_figure, width = 10, height = 2, units = "in", create.dir = T)



# Create a grid for these figures
gene_arrays <- plot_grid(
  PLA2_figure,
  SVSP_figure,
  SVMP_figure,
  ncol = 1,
  nrow = 3
)
gene_arrays
ggsave("/Users/ballardk/Library/CloudStorage/Dropbox/CastoeLabFolder/projects/Venom_Gene_Regulation/Venom_ncRNA/Figures/Gene_Arrays/Three_families_array.2024.09.17.pdf", gene_arrays, width = 10, height = 10, units = "in", create.dir = T)



