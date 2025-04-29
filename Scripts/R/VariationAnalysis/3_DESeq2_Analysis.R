# Last Edited: 2025/02/05

# Set up and Read Data ----

## Load in packages ----
library(compositions)
library(cowplot)
library(tidyverse)
library(arrow)
library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(readxl)
library(ggpmisc)
library(ggplot2)

## Load in DESeq2 ----
library(DESeq2)
library(regionReport)
library("IHW")
library('ashr')
library('vsn')
library('pheatmap')
# Import counts plotting function
source('Scripts/R/Functions/Plot_Gene_Counts_Function.R')


## Read the miRNA data
# Create variable for the fused data set.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Create a variable for the mRNA data
mRNA_protein_data <- 'Data/Merged/mRNA_Protein_Combined_Data_2025.01.22.parquet'

# Create a variable for the miRNA data
miRNA_data <- 'Data/miRNA/Proccessed/miRNA_counts_for_reference/miRNA_counts-rpm-vst-rlog-norm_counts.2025.01.22.parquet'

# Establish the path for the sample metadata
metadata <- 'Data/metadata/DESeq_Sample_metadata_2024.6.4.csv'

## Create a relationships DataFrame ----
relationships_df <- read_parquet(file = miRNA_mRNA_protein_data) %>% 
  select(
    genes, venom.family, miRNA.cluster, miRNA.cluster.original, feature.type, total.score, total.energy
  ) %>% 
  # Remove badly annotated genes
  dplyr::filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_')
  ) %>% 
  distinct()

## miRNA data frame ----
mirna_df <- read_parquet(file = miRNA_data) %>% 
  dplyr::select(
    sample.id, miRNA.cluster, miRNA.counts
  ) %>% 
  dplyr::distinct() %>% 
  pivot_wider(
    names_from = sample.id,
    values_from = miRNA.counts
  ) %>% 
  dplyr::rename(genes = miRNA.cluster)
rownames(mirna_df) <- NULL # Reset row names

## mRNA data frame ----
mrna_df <- read_parquet(file = mRNA_protein_data) %>% 
  dplyr::select(
    sample.id, genes, mRNA.counts
  ) %>%
  dplyr::distinct() %>%
  # Remove badly annotated genes
  dplyr::filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_')
  ) %>% 
  pivot_wider(
    names_from = sample.id,
    values_from = mRNA.counts
  )
rownames(mrna_df) <- NULL

## Read the sample data ----
sample_data <- read.csv(metadata, row.names = 1)
sample_data$species <- factor(sample_data$species) # Turn the sample data into factors
sample_data$sub_species <- factor(sample_data$sub_species) # Turn the sample data into factors
sample_data$state <- factor(sample_data$state)
sample_data$locality <- factor(sample_data$locality)
sample_data$sex <- factor(sample_data$sex)
sample_data$tissue <- factor(sample_data$tissue)

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


# DESeq2 Analysis ----

# Set species colors
viridis_color <- '#2B8B59'
oreganus_color <- '#D8A628'

# Concatenate the data frames
counts_matrix <- bind_rows(mrna_df, mirna_df) %>% 
  # Re-order to be in the same order as in the sample data
  select(genes, CV1081_viridis, CV0857_viridis, CV1086_viridis, CV1082_viridis, CV1087_viridis, CV0987_lutosus, CV0985_concolor)
glimpse(counts_matrix)

# Turn genes column into rownames
counts_matrix <- as.data.frame(counts_matrix)
rownames(counts_matrix) <- counts_matrix$genes
counts_matrix$genes <- NULL  # Remove the genes column after setting row names



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
  outdir = 'Figures/DESeq2/P-adjusted/Region_Report/DESeq2_miRNA-mRNA-Region_Report',
  output = 'Species_vs_Species'
)
browseURL(report1)


# Extract results for miRNAs and mRNAs separately
res_combined_df <- as.data.frame(res)
res_combined_df$genes <- rownames(res_combined_df) # Create new column based off of rownames


# Identify miRNA and mRNA results based on prefixes
# Filter for any miRNAs
res_mirna <- res_combined_df %>% filter(grepl("^cvi-|^Cluster_", genes)) %>% dplyr::rename(miRNA.cluster = genes)
rownames(res_mirna) <- NULL
# Filter out any miRNAs
res_mrna <- res_combined_df %>% filter(!grepl("^cvi-|^Cluster_", genes))
rownames(res_mrna) <- NULL

# Create a data frame for interaction data
interactions_df <- relationships_df %>% dplyr::select(miRNA.cluster, genes) %>% distinct()

# Join results into the same data frame again
merged_diff_df <- interactions_df %>% 
  left_join(res_mirna, by = 'miRNA.cluster') %>% 
  left_join(res_mrna, by = 'genes') %>% 
  dplyr::rename(
    baseMean.mi = baseMean.x, log2FoldChange.mi = log2FoldChange.x, lfcSE.mi = lfcSE.x, stat.mi = stat.x, pvalue.mi = pvalue.x, padj.mi = padj.x,
    baseMean.mr = baseMean.y, log2FoldChange.mr = log2FoldChange.y, lfcSE.mr = lfcSE.y, stat.mr = stat.y, pvalue.mr = pvalue.y, padj.mr = padj.y
  ) %>% 
  filter(!str_detect(genes, 'maker-scaffold|augustus|XP_')) # This should get rid of all of the weirdly annotated genes, make sure to note you did this in the methods section

# Filter down the above to only significant pairs
significant_diff_paired_df <- merged_diff_df %>% 
  filter(padj.mi < 0.05 & padj.mr < 0.05) %>% 
  filter(!is.na(padj.mi) & !is.na(padj.mr))
# Filter for each individually as well
significant_mi_diff_df <- merged_diff_df %>% filter(padj.mi < 0.05) %>% filter(!is.na(padj.mi))
significant_mr_diff_df <- merged_diff_df %>% filter(padj.mr < 0.05) %>% filter(!is.na(padj.mr))

# Extract normalized counts for significant pairs
significant_genes <- unique(c(significant_diff_paired_df$miRNA.cluster, significant_diff_paired_df$genes))
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

## PCA ----
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
  data = pca_df, aes(PC1, PC2, color = species, label = name)
) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps =  20) +
  xlab(paste0('PC1: ', percentVar[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar[2], '% variance')) +
  coord_fixed() +
  labs(title = 'PCA for differential expression between species') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) # Center the title
  # stat_ellipse(level = 0.95)
pca1
ggsave('Figures/DESeq2/P-adjusted/PCA/Species_DESeq_PCA.2025.02.05.pdf', pca1, create.dir = T)


## DESeq for miRNAs ----
# Concatenate the data frames
mirna_counts <- mirna_df %>% 
  dplyr::rename(miRNA.cluster = genes) %>%
  # Re-order to be in the same order as in the sample data
  select(miRNA.cluster, CV1081_viridis, CV0857_viridis, CV1086_viridis, CV1082_viridis, CV1087_viridis, CV0987_lutosus, CV0985_concolor)

# Turn genes column into rownames
mirna_counts <- as.data.frame(mirna_counts)
rownames(mirna_counts) <- mirna_counts$miRNA.cluster
mirna_counts$miRNA.cluster <- NULL  # Remove the genes column after setting row names



# Create deseq data set for miRNAs
dds_mirna <- DESeqDataSetFromMatrix(
  countData = mirna_counts,
  colData = sample_data,
  design = ~ species
)
dds_mirna

# Run DESeq for the miRNA data
dds_mirna <- DESeq(dds_mirna)
res_mirna2 <- results(dds_mirna)

# Trying an MA-plot on both
plotMA(res_mirna2, ylim = c(-2,2))
# Nothing really stands out

# Let's create reports for both
# Create report to do some initial analysis
report_mi <- DESeq2Report(
  dds = dds_mirna,
  res = res_mirna2,
  intgroup = c('species'),
  project = 'miRNA-mRNA Analysis',
  author = 'Kaas Ballard',
  outdir = 'Figures/DESeq2/P-adjusted/Region_Report/DESeq2_miRNA-Region_Report',
  output = 'Species_vs_Species'
)
browseURL(report_mi)

# Extract results for miRNAs and mRNAs separately
res_mi_df <- as.data.frame(res_mirna2)
res_mi_df$miRNA.cluster <- rownames(res_mi_df) # Create new column based off of rownames

### PCA ----

# Transform and extracting transformed values
vsd_mirna <- varianceStabilizingTransformation(dds_mirna, blind = FALSE)
# I have to use this because it gives an error with the vst function

# Plot different interactions
plotPCA(vsd_mirna, intgroup = c("sub_species", "sex"))
plotPCA(vsd_mirna, intgroup = c('species'))

# Extract data for ggplot PCA
pca_mi_df <- plotPCA(vsd_mirna, intgroup = c('species'), returnData = T)

# Create percent variation numbers
percentVar_mirna <- round(100 * attr(pca_mi_df, 'percentVar'))

# Plot using ggplot
pca_miRNA <- ggplot(
  data = pca_mi_df, aes(PC1, PC2, color = species, label = name)
) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps =  20) +
  xlab(paste0('PC1: ', percentVar_mirna[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar_mirna[2], '% variance')) +
  coord_fixed() +
  labs(title = 'PCA for differential expression between species') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title
  # stat_ellipse(level = 0.95)
pca_miRNA
ggsave('Figures/DESeq2/P-adjusted/PCA/Species_DESeq_PCA_miRNA_only.2025.02.05.pdf', pca_miRNA, create.dir = T)


## DESeq for mRNAs ----
# Concatenate the data frames
mrna_counts <- mrna_df %>% 
  # Re-order to be in the same order as in the sample data
  select(genes, CV1081_viridis, CV0857_viridis, CV1086_viridis, CV1082_viridis, CV1087_viridis, CV0987_lutosus, CV0985_concolor)

# Turn genes column into rownames
mrna_counts <- as.data.frame(mrna_counts)
rownames(mrna_counts) <- mrna_counts$genes
mrna_counts$genes <- NULL  # Remove the genes column after setting row names



# Create deseq data set for miRNAs
dds_mrna <- DESeqDataSetFromMatrix(
  countData = mrna_counts,
  colData = sample_data,
  design = ~ species
)
dds_mrna

# Run DESeq for the miRNA data
dds_mrna <- DESeq(dds_mrna)
res_mrna2 <- results(dds_mrna)

# Trying an MA-plot on both
plotMA(res_mrna2, ylim = c(-2,2))
# Nothing really stands out

# Let's create reports for both
# Create report to do some initial analysis
report_mr <- DESeq2Report(
  dds = dds_mrna,
  res = res_mrna2,
  intgroup = c('species'),
  project = 'miRNA-mRNA Analysis',
  author = 'Kaas Ballard',
  outdir = 'Figures/DESeq2/P-adjusted/Region_Report/DESeq2_mRNA-Region_Report',
  output = 'Species_vs_Species'
)
browseURL(report_mr)

# Extract results for miRNAs and mRNAs separately
res_mr_df <- as.data.frame(res_mrna2)
res_mr_df$genes <- rownames(res_mr_df) # Create new column based off of rownames

### PCA ----

# Transform and extracting transformed values
vsd_mrna <- vst(dds_mrna, blind = FALSE)
# I have to use this because it gives an error with the vst function

# Plot different interactions
plotPCA(vsd_mrna, intgroup = c("sub_species", "sex"))
plotPCA(vsd_mrna, intgroup = c('species'))

# Extract data for ggplot PCA
pca_mr_df <- plotPCA(vsd_mrna, intgroup = c('species'), returnData = T)

# Create percent variation numbers
percentVar_mrna <- round(100 * attr(pca_mr_df, 'percentVar'))

# Plot using ggplot
pca_mRNA <- ggplot(
  data = pca_mr_df, aes(PC1, PC2, color = species, label = name)
) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps =  20) +
  xlab(paste0('PC1: ', percentVar_mrna[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar_mrna[2], '% variance')) +
  coord_fixed() +
  labs(title = 'PCA for differential expression between species') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title
# stat_ellipse(level = 0.95)
pca_mRNA
ggsave('Figures/DESeq2/P-adjusted/PCA/Species_DESeq_PCA_mRNA_only.2025.02.05.pdf', pca_mRNA, create.dir = T)

## DESeq for Venom genes ----
# Concatenate the data frames
venom_counts <- mrna_df %>% 
  dplyr::filter(str_detect(genes, 'Venom_')) %>% 
  # Re-order to be in the same order as in the sample data
  select(genes, CV1081_viridis, CV0857_viridis, CV1086_viridis, CV1082_viridis, CV1087_viridis, CV0987_lutosus, CV0985_concolor)

# Turn genes column into rownames
venom_counts <- as.data.frame(venom_counts)
rownames(venom_counts) <- venom_counts$genes
venom_counts$genes <- NULL  # Remove the genes column after setting row names



# Create deseq data set for miRNAs
dds_venom <- DESeqDataSetFromMatrix(
  countData = venom_counts,
  colData = sample_data,
  design = ~ species
)
dds_venom

# Run DESeq for the miRNA data
dds_venom <- DESeq(dds_venom)
res_venom <- results(dds_venom)

# Trying an MA-plot on both
plotMA(res_venom, ylim = c(-2,2))
# Nothing really stands out

# Let's create reports for both
# Create report to do some initial analysis
report_venom <- DESeq2Report(
  dds = dds_venom,
  res = res_venom,
  intgroup = c('species'),
  project = 'miRNA-mRNA Analysis',
  author = 'Kaas Ballard',
  outdir = 'Figures/DESeq2/P-adjusted/Region_report/DESeq2_Venom-Region_Report',
  output = 'Species_vs_Species'
)
browseURL(report_venom)

### PCA ----

# Transform and extracting transformed values
vsd_venom <- varianceStabilizingTransformation(dds_venom, blind = FALSE)
# I have to use this because it gives an error with the vst function

# Plot different interactions
plotPCA(vsd_venom, intgroup = c("sub_species", "sex"))
plotPCA(vsd_venom, intgroup = c('species'))

# Extract data for ggplot PCA
pca_venom_df <- plotPCA(vsd_venom, intgroup = c('species'), returnData = T)

# Create percent variation numbers
percentVar_venom <- round(100 * attr(pca_venom_df, 'percentVar'))

# Plot using ggplot
pca_mRNA <- ggplot(
  data = pca_venom_df, aes(PC1, PC2, color = species, label = name)
) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps =  20) +
  xlab(paste0('PC1: ', percentVar_venom[1], '% variance')) +
  ylab(paste0('PC2: ', percentVar_venom[2], '% variance')) +
  coord_fixed() +
  labs(title = 'PCA for differential expression between species') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title
# stat_ellipse(level = 0.95)
pca_mRNA
ggsave('Figures/DESeq2/P-adjusted/PCA/Species_DESeq_PCA_venom_genes_only.2025.02.05.pdf', pca_mRNA, create.dir = T)




# Gene family plots ----

## SVMP test ----
plotCounts(dds, gene = 'Venom_SVMP1', intgroup = 'species')
# Extract data
svmp1_df <- plotCounts(dds, gene = 'Venom_SVMP1', intgroup = 'species', returnData = T)
svmp1_plot <- ggplot(svmp1_df, aes(x = species, y = count)) +
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  scale_y_log10(breaks = c(25, 100, 400)) +
  labs(y = 'Normalized Counts', x = 'Species', title = 'SVMP1 Differential Expression') +
  theme_classic2() +
  theme(plot.title = element_text(hjust = 0.5))
svmp1_plot



## DESeq Count Plots For Venom genes and Targeting miRNAs ----

# Create a merged diff for only venom genes
venom_merged_diff_df <- merged_diff_df %>% 
  filter(str_detect(genes, 'Venom_'))

# Create a new df for the genes only
venom_genes_only_df <- venom_merged_diff_df %>% 
  dplyr::select(-miRNA.cluster, -contains('.mi')) %>% 
  distinct()

# Create a new df for the miRNAs only
venom_mirna_only_df <- venom_merged_diff_df %>% 
  dplyr::select(-genes, -contains('.mr')) %>% 
  distinct()

# Create a list of genes to be run the function on based on the merged data frame
genes <- venom_genes_only_df$genes

# Run the function in a loop
for (gene in genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Venom/genes/Significant_and_Nonsignificant', date = '2025.02.05')
}

# Create a list of miRNAs
miRNAs <- venom_mirna_only_df$miRNA.cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Venom/miRNAs/Significant_and_Nonsignificant', date = '2025.02.05')
}


## DESeq Count Plots For Significant Venom genes and Targeting miRNA pairs ----

# Create a data frame of significant venom genes
venom_significant_diff_paired_df <- significant_diff_paired_df %>% 
  filter(str_detect(genes, 'Venom_')) %>% 
  distinct()

# Create a new df for the genes only
venom_sig_genes_paired_df <- venom_significant_diff_paired_df %>% 
  dplyr::select(-miRNA.cluster, -contains('.mi')) %>% 
  distinct()

# Create a new df for the miRNAs only
venom_sig_mi_paired_df <- venom_significant_diff_paired_df %>% 
  dplyr::select(-genes, -contains('.mr')) %>% 
  distinct()

# Create a list of genes that are significantly differentially expressed
sig_genes <- venom_sig_genes_paired_df$genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Venom/genes/Significant_Paired', date = '2025.02.05')
}

# Create a list of miRNAs
sig_miRNAs <- venom_sig_mi_paired_df$miRNA.cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Venom/miRNAs/Significant_Paired', date = '2025.02.05')
}


## DESeq Count Plots For Significant Venom genes and Targeting miRNAs that are not paired ----

# Create a data frame of significant venom genes (not paired)
venom_sig_mr_diff_df <- significant_mr_diff_df %>% 
  filter(str_detect(genes, 'Venom_')) %>% 
  dplyr::select(-miRNA.cluster, -contains('.mi')) %>% 
  distinct()

# Create a new df for the miRNAs only
venom_sig_mi_diff_df <- significant_mi_diff_df %>% 
  filter(str_detect(genes, 'Venom_')) %>% 
  dplyr::select(-genes, -contains('.mr')) %>% 
  distinct()


# Create a list of genes that are significantly differentially expressed
sig_genes <- venom_sig_mr_diff_df$genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Venom/genes/Significant', date = '2025.02.05')
}

# Create a list of miRNAs
sig_miRNAs <- venom_sig_mi_diff_df$miRNA.cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Venom/miRNAs/Significant', date = '2025.02.05')
}



## DESeq Count Plots For Non-Venom genes and Targeting miRNAs ----

# Create a merged diff for only nonvenom genes
non_venom_merged_diff_df <- merged_diff_df %>% 
  filter(!str_detect(genes, 'Venom_')) %>%
  distinct()

# Create a new df for the miRNAs only
non_venom_mirna_df <- non_venom_merged_diff_df %>% 
  dplyr::select(-genes, -contains('.mr')) %>% 
  filter(!is.na(miRNA.cluster)) %>% 
  distinct()

# I don't actually care about this, so I won't run it because it would take to long
# # Create a list of genes to be run the function on based on the merged data frame
# genes <- non_venom_genes_df$genes
# 
# # Run the function in a loop
# for (gene in genes) {
#   plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/New_Protein_Method/3UTR/DESeq_Counts/Species/Non_Venom/genes/Significant_and_Nonsignificant', date = '2024.08.02')
# }

# Create a list of miRNAs
miRNAs <- non_venom_mirna_df$miRNA.cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Non_Venom/miRNAs/Significant_and_Nonsignificant', date = '2025.02.05')
}



# Create a data frame of significant non-venom genes
non_venom_significant_diff_paired_df <- significant_diff_paired_df %>% 
  filter(!str_detect(genes, 'Venom_')) %>% 
  distinct()

# Create a new df for the genes only
non_venom_sig_genes_paired_df <- non_venom_significant_diff_paired_df %>% 
  dplyr::select(-miRNA.cluster, -contains('.mi')) %>% 
  filter(!str_detect(genes, 'maker-scaffold')) %>% 
  distinct()

# Create a new df for the miRNAs only
non_venom_sig_mi_paired_df <- non_venom_significant_diff_paired_df %>% 
  dplyr::select(-genes, -contains('.mr')) %>% 
  distinct()


# Create a list of genes that are significantly differentially expressed
sig_genes <- non_venom_sig_genes_paired_df$genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Non_Venom/genes/Significant_Paired', date = '2025.02.05')
}

# Create a list of miRNAs
sig_miRNAs <- non_venom_sig_mi_paired_df$miRNA.cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Non_Venom/miRNAs/Significant_Paired', date = '2025.02.05')
}


## DESeq Count Plots For Significant Non-Venom genes and Targeting miRNAs that are not paired ----

# Create a data frame of significant venom genes (not paired)
non_venom_sig_mr_diff_df <- significant_mr_diff_df %>% 
  filter(!str_detect(genes, 'Venom_')) %>% 
  dplyr::select(-miRNA.cluster, -contains('.mi')) %>% 
  distinct()

# Create a new df for the miRNAs only
non_venom_sig_mi_diff_df <- significant_mi_diff_df %>% 
  filter(!str_detect(genes, 'Venom_|PLA2GE.1')) %>% 
  dplyr::select(-genes, -contains('.mr')) %>% 
  distinct()


# Create a list of genes that are significantly differentially expressed
sig_genes <- non_venom_sig_mr_diff_df$genes

# Run the function in a loop
for (gene in sig_genes) {
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Non_Venom/genes/Significant', date = '2025.02.05')
}

# Create a list of miRNAs
sig_miRNAs <- non_venom_sig_mi_diff_df$miRNA.cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/P-adjusted/DESeq_Counts/Species/Non_Venom/miRNAs/Significant', date = '2025.02.05')
}



## Create nicer figures for significant miRNAs and genes ----

### SVSP7 Group ----
#### SVSP7 ----
# Data for SVSP7
svsp7 <- plotCounts(dds, gene = 'Venom_SVSP7', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
svsp7_stats <- svsp7 %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(svsp7_stats$color, svsp7_stats$species)

# Plot the gene
svsp7_plot <- ggplot(svsp7_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_fill_manual(values = colors) +
  labs(
    title = 'SVSP7 normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
svsp7_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVSP7_Group/SVSP7_DESeq2_normalized_counts_2025.02.05.pdf', plot = svsp7_plot, width = 4, height = 5, dpi = 900, create.dir = T)


#### cvi-miR-9-5p ----
# Data for cvi-miR-9-5p
cvi_miR_9_5p <- plotCounts(dds, gene = 'cvi-miR-9-5p', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cvi_miR_9_5p_stats <- cvi_miR_9_5p %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cvi_miR_9_5p_stats$color, cvi_miR_9_5p_stats$species)

# Plot the gene
cvi_miR_9_5p_plot <- ggplot(cvi_miR_9_5p_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'cvi-miR-9-5p normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cvi_miR_9_5p_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVSP7_Group/cvi-miR-9-5p_DESeq2_Targets_SVSP7_normalized_counts_2025.02.05.pdf', plot = cvi_miR_9_5p_plot, width = 4, height = 5, dpi = 900, create.dir = T)


#### Cluster_853 ----
# Data for Cluster_853
cluster_853 <- plotCounts(dds, gene = 'Cluster_853', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cluster_853_stats <- cluster_853 %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cluster_853_stats$color, cluster_853_stats$species)

# Plot the gene
cluster_853_plot <- ggplot(cluster_853_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'Cluster_853 normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cluster_853_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVSP7_Group/Cluster_853_DESeq2_Targets_SVSP7_normalized_counts_2025.02.05.pdf', plot = cluster_853_plot, width = 4, height = 5, dpi = 900, create.dir = T)




### SVMP7 Group ----
#### SVMP7 ----
# Data for SVMP7
svmp7 <- plotCounts(dds, gene = 'Venom_SVMP7', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
svmp7_stats <- svmp7 %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(svmp7_stats$color, svmp7_stats$species)

# Plot the gene
svmp7_plot <- ggplot(svmp7_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_fill_manual(values = colors) +
  labs(
    title = 'SVMP7 normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
svmp7_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVMP7_Group/SVMP7_DESeq2_normalized_counts_2025.02.05.pdf', plot = svmp7_plot, width = 4, height = 5, dpi = 900, create.dir = T)

#### Cluster_1718 ----
# Data for Cluster_1718
cluster_1718 <- plotCounts(dds, gene = 'Cluster_1718', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cluster_1718_stats <- cluster_1718 %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cluster_1718_stats$color, cluster_1718_stats$species)

# Plot the gene
cluster_1718_plot <- ggplot(cluster_1718_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'Cluster_1718 normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cluster_1718_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVMP7_Group/Cluster_1718_DESeq2_Targets_SVMP7_normalized_counts_2025.02.05.pdf', plot = cluster_1718_plot, width = 4, height = 5, dpi = 900, create.dir = T)


#### cvi-miR-24-2-5p ----
# Data for cvi-miR-24-2-5
cvi_miR_24_2_5p <- plotCounts(dds, gene = 'cvi-miR-24-2-5p', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cvi_miR_24_2_5p_stats <- cvi_miR_24_2_5p %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cvi_miR_24_2_5p_stats$color, cvi_miR_24_2_5p_stats$species)

# Plot the gene
cvi_miR_24_2_5p_plot <- ggplot(cvi_miR_24_2_5p_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'cvi-miR-24-2-5p normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cvi_miR_24_2_5p_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVMP7_Group/cvi-miR-24-2-5p_DESeq2_Targets_SVMP7_normalized_counts_2025.02.05.pdf', plot = cvi_miR_24_2_5p_plot, width = 4, height = 5, dpi = 900, create.dir = T)



#### cvi-miR-215-3p ----
# Data for cvi-miR-24-2-5
cvi_miR_215_3p <- plotCounts(dds, gene = 'cvi-miR-215-3p', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cvi_miR_215_3p_stats <- cvi_miR_215_3p %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cvi_miR_215_3p_stats$color, cvi_miR_215_3p_stats$species)

# Plot the gene
cvi_miR_215_3p_plot <- ggplot(cvi_miR_215_3p_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'cvi-miR-215-3p normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cvi_miR_215_3p_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVMP7_Group/cvi-miR-215-3p_DESeq2_Targets_SVMP7_normalized_counts_2025.02.05.pdf', plot = cvi_miR_215_3p_plot, width = 4, height = 5, dpi = 900, create.dir = T)


#### cvi-miR-34a-5p ----
# Data for cvi-miR-24-2-5
cvi_miR_34a_5p <- plotCounts(dds, gene = 'cvi-miR-34a-5p', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cvi_miR_34a_5p_stats <- cvi_miR_34a_5p %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cvi_miR_34a_5p_stats$color, cvi_miR_34a_5p_stats$species)

# Plot the gene
cvi_miR_34a_5p_plot <- ggplot(cvi_miR_34a_5p_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'cvi-miR-34a-5p normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cvi_miR_34a_5p_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVMP7_Group/cvi-miR-34a-5p_DESeq2_Targets_SVMP7_normalized_counts_2025.02.05.pdf', plot = cvi_miR_34a_5p_plot, width = 4, height = 5, dpi = 900, create.dir = T)




### Both SVMP7 and SVSP7 ----
#### Cluster_1395 ----
# Data for Cluster_1395
cluster_1395 <- plotCounts(dds, gene = 'Cluster_1395', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cluster_1395_stats <- cluster_1395 %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cluster_1395_stats$color, cluster_1395_stats$species)

# Plot the gene
cluster_1395_plot <- ggplot(cluster_1395_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'Cluster_1395 normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cluster_1395_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVMP7_Group/Cluster_1395_DESeq2_Targets_SVSP7_SVMP7_normalized_counts_2025.02.05.pdf', plot = cluster_1395_plot, width = 4, height = 5, dpi = 900, create.dir = T)


### CRISP2 group ----
#### CRISP 2 ----
# Data for CRISP2
crisp2 <- plotCounts(dds, gene = 'Venom_CRISP2', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
crisp2_stats <- crisp2 %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(crisp2_stats$color, crisp2_stats$species)

# Plot the gene
crisp2_plot <- ggplot(crisp2_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_fill_manual(values = colors) +
  labs(
    title = 'CRISP2 normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
crisp2_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/CRISP2_Group/CRISP2_DESeq2_normalized_counts_2025.02.05.pdf', plot = crisp2_plot, width = 4, height = 5, dpi = 900, create.dir = T)



#### cvi-miR-34a-5p ----
# Data for cvi-miR-34a-5p
cvi_mir_34a_5p <- plotCounts(dds, gene = 'cvi-miR-34a-5p', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
cvi_mir_34a_5p_stats <- cvi_mir_34a_5p %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(cvi_mir_34a_5p_stats$color, cvi_mir_34a_5p_stats$species)

# Plot the gene
cvi_mir_34a_5p_plot <- ggplot(cvi_mir_34a_5p_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_fill_manual(values = colors) +
  labs(
    title = 'cvi-miR-34a-5p normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
cvi_mir_34a_5p_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/CRISP2_Group/cvi-miR-34a-5p_DESeq2_Targets_CRISP2_normalized_counts_2025.02.05.pdf', plot = cvi_mir_34a_5p_plot, width = 4, height = 5, dpi = 900, create.dir = T)

### PLA2K Group ----
#### PLA2K ----
# Data for cvi-miR-34a-5p
pla2k <- plotCounts(dds, gene = 'Venom_PLA2K', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
pla2k_stats <- pla2k %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(pla2k_stats$color, pla2k_stats$species)

# Plot the gene
pla2k_plot <- ggplot(pla2k_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_fill_manual(values = colors) +
  labs(
    title = 'PLA2K normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
pla2k_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/PLA2K_Group/PLA2K_DESeq2_normalized_counts_2025.02.05.pdf', plot = pla2k_plot, width = 4, height = 5, dpi = 900, create.dir = T)


#### Cluster_182 ----
# Data for cvi-miR-24-2-5
Cluster_182 <- plotCounts(dds, gene = 'Cluster_182', intgroup = 'species', returnData = T)

# Calculate mean and standard deviation for each group
Cluster_182_stats <- Cluster_182 %>% 
  group_by(species) %>% 
  summarise(
    mean_count = mean(count),
    sd_count = sd(count),
    n = n(),
    se_count = sd_count / sqrt(n) # Standard error for variance intervals
  ) %>% 
  mutate(
    color = case_when(
      grepl('Crotalus_oreganus', species) ~ oreganus_color,
      grepl('Crotalus_viridis', species) ~ viridis_color
    )
  ) %>% 
  mutate(
    species = case_when(
      grepl('Crotalus_oreganus', species) ~ 'C. oreganus',
      grepl('Crotalus_viridis', species) ~ 'C. viridis'
    )
  )

# Set a vector for the colors
colors <- setNames(Cluster_182_stats$color, Cluster_182_stats$species)

# Plot the gene
Cluster_182_plot <- ggplot(Cluster_182_stats, aes(x = species, y = mean_count, fill = species)) +
  geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) + 
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
  scale_y_continuous() +
  scale_fill_manual(values = colors) +
  labs(
    title = 'Cluster_182 normalized counts',
    x = 'Species',
    y = 'Mean normalized counts'
  ) +
  theme_classic2() +
  theme(
    plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.position = 'none'
  )
Cluster_182_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/PLA2K_Group/Cluster_182_DESeq2_Targets_PLA2_normalized_counts_2025.02.05.pdf', plot = Cluster_182_plot, width = 4, height = 5, dpi = 900, create.dir = T)




### Combine Plots ----
# Combine the plots
svmp7_and_miRNA <- plot_grid(
  svmp7_plot,
  cluster_1395_plot,
  cvi_miR_215_3p_plot,
  cvi_miR_34a_5p_plot,
  cvi_miR_24_2_5p_plot,
  cluster_1395_plot,
  cluster_1718_plot,
  ncol = 3,
  nrow = 3,
  align = 'v'
)
svmp7_and_miRNA
ggsave('Figures/1_Main_Figures/Figure_z/Figure_z_2025.02.05.pdf', plot = svmp7_and_miRNA, width = 8, height = 12, dpi = 900, create.dir = T)

svsp7_and_miRNA <- plot_grid(
  svsp7_plot,
  cluster_1395_plot,
  cluster_853_plot,
  cvi_miR_9_5p_plot,
  nrow = 2,
  ncol = 2,
  align = 'v'
)
svsp7_and_miRNA
ggsave('Figures/1_Main_Figures/Figure_ab/Figure_ab_2025.02.05.pdf', plot = svsp7_and_miRNA, width = 5.5, height = 8, dpi = 900, create.dir = T)


crisp2_and_miRNA <- plot_grid(
  crisp2_plot,
  cvi_miR_34a_5p_plot,
  nrow = 1,
  ncol = 2,
  align = 'v'
)
crisp2_and_miRNA
ggsave('Figures/1_Main_Figures/Figure_ac/Figure_ac_2025.02.05.pdf', plot = crisp2_and_miRNA, width = 5.5, height = 4, dpi = 900, create.dir = T)

pla2k_and_miRNA <- plot_grid(
  pla2k_plot,
  Cluster_182_plot,
  nrow = 1,
  ncol = 2,
  align = 'v'
)
pla2k_and_miRNA
ggsave('Figures/1_Main_Figures/Figure_ad/Figure_ad_2025.02.05.pdf', plot = pla2k_and_miRNA, width = 5, height = 4, dpi = 900, create.dir = T)


### Create Plots Containing Venom Genes and their paired miRNAs (sig DE) ----

#### Format data ----
# Get VST
vst <- as.data.frame(assay(vsd))

# Convert row names to column
vst <- rownames_to_column(vst, var = 'genes') %>% 
  pivot_longer(
    cols = matches('CV'),
    names_to = 'sample.id',
    values_to = 'vst'
  )

# Convert the dds results to a more managable format and filter out insignificant genes
dds_results <- remove_rownames(res_combined_df) %>% 
  filter(
    # Filter out insignificant results
    padj < 0.05,
    # Filter out anything that isn't Venom or miRNA
    str_detect(genes, 'Venom_|cvi-|Cluster')
  ) %>% 
  # Join the vst results to this
  left_join(
    vst,
    by = 'genes'
  )

# Get the significant miRNAs
sig_miRNA_df <- dds_results %>% 
  filter(!str_detect(genes, 'Venom')) %>% 
  rename('genes' = 'miRNA.cluster') %>% 
  rename_with(
    ~paste0('mi.', .), -c('miRNA.cluster', 'sample.id')
  )

# Get the significant venoms
sig_genes_df <- dds_results %>% 
  filter(str_detect(genes, 'Venom')) %>% 
  rename_with(
    ~paste0('mr.', .), -c('genes', 'sample.id')
  )

# Fuse back together with the relationships_df at the core
sig_df <- relationships_df %>% 
  full_join(
    sig_genes_df,
    by = c('genes'),
    relationship = 'many-to-many'
  ) %>% 
  filter(!is.na(mr.vst)) %>%
  inner_join(
    sig_miRNA_df,
    by = c('miRNA.cluster', 'sample.id')
  ) %>% 
  distinct(
    sample.id, genes, venom.family, miRNA.cluster, mr.pvalue, mi.pvalue, mr.vst, mi.vst,
  ) %>% 
  mutate(
    species = case_when(
      str_detect(sample.id, 'lutosus|concolor') ~ 'C. oreganus group',
      TRUE ~ 'C. v. viridis'
    )
  )

# Create a vector for the colors
species_colors <- c(
  'C. oreganus group' = oreganus_color,
  'C. v. viridis' = viridis_color
)

#### Figures ----

##### CRISP2 ----

# Create a data frame for the gene in quesion
crisp2_df <- sig_df %>% 
  filter(str_detect(genes, 'CRISP'))

# Get it's miRNAs
mi_crisp2_df <- crisp2_df %>% 
  select(sample.id, genes = miRNA.cluster, pvalue = mi.pvalue, vst = mi.vst, species)

# Format the data frame for joining
crisp2_df <- crisp2_df %>% 
  distinct(sample.id, genes, pvalue = mr.pvalue, vst = mr.vst, species) %>% 
  full_join(mi_crisp2_df)
  
crisp2_summary_df <- crisp2_df %>% 
  group_by(genes, species) %>%
  summarise(
    mean_vst = mean(vst),
    sd_vst = sd(vst),
    se_vst = sd(vst) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Create a plot combining bar chart and error bars with individual points
crisp2_bar_plot <- ggplot(crisp2_summary_df, aes(x = species, y = mean_vst, fill = species)) +
  # Add bars for the means
  geom_bar(stat = "identity", width = 0.7) +
  # Add error bars for standard error
  geom_errorbar(aes(ymin = mean_vst - se_vst, ymax = mean_vst + se_vst), 
                width = 0.2, position = position_dodge(0.7)) +
  # Add individual data points
  geom_point(data = crisp2_df, aes(x = species, y = vst, group = species),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
             size = 2, alpha = 0.8, shape = 16) +
  # Use your color scheme
  scale_fill_manual(values = species_colors) +
  # Labels
  labs(
    x = 'Species group',
    y = 'Expression (VST)',
    fill = 'Species group'
  ) +
  # Theme
  theme_linedraw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold", color = 'black')
  ) +  # Create facets for the different genes/miRNAs
  facet_grid(cols = vars(genes), scales = "free_y")
crisp2_bar_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/CRISP2_Group/CRISP_and_targeting_miRNAs_2025.04.01.pdf', plot = crisp2_bar_plot, width = 5, height = 3, dpi = 900, create.dir = T)


##### PLA2K ----

# Create a data frame for the gene in quesion
pla2k_df <- sig_df %>% 
  filter(str_detect(genes, 'PLA2K'))

# Get it's miRNAs
mi_pla2k_df <- pla2k_df %>% 
  select(sample.id, genes = miRNA.cluster, pvalue = mi.pvalue, vst = mi.vst, species)

# Format the data frame for joining
pla2k_df <- pla2k_df %>% 
  distinct(sample.id, genes, pvalue = mr.pvalue, vst = mr.vst, species) %>% 
  full_join(mi_pla2k_df)

pla2k_summary_df <- pla2k_df %>% 
  group_by(genes, species) %>%
  summarise(
    mean_vst = mean(vst),
    sd_vst = sd(vst),
    se_vst = sd(vst) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Create a plot combining bar chart and error bars with individual points
pla2k_bar_plot <- ggplot(pla2k_summary_df, aes(x = species, y = mean_vst, fill = species)) +
  # Add bars for the means
  geom_bar(stat = "identity", width = 0.5) +
  # Add error bars for standard error
  geom_errorbar(aes(ymin = mean_vst - se_vst, ymax = mean_vst + se_vst), 
                width = 0.2, position = position_dodge(0.7)) +
  # Add individual data points
  geom_point(data = pla2k_df, aes(x = species, y = vst, group = species),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
             size = 2, alpha = 0.8, shape = 16) +
  # Use your color scheme
  scale_fill_manual(values = species_colors) +
  # Labels
  labs(
    x = 'Species group',
    y = 'Expression (VST)',
    fill = 'Species group'
  ) +
  # Theme
  theme_linedraw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold", color = 'black')
  ) +
  # Create facets for the different genes/miRNAs
  facet_grid(cols = vars(genes), scales = "free_y")
pla2k_bar_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/PLA2K_Group/PLA2K_and_targeting_miRNAs_2025.04.01.pdf', plot = pla2k_bar_plot, width = 4, height = 3, dpi = 900, create.dir = T)


##### SVMP7 ----

# Create a data frame for the gene in quesion
svmp7_df <- sig_df %>% 
  filter(str_detect(genes, 'SVMP7'))

# Get it's miRNAs
mi_svmp7_df <- svmp7_df %>% 
  select(sample.id, genes = miRNA.cluster, pvalue = mi.pvalue, vst = mi.vst, species)

# Format the data frame for joining
svmp7_df <- svmp7_df %>% 
  distinct(sample.id, genes, pvalue = mr.pvalue, vst = mr.vst, species) %>% 
  full_join(mi_svmp7_df)

svmp7_summary_df <- svmp7_df %>% 
  group_by(genes, species) %>%
  summarise(
    mean_vst = mean(vst),
    sd_vst = sd(vst),
    se_vst = sd(vst) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Create a plot combining bar chart and error bars with individual points
svmp7_bar_plot <- ggplot(svmp7_summary_df, aes(x = species, y = mean_vst, fill = species)) +
  # Add bars for the means
  geom_bar(stat = "identity", width = 0.5) +
  # Add error bars for standard error
  geom_errorbar(aes(ymin = mean_vst - se_vst, ymax = mean_vst + se_vst), 
                width = 0.1, position = position_dodge(0.7)) +
  # Add individual data points
  geom_point(data = svmp7_df, aes(x = species, y = vst, group = species),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
             size = 1, alpha = 0.8, shape = 16) +
  # Use your color scheme
  scale_fill_manual(values = species_colors) +
  # Labels
  labs(
    y = 'Expression (VST)',
    fill = 'Species group'
  ) +
  # Theme
  theme_linedraw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold", color = 'black')
  ) +
  # Create facets for the different genes/miRNAs
  facet_grid(cols = vars(genes), scales = "free_y")
svmp7_bar_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVMP7_Group/SVMP7_and_targeting_miRNAs_2025.04.01.pdf', plot = svmp7_bar_plot, width = 10, height = 4, dpi = 900, create.dir = T)



##### SVSP7 ----

# Create a data frame for the gene in quesion
svsp7_df <- sig_df %>% 
  filter(str_detect(genes, 'SVSP7'))

# Get it's miRNAs
mi_svsp7_df <- svsp7_df %>% 
  select(sample.id, genes = miRNA.cluster, pvalue = mi.pvalue, vst = mi.vst, species)

# Format the data frame for joining
svsp7_df <- svsp7_df %>% 
  distinct(sample.id, genes, pvalue = mr.pvalue, vst = mr.vst, species) %>% 
  full_join(mi_svsp7_df)

svsp7_summary_df <- svsp7_df %>% 
  group_by(genes, species) %>%
  summarise(
    mean_vst = mean(vst),
    sd_vst = sd(vst),
    se_vst = sd(vst) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Create a plot combining bar chart and error bars with individual points
svsp7_bar_plot <- ggplot(svsp7_summary_df, aes(x = species, y = mean_vst, fill = species)) +
  # Add bars for the means
  geom_bar(stat = "identity", width = 0.5) +
  # Add error bars for standard error
  geom_errorbar(aes(ymin = mean_vst - se_vst, ymax = mean_vst + se_vst), 
                width = 0.1, position = position_dodge(0.7)) +
  # Add individual data points
  geom_point(data = svsp7_df, aes(x = species, y = vst, group = species),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7),
             size = 1, alpha = 0.8, shape = 16) +
  # Use your color scheme
  scale_fill_manual(values = species_colors) +
  # Labels
  labs(
    y = 'Expression (VST)',
    fill = 'Species group'
  ) +
  # Theme
  theme_linedraw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(face = "bold", color = 'black')
  ) +
  # Create facets for the different genes/miRNAs
  facet_grid(cols = vars(genes), scales = "free_y")
svsp7_bar_plot
ggsave('Figures/DESeq2/P-adjusted/DESeq_Counts/Nicer_Plots/SVSP7_Group/SVSP7_and_targeting_miRNAs_2025.04.01.pdf', plot = svsp7_bar_plot, width = 8, height = 4, dpi = 900, create.dir = T)


