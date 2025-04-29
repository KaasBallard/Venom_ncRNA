# Last Edited: 2025/02/06

# Set up and Read Data ----

## Load in packages ----
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
library(scales)

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

## Prepare data ----
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

# Create a data frame for the venom interactions only
venom_inter_df <- interactions_df %>% dplyr::filter(str_detect(genes, 'Venom_')) %>% distinct()

# Join results into the same data frame again
merged_diff_df <- interactions_df %>% 
  left_join(res_mirna, by = 'miRNA.cluster') %>% 
  left_join(res_mrna, by = 'genes') %>% 
  dplyr::rename(
    baseMean.mi = baseMean.x, log2FoldChange.mi = log2FoldChange.x, lfcSE.mi = lfcSE.x, stat.mi = stat.x, pvalue.mi = pvalue.x, padj.mi = padj.x,
    baseMean.mr = baseMean.y, log2FoldChange.mr = log2FoldChange.y, lfcSE.mr = lfcSE.y, stat.mr = stat.y, pvalue.mr = pvalue.y, padj.mr = padj.y
  ) %>% 
  filter(!str_detect(genes, 'maker-scaffold|augustus|XP_')) # This should get rid of all of the weirdly annotated genes, make sure to note you did this in the methods section


# Filter down the above using the p-values rather than the adjusted p-values
significant_diff_paired_df <- merged_diff_df %>% 
  filter(pvalue.mi < 0.05 & pvalue.mr < 0.05) %>% 
  filter(!is.na(pvalue.mi) & !is.na(pvalue.mr))

# Filter for each individually as well
significant_mi_diff_df <- merged_diff_df %>% filter(pvalue.mi < 0.05) %>% filter(!is.na(pvalue.mi))
significant_mr_diff_df <- merged_diff_df %>% filter(pvalue.mr < 0.05) %>% filter(!is.na(pvalue.mr))

## DESeq count plots for significant venom genes and targeting miRNA pairs ----

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
  plot_gene_counts(dds, gene = gene, intgroup = 'species', saveDir = 'Figures/DESeq2/P-unadjusted/DESeq_Counts/Species/Venom/genes/Significant_Paired', date = '2025.02.07')
}

# Create a list of miRNAs
sig_miRNAs <- venom_sig_mi_paired_df$miRNA.cluster

# Run the function in a loop for the miRNAs as well
for (miRNA in sig_miRNAs) {
  plot_gene_counts(dds, gene = miRNA, intgroup = 'species', saveDir = 'Figures/DESeq2/P-unadjusted/DESeq_Counts/Species/Venom/miRNAs/Significant_Paired', date = '2025.02.07')
}

