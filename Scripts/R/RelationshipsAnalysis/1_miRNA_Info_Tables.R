# Last Edited: 2025/02/04

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
library(arrow)
# Source my functions
source('/Users/ballardk/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')
# source('/Users/kaasballard/Library/CloudStorage/OneDrive-UTArlington/Bin/R/MyFunctions/MyFunctions.R')

# Create variable for the fused data set.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read in the data
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)

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


# Check # of miRNAs ----

# Check all miRNAs 
miRNA_only_df <- mi_df %>% 
  dplyr::select(miRNA.cluster) %>% distinct() %>% rowwise() %>% mutate(index = row_number()) %>% ungroup()

# Check venom targeting miRNAs
venom_miRNA_only_df <- mi_df %>% 
  filter(str_starts(genes, 'Venom'), !str_detect(genes, 'ADAM')) %>% 
  dplyr::select(miRNA.cluster) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(index = row_number()) %>% 
  ungroup()


# Tables by feature type ----

## Gene table for 3UTR ----

# Create table of genes regulated by miRNAs
gene_miRNA_table_3utr <- mi_df %>%
  filter(feature.type == 'three_prime_utr') %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )

# Now save that table as a csv
write_csv3(gene_miRNA_table_3utr, file = 'Tables/Genes_and_Putative_miRNAs/Genes_and_putative_targeting_miRNAs_in_3UTR_2025.02.04.csv')

# Now create the same thing but only for venom genes.
# Create table of genes regulated by miRNAs
venom_gene_miRNA_table_3utr <- mi_df %>%
  filter(feature.type == 'three_prime_utr') %>%
  filter(str_starts(genes, 'Venom_'), !str_detect(genes, 'ADAM')) %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
View(venom_gene_miRNA_table_3utr)
# Now save that table as a csv
write_csv3(venom_gene_miRNA_table_3utr, file = 'Tables/Genes_and_Putative_miRNAs/Venom_genes_and_putative_targeting_miRNAs_in_3UTR_2025.02.04.csv')

# Now create the same thing but only for non-venom genes.
# Create table of genes regulated by miRNAs
non_venom_gene_miRNA_table_3utr <- mi_df %>%
  filter(feature.type == 'three_prime_utr') %>%
  filter(!str_starts(genes, 'Venom')) %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
# Now save that table as a csv
write_csv3(non_venom_gene_miRNA_table_3utr, file = 'Tables/Genes_and_Putative_miRNAs/Non-Venom_genes_and_putative_targeting_miRNAs_in_3UTR_2025.02.04.csv')




## Gene table for 5UTR ----

# Create table of genes regulated by miRNAs
gene_miRNA_table_5utr <- mi_df %>%
  filter(feature.type == 'five_prime_utr') %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )

# Now save that table as a csv
write_csv3(gene_miRNA_table_5utr, file = 'Tables/Genes_and_Putative_miRNAs/Genes_and_putative_targeting_miRNAs_in_5UTR_2025.02.04.csv')

# Now create the same thing but only for venom genes.
# Create table of genes regulated by miRNAs
venom_gene_miRNA_table_5utr <- mi_df %>%
  filter(feature.type == 'five_prime_utr') %>%
  filter(str_starts(genes, 'Venom_|PLA2G2E.1')) %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
View(venom_gene_miRNA_table_5utr)
# Now save that table as a csv
write_csv3(venom_gene_miRNA_table_5utr, file = 'Tables/Genes_and_Putative_miRNAs/Venom_genes_and_putative_targeting_miRNAs_in_5UTR_2025.02.04.csv')

# Now create the same thing but only for non-venom genes.
# Create table of genes regulated by miRNAs
non_venom_gene_miRNA_table_5utr <- mi_df %>%
  filter(feature.type == 'five_prime_utr') %>%
  filter(!str_starts(genes, 'Venom')) %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
View(non_venom_gene_miRNA_table_5utr)
# Now save that table as a csv
write_csv3(non_venom_gene_miRNA_table_5utr, file = 'Tables/Genes_and_Putative_miRNAs/Non-Venom_genes_and_putative_targeting_miRNAs_in_5UTR_2025.02.04.csv')




## Gene table for CDS ----

# Create table of genes regulated by miRNAs
gene_miRNA_table_cds <- mi_df %>%
  filter(feature.type == 'CDS') %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )

# Now save that table as a csv
write_csv3(gene_miRNA_table_cds, file = 'Tables/Genes_and_Putative_miRNAs/Genes_and_putative_targeting_miRNAs_in_CDS_2025.02.04.csv')

# Now create the same thing but only for venom genes.
# Create table of genes regulated by miRNAs
venom_gene_miRNA_table_cds <- mi_df %>%
  filter(feature.type == 'CDS') %>%
  filter(str_starts(genes, 'Venom_|PLA2G2E.1')) %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
# Now save that table as a csv
write_csv3(venom_gene_miRNA_table_cds, file = 'Tables/Genes_and_Putative_miRNAs/Venom_genes_and_putative_targeting_miRNAs_in_CDS_2025.02.04.csv')

# Now create the same thing but only for non-venom genes.
# Create table of genes regulated by miRNAs
non_venom_gene_miRNA_table_cds <- mi_df %>%
  filter(feature.type == 'CDS') %>%
  filter(!str_starts(genes, 'Venom')) %>% 
  distinct(genes, miRNA.cluster.original, miRNA.cluster) %>%
  group_by(genes) %>%
  summarize(
    miRNA.cluster = paste(miRNA.cluster, collapse = ", "),
    Number.of.miRNAs.per.Gene = n()
  ) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Number of miRNAs per Gene' = 'Number.of.miRNAs.per.Gene'
  )
# Now save that table as a csv
write_csv3(non_venom_gene_miRNA_table_cds, file = 'Tables/Genes_and_Putative_miRNAs/Non-Venom_genes_and_putative_targeting_miRNAs_in_CDS_2025.02.04.csv')


# miRNAs and putatively targeted genes -----

## miRNAs and putatively targeted genes in the 3UTR ----

# Now I want to create the opposite thing, clusters with their genes
miRNA_gene_table_3utr <- mi_df %>% 
  filter(feature.type == 'three_prime_utr') %>%
  distinct(miRNA.cluster, miRNA.sequence, hairpin.sequence, miRNA.cluster.original, genes) %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    genes = paste(genes, collapse = ", "),
    miRNA.sequence = paste(unique(miRNA.sequence), collapse = ", "),
    hairpin.sequence = paste(unique(hairpin.sequence), collapse = ", "),
    Number.of.genes.per.miRNA = n()
  ) %>% 
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Mature miRNA Sequence' = 'miRNA.sequence',
    'miRNA Hairpin Sequence' = 'hairpin.sequence',
    'Target genes' = 'genes',
    'Number of genes per miRNA' = 'Number.of.genes.per.miRNA'
  ) %>% 
  dplyr::select(
    'miRNA Locus', 'Mature miRNA Sequence', 'miRNA Hairpin Sequence', 'Target genes', 'Number of genes per miRNA'
  )
# Now save that table as a csv
write_csv3(miRNA_gene_table_3utr, file = 'Tables/miRNAs_and_Targets/miRNAs_and_targeted_genes_in_3UTR_2025.02.04.csv')

# Create a table for venom targeting miRNAs and the targeted genes
venom_miRNA_gene_table_3utr <- mi_df %>% 
  filter(feature.type == 'three_prime_utr') %>%
  filter(str_starts(genes, 'Venom_'), !str_detect(genes, 'ADAM')) %>% 
  distinct(miRNA.cluster, miRNA.sequence, hairpin.sequence, miRNA.cluster.original, genes) %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    genes = paste(genes, collapse = ", "),
    miRNA.sequence = paste(unique(miRNA.sequence), collapse = ", "),
    hairpin.sequence = paste(unique(hairpin.sequence), collapse = ", "),
    Number.of.genes.per.miRNA = n()
  ) %>% 
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Mature miRNA Sequence' = 'miRNA.sequence',
    'miRNA Hairpin Sequence' = 'hairpin.sequence',
    'Target genes' = 'genes',
    'Number of genes per miRNA' = 'Number.of.genes.per.miRNA'
  ) %>% 
  dplyr::select(
    'miRNA Locus', 'Mature miRNA Sequence', 'miRNA Hairpin Sequence', 'Target genes', 'Number of genes per miRNA'
  )
# Now save that table as a csv
write_csv3(venom_miRNA_gene_table_3utr, file = 'Tables/miRNAs_and_Targets/miRNAs_and_targeted_Venom_genes_in_3UTR_2025.02.04.csv')


## miRNAs and putatively targeted genes in the 5UTR ----

# Now I want to create the opposite thing, clusters with their genes
miRNA_gene_table_5utr <- mi_df %>% 
  filter(feature.type == 'five_prime_utr') %>%
  distinct(miRNA.cluster, miRNA.sequence, hairpin.sequence, miRNA.cluster.original, genes) %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    genes = paste(genes, collapse = ", "),
    miRNA.sequence = paste(unique(miRNA.sequence), collapse = ", "),
    hairpin.sequence = paste(unique(hairpin.sequence), collapse = ", "),
    Number.of.genes.per.miRNA = n()
  ) %>% 
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Mature miRNA Sequence' = 'miRNA.sequence',
    'miRNA Hairpin Sequence' = 'hairpin.sequence',
    'Target genes' = 'genes',
    'Number of genes per miRNA' = 'Number.of.genes.per.miRNA'
  ) %>% 
  dplyr::select(
    'miRNA Locus', 'Mature miRNA Sequence', 'miRNA Hairpin Sequence', 'Target genes', 'Number of genes per miRNA'
  )
# Now save that table as a csv
write_csv3(miRNA_gene_table_5utr, file = 'Tables/miRNAs_and_Targets/miRNAs_and_targeted_genes_in_5UTR_2025.02.04.csv')

# Create a table for venom targeting miRNAs and the targeted genes
venom_miRNA_gene_table_5utr <- mi_df %>% 
  filter(feature.type == 'five_prime_utr') %>%
  filter(str_starts(genes, 'Venom_'), !str_detect(genes, 'ADAM')) %>% 
  distinct(miRNA.cluster, miRNA.sequence, hairpin.sequence, miRNA.cluster.original, genes) %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    genes = paste(genes, collapse = ", "),
    miRNA.sequence = paste(unique(miRNA.sequence), collapse = ", "),
    hairpin.sequence = paste(unique(hairpin.sequence), collapse = ", "),
    Number.of.genes.per.miRNA = n()
  ) %>% 
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Mature miRNA Sequence' = 'miRNA.sequence',
    'miRNA Hairpin Sequence' = 'hairpin.sequence',
    'Target genes' = 'genes',
    'Number of genes per miRNA' = 'Number.of.genes.per.miRNA'
  ) %>% 
  dplyr::select(
    'miRNA Locus', 'Mature miRNA Sequence', 'miRNA Hairpin Sequence', 'Target genes', 'Number of genes per miRNA'
  )
# Now save that table as a csv
write_csv3(venom_miRNA_gene_table_5utr, file = 'Tables/miRNAs_and_Targets/miRNAs_and_targeted_Venom_genes_in_5UTR_2025.02.04.csv')


## miRNAs and putatively targeted genes in the CDS ----

# Now I want to create the opposite thing, clusters with their genes
miRNA_gene_table_cds <- mi_df %>% 
  filter(feature.type == 'CDS') %>%
  distinct(miRNA.cluster, miRNA.sequence, hairpin.sequence, miRNA.cluster.original, genes) %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    genes = paste(genes, collapse = ", "),
    miRNA.sequence = paste(unique(miRNA.sequence), collapse = ", "),
    hairpin.sequence = paste(unique(hairpin.sequence), collapse = ", "),
    Number.of.genes.per.miRNA = n()
  ) %>% 
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Mature miRNA Sequence' = 'miRNA.sequence',
    'miRNA Hairpin Sequence' = 'hairpin.sequence',
    'Target genes' = 'genes',
    'Number of genes per miRNA' = 'Number.of.genes.per.miRNA'
  ) %>% 
  dplyr::select(
    'miRNA Locus', 'Mature miRNA Sequence', 'miRNA Hairpin Sequence', 'Target genes', 'Number of genes per miRNA'
  )
# Now save that table as a csv
write_csv3(miRNA_gene_table_cds, file = 'Tables/miRNAs_and_Targets/miRNAs_and_targeted_genes_in_CDS_2025.02.04.csv')

# Create a table for venom targeting miRNAs and the targeted genes
venom_miRNA_gene_table_cds <- mi_df %>% 
  filter(feature.type == 'CDS') %>%
  filter(str_starts(genes, 'Venom_'), !str_detect(genes, 'ADAM')) %>% 
  distinct(miRNA.cluster, miRNA.sequence, hairpin.sequence, miRNA.cluster.original, genes) %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    genes = paste(genes, collapse = ", "),
    miRNA.sequence = paste(unique(miRNA.sequence), collapse = ", "),
    hairpin.sequence = paste(unique(hairpin.sequence), collapse = ", "),
    Number.of.genes.per.miRNA = n()
  ) %>% 
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster',
    'Mature miRNA Sequence' = 'miRNA.sequence',
    'miRNA Hairpin Sequence' = 'hairpin.sequence',
    'Target genes' = 'genes',
    'Number of genes per miRNA' = 'Number.of.genes.per.miRNA'
  ) %>% 
  dplyr::select(
    'miRNA Locus', 'Mature miRNA Sequence', 'miRNA Hairpin Sequence', 'Target genes', 'Number of genes per miRNA'
  )
# Now save that table as a csv
write_csv3(venom_miRNA_gene_table_cds, file = 'Tables/miRNAs_and_Targets/miRNAs_and_targeted_Venom_genes_in_CDS_2025.02.04.csv')



# Putative names based on BLAST ----

# Create a table for the best blast hits for each miRNA
putative_miRNA_names <- mi_df %>% 
  distinct(miRNA.cluster.original, best.miRNA.ortholog, miRNA.cluster) %>%
  group_by(miRNA.cluster.original) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster.original',
    'Best Blast Hit' = 'best.miRNA.ortholog',
    'Putative miRNA Name' = 'miRNA.cluster'
  )
write_csv3(putative_miRNA_names, file = 'Tables/Putative_miRNA_Names/Putative_miRNA_Names_2025.02.04.csv')

# Create a table for the best blast hits for each miRNA that targets a venom gene
venom_putative_miRNA_names <- mi_df %>% 
  filter(str_starts(genes, 'Venom_'), !str_detect(genes, 'ADAM')) %>% 
  distinct(miRNA.cluster.original, best.miRNA.ortholog, miRNA.cluster) %>%
  group_by(miRNA.cluster.original) %>%
  dplyr::rename(
    'miRNA Locus' = 'miRNA.cluster.original',
    'Best Blast Hit' = 'best.miRNA.ortholog',
    'Putative miRNA Name' = 'miRNA.cluster'
  )
write_csv3(venom_putative_miRNA_names, file = 'Tables/Putative_miRNA_Names/Venom_targeting_putative_miRNA_Names_2025.02.04.csv')

