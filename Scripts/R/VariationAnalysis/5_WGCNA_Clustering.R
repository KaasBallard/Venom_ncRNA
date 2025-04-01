# Last Edited: 2025/02/25

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
library(WGCNA)
library(ggdendro)


## Read in the data ----

# Set the path for the data
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_2025.01.22.parquet'

# Read in the data
mi_df <- read_parquet(file = miRNA_mRNA_protein_data) %>% 
  filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM'),
    str_detect(genes, 'Venom_'),
    # I am going to remove CTL6 since it wasn't expressed by any sample
    !str_detect(genes, 'Venom_CTL6')
  )
glimpse(mi_df)

samples <- mi_df %>% distinct(sample.id)
## Create a relationships DataFrame ----
relationships_df <- read_parquet(file = miRNA_mRNA_protein_data) %>% 
  select(
    genes, venom.family, miRNA.cluster, miRNA.cluster.original, feature.type, total.score, total.energy
  ) %>% 
  # Remove badly annotated genes
  dplyr::filter(
    !str_detect(genes, 'maker-scaffold|augustus|XP_|ADAM'),
    str_detect(genes, 'Venom_'),
    # I am going to remove CTL6 since it wasn't expressed by any sample
    !str_detect(genes, 'Venom_CTL6')
  ) %>% 
  distinct()

## Create a color scheme ----

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

# Create color scheme for the venom genes
venom_colors <- c(
  SVMP = SVMP_color,
  ADAM = ADAM_color,
  SVSP = SVSP_color,
  PLA2 = PLA2_color,
  miRNA = miRNA_color,
  VEGF = VEGF_color,
  Ohanin = ohanin_color,
  Myotoxin = myotoxin_color,
  vQC = vQC_color,
  CRISP = CRISP_color,
  CTL = CTL_color,
  EXO = EXO_color,
  LAAO = LAAO_color,
  BPP = BPP_color,
  others = other_color
)


# Format the data for WGCNA ----

# Create a data frame for the genes specifically
genes_df <- mi_df %>% 
  select(
    sample.id, genes, mRNA.vst
  ) %>%
  distinct() %>% 
  pivot_wider(
    names_from = genes,
    values_from = mRNA.vst
  )

# Set the sample ids as the row names
gene_matrix <- genes_df %>% 
  column_to_rownames(var = 'sample.id') %>% 
  as.matrix()


# Create a data frame for the miRNA data
miRNA_df <- mi_df %>% 
  select(
    sample.id, miRNA.cluster, miRNA.vst
  ) %>% 
  distinct() %>% 
  pivot_wider(
    names_from = miRNA.cluster,
    values_from = miRNA.vst
  )

# Set sample IDs as row names for miRNA data
miRNA_matrix <- miRNA_df %>%
  column_to_rownames(var = 'sample.id') %>%
  as.matrix()


# Check to make sure the sample row names are the the same
all(rownames(gene_matrix) == rownames(miRNA_matrix))

# Since the above returned TRUE, I can bind the rows
miRNA_gene_matrix <- cbind(gene_matrix, miRNA_matrix)


# WGCNA ----

## Set up and pick power ----
# Allow multi-threading
allowWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Run the network topology analysis function
# The aim of the function is to help the user pick an appropriate soft-thresholding power for network construction.
sft <- pickSoftThreshold(
  # Input data
  miRNA_gene_matrix,
  powerVector = powers,
  verbose = 5
)
# Output:

# pickSoftThreshold: will use block size 173.
# pickSoftThreshold: calculating connectivity for given powers...
# ..working on genes 1 through 173 of 173
# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1   0.1600  4.160         0.6540   70.00     68.40  91.40
# 2      2   0.2520  4.250         0.2330   40.10     37.60  62.50
# 3      3   0.2360  3.110         0.0294   26.60     24.10  47.40
# 4      4   0.4060  0.953         0.4310   19.10     16.50  38.10
# 5      5   0.2860  0.622         0.8430   14.50     12.10  31.80
# 6      6   0.1220  0.416         0.8060   11.50      9.26  27.30
# 7      7   0.0117  0.112         0.7540    9.31      7.40  23.80
# 8      8   0.1900 -0.429         0.6170    7.74      6.15  21.10
# 9      9   0.3980 -0.608         0.7670    6.54      5.25  18.90
# 10    10   0.5640 -0.738         0.8140    5.62      4.54  17.10
# 11    12   0.7090 -0.897         0.8780    4.28      3.46  14.30
# 12    14   0.7740 -1.120         0.8210    3.39      2.63  12.30
# 13    16   0.8650 -1.170         0.9440    2.75      1.97  10.70
# 14    18   0.8560 -1.240         0.9260    2.28      1.61   9.48
# 15    20   0.8580 -1.200         0.9270    1.93      1.31   8.51

# Create a diagnostic plot to pick a threshold value
par(mfrow = c(1,2));
cex1 = 0.9;

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1, col = "red"
)




# Set the power based on the above plot
soft_power <- 16


## Create network ----

# Run the blockwiseModules function to create the TOM (topological overlap matrix)
network <- blockwiseModules(
  # Input data
  miRNA_gene_matrix,
  
  # Adjacency function
  power = soft_power,
  networkType = 'signed',
  
  # Tree and block options
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  # detectCutHeight = 0.75,
  minModuleSize = 30,
  maxBlockSize = 4000,
  
  # Module adjustments
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  
  # TOM -- Archive the run results in TOM file (saves time)
  saveTOMs = TRUE,
  saveTOMFileBase = "Data/WGCNA/WGCNA_TOM_2025.02.25",
  
  # Output options
  numericLabels = TRUE,
  verbose = 3
)
# Calculating module eigengenes block-wise from all genes
# Flagging genes and samples with too many missing values...
# ..step 1
# ..Working on block 1 .
# TOM calculation: adjacency..
# ..will use 8 parallel threads.
# Fraction of slow calculations: 0.000000
# ..connectivity..
# ..matrix multiplication (system BLAS)..
# ..normalization..
# ..done.
# ..saving TOM for block 1 into file Scripts/R/VariationAnalysis/5_WGCNA_miRNA_and_venom_genes_2025.02.25-block.1.RData
# ....clustering..
# ....detecting modules..
# ....calculating module eigengenes..
# ....checking kME in modules..
# ..removing 4 genes from module 1 because their KME is too low.
# ..removing 10 genes from module 2 because their KME is too low.
# ..removing 3 genes from module 3 because their KME is too low.
# ..merging modules that are too close..
# mergeCloseModules: Merging modules whose distance is less than 0.25
# Calculating new MEs...

# Save the entire network object
save(network, file = "Data/WGCNA/Complete_Network_2025.02.25.RData")

# Convert labels to colors for plotting
merged_colors = labels2colors(network$colors)

# Get gene names
gene_names <- colnames(miRNA_gene_matrix)

# Get the first block's dendrogram and gene names
dendro <- network$dendrograms[[1]]
blockGenes <- network$blockGenes[[1]]
blockGeneNames <- gene_names[blockGenes]

# Make sure the number of labels matches the number of genes in the block
length(blockGeneNames)
length(dendro$order)

# Save the dendrogram
pdf(file = 'Figures/WGCNA/WGCNA_Dendrogram_2025.02.27.pdf', width = 25, height = 15, pointsize = 14)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  network$dendrograms[[1]],
  merged_colors[network$blockGenes[[1]]],
  "Module colors",
  dendroLabels = blockGeneNames,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 
)
dev.off()


### ggplot version of the dendrogram ----

# Extract the dendrogram data
dendro <- network$dendrograms[[1]]
dendro_data <- as.dendrogram(dendro)
dendro_df <- ggdendro::dendro_data(dendro_data)

# Create a data frame for the modules, allowing us to relate treatment to module color
module_df <- data.frame(
  genes = names(network$colors),
  module = labels2colors(network$colors)
)
write_parquet(module_df, sink = 'Data/WGCNA/WGCNA_Loci_Clusters_2025.03.26.parquet')



# Create a smaller version of the relationships data frame so that I can make fontfaces for the module_df
relationships_df2 <- relationships_df %>% 
  distinct(
    genes, miRNA.cluster
  ) %>% 
  rename(miRNA.target.cluster = miRNA.cluster)

# Create a data frame to contain the modules for the miRNAs
miRNA_module_df <- module_df %>% 
  filter(!str_detect(genes, 'Venom_')) %>% 
  rename(miRNA.cluster = genes)

# Create a data frame to contain the modules for the venom genes
gene_module_df <- module_df %>% 
  filter(str_detect(genes, 'Venom_'))

# Fuse the modules together
gene_miRNA_module_df <- full_join(
  miRNA_module_df,
  gene_module_df,
  by = 'module',
  relationship = 'many-to-many'
) %>% 
  full_join(
    relationships_df2,
    by = 'genes',
    relationship = 'many-to-many'
  ) %>% 
  mutate(
    # Set the font to bold if the miRNA or gene clusters together based on WGCNA and share targeting based on miRanda
    font = ifelse(miRNA.target.cluster == miRNA.cluster, 'bold', 'plain')
  )

# Resolve conflicting font values for genes/miRNAs by assigning 'bold' if any 'bold' font exists for a gene/miRNA
gene_miRNA_module_df <- gene_miRNA_module_df %>%
  group_by(genes, module) %>%
  mutate(font = ifelse(any(font == 'bold'), 'bold', 'plain')) %>%
  ungroup() %>% 
  group_by(miRNA.cluster, module) %>%
  mutate(font = ifelse(any(font == 'bold'), 'bold', 'plain')) %>%
  ungroup()

# Split the data again
# Remap the miRNA moduled data frame so that it contains the font data
miRNA_module_df2 <- gene_miRNA_module_df %>% 
  select(
    miRNA.cluster, module, font
  ) %>% 
  distinct() %>% 
  rename(genes = miRNA.cluster)

# Remap the gene moduled data frame so that it contains the font data
gene_module_df2 <- gene_miRNA_module_df %>% 
  select(
    genes, module, font
  ) %>% 
  distinct()

# Fuse back to the module data frame
module_df2 <- module_df %>% 
  full_join(
    gene_module_df2,
    by = c('genes', 'module')
  ) %>% 
  full_join(
    miRNA_module_df2,
    by = c('genes', 'module', 'font')
  ) %>% 
  filter(!is.na(font))

# Extract the correct gene order from the dendrogram structure
dendro_gene_order <- dendro$order

# Get the gene names in the order they appear in the dendrogram
ordered_gene_names <- gene_names[dendro_gene_order]

# Reorder module_df to match the new labels
ordered_modules <- module_df2 %>%
  filter(genes %in% ordered_gene_names) %>%
  arrange(match(genes, ordered_gene_names)) 

# Replace the numerical labels with gene names
dendro_df$labels$label <- ordered_gene_names

# Fuse the ordered_modules data to the dendro_df
dendro_df$labels <- dendro_df$labels %>% 
  left_join(
    ordered_modules,,
    by = c('label' = 'genes')
  )

# Set a colors vector
colors <- setNames(unique(ordered_modules$module), unique(dendro_df$labels$module))

# Create plot
cluster_dendrogram <- ggplot() +
  geom_segment(
    data = dendro_df$segments,
    aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  geom_text(
    data = dendro_df$labels,
    aes(x = x, y = y, label = label, color = module, fontface = font),
    hjust = 1,
    size = 2,
    angle = 90
  ) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = 'bold'),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "lines")  # Increase right margin
  ) +
  scale_y_continuous(expand = c(0.2, 0)) +
  labs(
    title = "Clustering of venom genes and miRNAs",
    y = "Height",
    x = "Genes and miRNAs",
    color = "Modules",
  ) 
cluster_dendrogram
ggsave(filename = 'Figures/WGCNA/WGCNA_Dendrogram_ggplot_version_2025.02.27.pdf', plot = cluster_dendrogram, create.dir = TRUE, height = 10, width = 15, dpi = 900)


### Check if the the miRNAs are in modules with at least 1 venom gene they target in the 3UTR and 5UTR ----

# Create a smaller version of the relationships data frame so that I can make fontfaces for the module_df
relationships_3utr_and_5utr_df <- relationships_df %>% 
  filter(
    feature.type == c('three_prime_utr', 'five_prime_utr')
  ) %>% 
  distinct(
    genes, miRNA.cluster
  ) %>% 
  rename(miRNA.target.cluster = miRNA.cluster)

# Fuse the modules together
gene_miRNA_3_5utr_module_df <- full_join(
  miRNA_module_df,
  gene_module_df,
  by = 'module',
  relationship = 'many-to-many'
) %>% 
  # Right join to filter out any miRNAs that don't exist anymore because they targeted something only in the CDS region
  right_join(
    relationships_3utr_and_5utr_df,
    by = 'genes',
    relationship = 'many-to-many'
  ) %>% 
  mutate(
    # Set the font to bold if the miRNA or gene clusters together based on WGCNA and share targeting based on miRanda
    same.module = ifelse(miRNA.target.cluster == miRNA.cluster, 'yes', 'no')
  )

# Resolve conflicting font values for genes/miRNAs by assigning 'bold' if any 'bold' font exists for a gene/miRNA
gene_miRNA_3_5utr_module_df <- gene_miRNA_3_5utr_module_df %>%
  group_by(genes, module) %>%
  mutate(same.module = ifelse(any(same.module == 'yes'), 'yes', 'no')) %>%
  ungroup() %>%
  group_by(miRNA.cluster, module) %>%
  mutate(same.module = ifelse(any(same.module == 'yes'), 'yes', 'no')) %>%
  ungroup()
# YEAH, all of them still target at least one gene they are in the same module as a venom gene they target!!!!!!


### Check if the the miRNAs are in modules with at least 1 venom gene they target in the 3UTR only ----

# Create a smaller version of the relationships data frame so that I can make fontfaces for the module_df
relationships_3utr_df <- relationships_df %>% 
  filter(
    feature.type == c('three_prime_utr')
  ) %>% 
  distinct(
    genes, miRNA.cluster
  ) %>% 
  rename(miRNA.target.cluster = miRNA.cluster)

# Fuse the modules together
gene_miRNA_3utr_module_df <- full_join(
  miRNA_module_df,
  gene_module_df,
  by = 'module',
  relationship = 'many-to-many'
) %>% 
  # Right join to filter out any miRNAs that don't exist anymore because they targeted something only in the CDS region
  right_join(
    relationships_3utr_df,
    by = 'genes',
    relationship = 'many-to-many'
  ) %>% 
  mutate(
    # Set the font to bold if the miRNA or gene clusters together based on WGCNA and share targeting based on miRanda
    same.module = ifelse(miRNA.target.cluster == miRNA.cluster, 'yes', 'no')
  ) %>% 
  select(
    -miRNA.cluster
  )

# Resolve conflicting font values for genes/miRNAs by assigning 'bold' if any 'bold' font exists for a gene/miRNA
gene_miRNA_3utr_module_df <- gene_miRNA_3utr_module_df %>%
  group_by(genes, module) %>%
  mutate(same.module = ifelse(any(same.module == 'yes'), 'yes', 'no')) %>%
  ungroup() %>%
  group_by(miRNA.target.cluster, module) %>%
  mutate(same.module = ifelse(any(same.module == 'yes'), 'yes', 'no')) %>%
  ungroup()
# YEAH, all of them still target at least one gene they are in the same module as a venom gene they target!!!!!!





## Calculate Eigengenes ----

# WGCNA can calculate the hypothetical central gene for each module
# Get module Eigengenes per cluster
module_eigenegenes <- moduleEigengenes(miRNA_gene_matrix, merged_colors)$eigengenes

# Reorder modules so similar modules are next to each other
module_eigenegenes <- orderMEs(module_eigenegenes)

# Get module order
module_order <- names(module_eigenegenes) %>% gsub("ME", "", .)
print(module_order)
# [1] "turquoise" "blue"      "brown"     "grey"   
# That is not many modules

# Add treatment names
module_eigenegenes$treatment <- row.names(module_eigenegenes)

# Create data frame from the module information
eigengenes_df <- module_eigenegenes %>% 
  pivot_longer(-treatment) %>% 
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

# Create a plot of module-trait relationships
module_trait_plot <- ggplot(eigengenes_df, aes(x = treatment, y = name, fill = value)) +
  geom_tile() +
  theme_linedraw() +
  scale_fill_gradient2(
    low = 'blue',
    mid = 'white',
    high = 'red',
    midpoint = 0,
    limit = c(-1, 1)
  ) +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  labs(
    title = 'Module-trait Relationships',
    y = 'Modules',
    fill = 'Correlation'
  )
module_trait_plot


## Examine module expression ----
# I will just use the module order

# Get miRNA data
miRNA_df <- mi_df %>% 
  select(
    sample.id, miRNA.cluster, miRNA.vst
  ) %>% 
  distinct() %>% 
  rename(genes = miRNA.cluster, expression = miRNA.vst)

# Get gene data
gene_df <- mi_df %>% 
  select(
    sample.id, genes, mRNA.vst
  ) %>% 
  distinct() %>% 
  rename(expression = mRNA.vst)

# Fuse the data frames to create a expression data frame with module information
module_expression_df <- full_join(
  miRNA_df,
  gene_df,
  by = names(gene_df)
) %>% 
  full_join(
    module_df,
    by = 'genes'
  ) 

# Create a plot of the module
module_expression_plot <- ggplot(module_expression_df, aes(x = sample.id, y = expression, group = genes)) +
  geom_line(aes(color = module), alpha = 0.2) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(
    rows = vars(module)
  ) +
  labs(
    x = "Sample",
    y = 'Normalized expression'
  )
module_expression_plot


## TOM plot ----

# Extract the gene dendrogram as geneTree
geneTree <- network$dendrograms[[1]]

# Extract the module colors
module_colors <- labels2colors(network$colors)

# Calculate adjacency and TOM
adjacency_matrix <- adjacency(miRNA_gene_matrix, power = soft_power)
TOM <- TOMsimilarity(adjacency_matrix)

# Calculate the dissimilarity from the loaded TOM
dissTOM <- 1 - TOM

# Create the TOMplot
TOMplot(
  dissTOM,
  geneTree,
  module_colors,
  main = 'TOM Plot'
)


# Custom TOMplot function with dendrogram labels
customTOMplot <- function(TOM, dendro, Colors, labels = NULL, main = "Network heatmap plot", 
                          labelSize = 0.5, labelRotation = 90, ...) {
  # Convert the TOM to a distance
  dissTOM = 1 - TOM
  
  # Set up the plotting area 
  par(mar = c(12, 5, 5, 5))
  
  # Call the standard TOMplot function
  TOMplot(dissTOM, dendro, Colors, main = main, ...)
  
  # If labels are provided, add them to the plot
  if (!is.null(labels)) {
    # Get the order from the dendrogram
    order <- dendro$order
    
    # Add the labels
    axis(1, at = 1:length(labels[order]), labels = labels[order], 
         las = 2, cex.axis = labelSize, srt = labelRotation)
  }
}

# Save the below plot
pdf(file = 'Figures/WGCNA/TOMplot_2025.02.27.pdf', width = 20, height = 20)

# Example usage:
customTOMplot(
  TOM,
  geneTree,
  module_colors,
  labels = gene_names,
  main = 'TOM Plot'
)
dev.off()




## Create a cytoscape data set ----

# First, set a threshold for the TOM values to determine which connections to include
TOM_threshold <- 0.05  # Adjust this based on your network's connectivity

# Get the module assignments for each gene and miRNA
node_df <- data.frame(
  NodeID = colnames(miRNA_gene_matrix),
  Module = labels2colors(network$colors),
  stringsAsFactors = FALSE
)

# Add information about whether each node is an miRNA or a gene
node_df <- node_df %>%
  mutate(
    NodeType = ifelse(str_detect(NodeID, "Venom_"), "gene", "miRNA")
  )

# Add venom family information for genes
node_df <- node_df %>%
  left_join(
    relationships_df %>% select(genes, venom.family) %>% distinct(),
    by = c("NodeID" = "genes")
  )

# Create edge data from the TOM matrix
# Convert TOM to a data frame with node pairs and weights
edge_data <- matrix(0, nrow = sum(TOM > TOM_threshold), ncol = 3)
colnames(edge_data) <- c("from", "to", "weight")
counter <- 1

# Loop through the TOM matrix to extract edges above the threshold
for (i in 1:(nrow(TOM) - 1)) {
  for (j in (i + 1):ncol(TOM)) {
    if (TOM[i, j] > TOM_threshold) {
      edge_data[counter, ] <- c(i, j, TOM[i, j])
      counter <- counter + 1
    }
  }
}

# Convert to a data frame
edge_data <- as.data.frame(edge_data)

# Replace indices with actual gene/miRNA names
edge_data$from <- colnames(TOM)[edge_data$from]
edge_data$to <- colnames(TOM)[edge_data$to]

# Save node and edge data for Cytoscape
write.table(node_df, file = "Data/Cytoscape/WGCNA/WGCNA_NodeData_2025.02.27.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(edge_data, file = "Data/Cytoscape/WGCNA/WGCNA_EdgeData_2025.02.27.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Alternatively, you can use the WGCNA built-in function to export to Cytoscape
# This creates a more comprehensive export with additional network statistics

# First, retrieve module-specific nodes
modules <- unique(node_df$Module)
modules <- modules[modules != "grey"]  # Typically, grey is reserved for unassigned genes

# Create a list to store module data
cytoscape_output <- list()

for (module in modules) {
  # Get genes in this module
  module_genes <- node_df$NodeID[node_df$Module == module]
  
  # Extract the subnetwork for this module
  module_TOM <- TOM[module_genes, module_genes]
  
  # Set up the Cytoscape export using WGCNA's function
  cytoscape_output[[module]] <- exportNetworkToCytoscape(
    adjMat = module_TOM,
    edgeFile = paste0("Data/Cytoscape/CytoscapeEdges-", module, "_2025.02.27.txt"),
    nodeFile = paste0("Data/Cytoscape/CytoscapeNodes-", module, "_2025.02.27.txt"),
    weighted = TRUE,
    threshold = TOM_threshold,
    nodeNames = module_genes,
    nodeAttr = node_df %>% 
      filter(NodeID %in% module_genes) %>% 
      select(Module, NodeType, venom.family) %>% 
      as.data.frame()
  )
}

# You can also export the whole network if desired
whole_network <- exportNetworkToCytoscape(
  adjMat = TOM,
  edgeFile = "Data/Cytoscape/CytoscapeEdges-WholeNetwork_2025.02.27.txt",
  nodeFile = "Data/Cytoscape/CytoscapeNodes-WholeNetwork_2025.02.27.txt",
  weighted = TRUE,
  threshold = TOM_threshold,
  nodeNames = colnames(TOM),
  nodeAttr = node_df %>% select(Module, NodeType, venom.family) %>% as.data.frame()
)

# Save the complete Cytoscape output for reference
saveRDS(cytoscape_output, file = "Data/Cytoscape/WGCNA_CytoscapeOutput_2025.02.27.rds")


