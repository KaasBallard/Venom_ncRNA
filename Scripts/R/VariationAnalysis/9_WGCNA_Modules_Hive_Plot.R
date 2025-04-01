# Last Edited: 2025/03/27

# Set up and Read Data ----

## Load in packages ----
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
library(ggraph)
library(tidygraph)
library(igraph)
library(tidygraph)
library(ggtree)
library(scales)
library(arrow)
library(ggforce)

## Read data ----
# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_Venom_2025.01.22.parquet'
# Create a variable for the mRNA and protein data
mRNA_protein_data <- 'Data/Merged/mRNA_Protein_Combined_Data_2025.01.22.parquet'

# Create a df name and do some minor tweaks to it's structure for readability
mi_df <- read_parquet(file = miRNA_mRNA_protein_data) %>% 
  dplyr::filter(
    !sample.id == 'CV1082_viridis',
    str_detect(genes, 'Venom_'),
    !str_detect(genes, 'ADAM|Venom_CTL6'), # CTL6 wasn't expressed and had to be removed from WGCNA
    # feature.type == 'three_prime_utr'
  ) %>% 
  dplyr::select(
    genes, venom.family, miRNA.cluster, total.energy, total.score, feature.type
  ) %>% 
  dplyr::distinct()
glimpse(mi_df)

# Edge colors for miRNA targeting
edge_colors <- c(
  Same = '#FF5F1F',
  Different = '#ADD8E6'
)

## WGCNA Module Data ---
# Set the path to the WGCNA modules
wgcna_modules <- 'Data/WGCNA/WGCNA_Loci_Clusters_2025.03.26.parquet'

# Read in the WGNCA modules
wgcna_df <- read_parquet(file = wgcna_modules)

# Get miRNA modules
miRNA_wgcna_df <- wgcna_df %>% 
  dplyr::filter(str_detect(genes, 'cvi-|Cluster_')) %>% 
  dplyr::rename(miRNA.cluster = genes, miRNA.module = module)

# Get gene modules
gene_wgcna_df <- wgcna_df %>%
  dplyr::filter(str_detect(genes, 'Venom')) %>% 
  dplyr::rename(gene.module = module)

# WGCNA module colors
module_colors <- c(
  blue_Gene = '#0000FF',
  blue_miRNA = '#1A5EFF',
  brown_Gene = '#A52A2A',
  brown_miRNA = '#B84545',
  grey_Gene = '#C0C0C0',
  grey_miRNA = '#D3D3D3',
  turquoise_Gene = '#40E0D0',
  turquoise_miRNA = '#77F2E0'
)

# Format data for hive plot ----

# Create an interaction table
# Note that I am using left join so that only genes that were targeted by something in the 3UTR are kept
interaction_table <- dplyr::left_join( 
  mi_df,
  gene_wgcna_df,
  by = c('genes')
) %>%
  # Join the miRNA module data, again doing a left_join so that only miRNAs that target something in the 
  dplyr::left_join(
    miRNA_wgcna_df, by = c('miRNA.cluster')
  ) %>% 
  distinct() %>% 
  # Categorize miRNAs and genes as in the same module or not
  mutate(
    same.module = if_else(gene.module == miRNA.module, "Same", "Different")
  )

# Create nodes
nodes <- bind_rows(
  interaction_table %>% 
    select(name = miRNA.cluster, color = miRNA.module) %>% 
    mutate(type = 'miRNA'),
  interaction_table %>% 
    select(name = genes, color = gene.module) %>% 
    mutate(type = 'Gene')
) %>% 
  distinct() %>% 
  mutate(
    module = str_c(color, type, sep = "_")
  )
  
# Create edges
edges <- interaction_table %>% 
  select(from = miRNA.cluster, to = genes, same.module, total.energy) %>% 
  distinct()

# Network data
network_graph <- tbl_graph(
  nodes = nodes,
  edges = edges,
  directed = FALSE
)
glimpse(network_graph)

## Hive Plot ----
# Create hive plot
wgcna_hive_plot <- ggraph::ggraph(network_graph, 'hive', axis = module, sort.by = 'degreee') +
  # Add edges first
  geom_edge_hive(
    aes(
      colour = same.module,
      alpha = total.energy
    ),
    show.legend = TRUE
  ) +
  # Add nodes with module-based coloring
  geom_node_point(
    aes(
      colour = module
    ),
    alpha = 0.8,
    size = 2
  ) +
  # # 81 labels can't be drawn, so I will comment it out
  # geom_node_text(
  #   aes(label = name),
  #   size = 2,
  #   repel = TRUE,
  #   max.overlaps = 10
  # ) +
  scale_edge_color_manual(
    values = edge_colors,
    name = 'Module interaction'
  ) +
  scale_edge_alpha_continuous(
    name = 'Binding Energy (kcal/mol)',
    range = c(0.05, 1),
    trans = 'reverse'
  ) +
  scale_color_manual(
    values = module_colors,
    name = 'Node Module'
  ) +
  coord_fixed() +
  theme_graph() +
  theme(
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5)
  )
  # theme_no_axes(base.theme = theme_linedraw())
wgcna_hive_plot

# With facet
wgcna_hive_plot_facet <- wgcna_hive_plot + facet_edges(~same.module)
wgcna_hive_plot_facet


## Quantifying Shared vs Different relationships ----

# Count the number of existing miRNA gene relationships (pairs)
miRNA_gene_pairs_df <- interaction_table %>% 
  group_by(miRNA.module, gene.module, same.module) %>% 
  summarise(
    pairs = n()
  )

# Calculate module sizes for miRNAs
miRNA_module_size_df <- interaction_table %>% 
  distinct(miRNA.cluster, miRNA.module) %>% 
  group_by(miRNA.module) %>% 
  summarise(
    miRNA.module.size = n()
  )

# Calculate the sizes for the gene modules
gene_module_size_df <- interaction_table %>% 
  distinct(genes, gene.module) %>% 
  group_by(gene.module) %>% 
  summarise(
    gene.module.size = n()
  )

# Join the data and normalize for the same pairs
same_count_normalized_df <- miRNA_gene_pairs_df %>% 
  full_join(
    miRNA_module_size_df, by = 'miRNA.module'
  ) %>% 
  full_join(
    gene_module_size_df, by = 'gene.module'
  ) %>% 
  filter(same.module == 'Same') %>% 
  mutate(
    possible.pairs = miRNA.module.size * gene.module.size
  ) %>% 
  # Calculate the normalized value
  mutate(
    norm.count = pairs / possible.pairs
  )

# Join the data and normalize for the different pairs
different_count_normalized_df <- miRNA_gene_pairs_df %>% 
  full_join(
    miRNA_module_size_df, by = 'miRNA.module'
  ) %>% 
  full_join(
    gene_module_size_df, by = 'gene.module'
  ) %>% 
  filter(same.module == 'Different') %>% 
  mutate(
    possible.pairs = miRNA.module.size * gene.module.size
  ) %>% 
  # Calculate the normalized value
  mutate(
    norm.count = pairs / possible.pairs
  )

# Full set of normalized information
norm_count_df <- bind_rows(same_count_normalized_df, different_count_normalized_df)




# Grey WGCNA module removed ----

# Remove the grey module entirely
no_grey_wgcna_df <- wgcna_df %>%
  filter(!module == 'grey')


# Get miRNA modules
no_grey_miRNA_wgcna_df <- no_grey_wgcna_df %>% 
  dplyr::filter(str_detect(genes, 'cvi-|Cluster_')) %>% 
  dplyr::rename(miRNA.cluster = genes, miRNA.module = module)

# Get gene modules
no_grey_gene_wgcna_df <- no_grey_wgcna_df %>%
  dplyr::filter(str_detect(genes, 'Venom')) %>% 
  dplyr::rename(gene.module = module)

# Create an interaction table
# Note that I am using inner join to get instances where all genes are in the data
no_grey_interaction_table <- dplyr::inner_join( 
  # mi_df %>% filter(feature.type == c("three_prime_utr", 'five_prime_utr')),
  mi_df,
  no_grey_gene_wgcna_df,
  by = c('genes')
) %>%
  # Join the miRNA module data, again doing a inner_join so that only miRNAs that target something in the 
  dplyr::inner_join(
    no_grey_miRNA_wgcna_df, by = c('miRNA.cluster')
  ) %>% 
  distinct() %>% 
  # Categorize miRNAs and genes as in the same module or not
  mutate(
    same.module = if_else(gene.module == miRNA.module, "Same", "Different")
  )

# Create nodes
no_grey_nodes <- bind_rows(
  no_grey_interaction_table %>% 
    select(name = miRNA.cluster, color = miRNA.module) %>% 
    mutate(type = 'miRNA'),
  no_grey_interaction_table %>% 
    select(name = genes, color = gene.module) %>% 
    mutate(type = 'Gene')
) %>% 
  distinct() %>% 
  mutate(
    module = str_c(color, type, sep = "_")
  )

# Create edges
no_grey_edges <- no_grey_interaction_table %>% 
  select(from = miRNA.cluster, to = genes, same.module, total.energy) %>% 
  distinct()

# Network data
no_grey_network_graph <- tbl_graph(
  nodes = no_grey_nodes,
  edges = no_grey_edges,
  directed = FALSE
)
glimpse(no_grey_network_graph)

## Hive Plot ----
# Create hive plot
no_grey_wgcna_hive_plot <- ggraph::ggraph(no_grey_network_graph, 'hive', axis = module, sort.by = 'degreee') +
  # Add edges first
  geom_edge_hive(
    aes(
      colour = same.module,
      alpha = total.energy
    ),
    show.legend = TRUE
  ) +
  # Add nodes with module-based coloring
  geom_node_point(
    aes(
      colour = module
    ),
    alpha = 0.8,
    size = 2
  ) +
  # # 81 labels can't be drawn, so I will comment it out
  # geom_node_text(
  #   aes(label = name),
  #   size = 2,
  #   repel = TRUE,
  #   max.overlaps = 10
  # ) +
  scale_edge_color_manual(
    values = edge_colors,
    name = 'Module interaction'
  ) +
  scale_edge_alpha_continuous(
    name = 'Binding Energy (kcal/mol)',
    range = c(0.05, 1),
    trans = 'reverse'
  ) +
  scale_color_manual(
    values = module_colors,
    name = 'Node Module'
  ) +
  coord_fixed() +
  theme_graph() +
  theme(
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5)
  )
# theme_no_axes(base.theme = theme_linedraw())
no_grey_wgcna_hive_plot
ggsave(filename = 'Figures/WGCNA/WGCNA_Hive_plot_no_grey_module_2025.03.31.png', plot = no_grey_wgcna_hive_plot, create.dir = TRUE, width = 10, height = 8)


# With facet
no_grey_wgcna_hive_plot_facet <- no_grey_wgcna_hive_plot + facet_edges(~same.module)
no_grey_wgcna_hive_plot_facet
ggsave(filename = 'Figures/WGCNA/WGCNA_Hive_plot_no_grey_module_facet_2025.03.31.png', plot = no_grey_wgcna_hive_plot_facet, create.dir = TRUE, width = 10, height = 8)


## Quantifying Shared vs Different relationships ----

# Count the number of existing miRNA gene relationships (pairs)
miRNA_gene_pairs_df <- no_grey_interaction_table %>% 
  group_by(miRNA.module, gene.module, same.module) %>% 
  summarise(
    pairs = n()
  )

# Calculate module sizes for miRNAs
miRNA_module_size_df <- no_grey_interaction_table %>% 
  distinct(miRNA.cluster, miRNA.module) %>% 
  group_by(miRNA.module) %>% 
  summarise(
    miRNA.module.size = n()
  )

# Calculate the sizes for the gene modules
gene_module_size_df <- no_grey_interaction_table %>% 
  distinct(genes, gene.module) %>% 
  group_by(gene.module) %>% 
  summarise(
    gene.module.size = n()
  )

# Join the data and normalize for the same pairs
same_count_normalized_df <- miRNA_gene_pairs_df %>% 
  full_join(
    miRNA_module_size_df, by = 'miRNA.module'
  ) %>% 
  full_join(
    gene_module_size_df, by = 'gene.module'
  ) %>% 
  filter(same.module == 'Same') %>% 
  mutate(
    possible.pairs = miRNA.module.size * gene.module.size
  ) %>% 
  # Calculate the normalized value
  mutate(
    norm.count = pairs / possible.pairs
  )

# Join the data and normalize for the different pairs
different_count_normalized_df <- miRNA_gene_pairs_df %>% 
  full_join(
    miRNA_module_size_df, by = 'miRNA.module'
  ) %>% 
  full_join(
    gene_module_size_df, by = 'gene.module'
  ) %>% 
  filter(same.module == 'Different') %>% 
  mutate(
    possible.pairs = miRNA.module.size * gene.module.size
  ) %>% 
  # Calculate the normalized value
  mutate(
    norm.count = pairs / possible.pairs
  )

# Full set of normalized information
norm_count_df <- bind_rows(same_count_normalized_df, different_count_normalized_df)


# Create a box plot that show whether s
norm_count_boxplot <- ggplot(
  norm_count_df,
  aes(
    x = same.module, y = norm.count, fill = same.module
  )
) +
  geom_boxplot() +
  scale_color_manual(values = c('blue', 'red')) + 
  labs(
    y = 'Normalized Count',
    x = 'Relationship',
    fill = 'Relationship'
  ) + 
  ggpubr::stat_compare_means(method = 't.test') +
  theme_linedraw()
norm_count_boxplot
ggsave(filename = 'Figures/WGCNA/Module_Comparison_Boxplot_2025.03.31.png', plot = norm_count_boxplot, create.dir = TRUE, width = 6, height = 8)




# Grey WGCNA module added to other modules ----

# Change the greys to the coexpression module they should belong to
alt_wgcna_df <- wgcna_df %>% 
  # Manually change the clustering
  mutate(
    module = case_when(
      genes == 'Venom_CTL2' ~ 'blue', # Not sure about this one
      genes == 'Cluster_196' ~ 'blue', # Not sure about this one either
      genes == 'Venom_SVSP2' ~ 'blue',
      genes == 'Venom_SVSP6' ~ 'blue',
      genes == 'Cluster_1292' ~ 'blue',
      genes == 'Cluster_807' ~ 'blue',
      genes == 'cvi-miR-23a-5p' ~ 'blue',
      genes == 'Venom_SVSP3' ~ 'blue',
      genes == 'cvi-let-7a-5p' ~ 'blue',
      genes == 'cvi-miR-365a-3p' ~ 'blue',
      genes == 'Venom_SVMP5' ~ 'brown',
      genes == 'Venom_vQC1' ~ 'brown',
      genes == 'Cluster_988' ~ 'brown',
      genes == 'Cluster_1428' ~ 'turquoise',
      genes == 'cvi-miR-737-5p' ~ 'turquoise',
      genes == 'Venom_LAAO2' ~ 'turquoise',
      genes == 'Cluster_919' ~ 'turquoise',
      TRUE ~ module
    )
  )


# Get miRNA modules
alt_miRNA_wgcna_df <- alt_wgcna_df %>% 
  dplyr::filter(str_detect(genes, 'cvi-|Cluster_')) %>% 
  dplyr::rename(miRNA.cluster = genes, miRNA.module = module)

# Get gene modules
alt_gene_wgcna_df <- alt_wgcna_df %>%
  dplyr::filter(str_detect(genes, 'Venom')) %>% 
  dplyr::rename(gene.module = module)

# Create an interaction table
# Note that I am using left join so that only genes that were targeted by something in the 3UTR are kept
alt_interaction_table <- dplyr::left_join( 
  # mi_df %>% filter(feature.type == c("three_prime_utr", 'five_prime_utr')),
  mi_df,
  alt_gene_wgcna_df,
  by = c('genes')
) %>%
  # Join the miRNA module data, again doing a left_join so that only miRNAs that target something in the 
  dplyr::left_join(
    alt_miRNA_wgcna_df, by = c('miRNA.cluster')
  ) %>% 
  distinct() %>% 
  # Categorize miRNAs and genes as in the same module or not
  mutate(
    same.module = if_else(gene.module == miRNA.module, "Same", "Different")
  )

# Create nodes
alt_nodes <- bind_rows(
  alt_interaction_table %>% 
    select(name = miRNA.cluster, color = miRNA.module) %>% 
    mutate(type = 'miRNA'),
  alt_interaction_table %>% 
    select(name = genes, color = gene.module) %>% 
    mutate(type = 'Gene')
) %>% 
  distinct() %>% 
  mutate(
    module = str_c(color, type, sep = "_")
  )

# Create edges
alt_edges <- alt_interaction_table %>% 
  select(from = miRNA.cluster, to = genes, same.module, total.energy) %>% 
  distinct()

# Network data
alt_network_graph <- tbl_graph(
  nodes = alt_nodes,
  edges = alt_edges,
  directed = FALSE
)
glimpse(alt_network_graph)

## Hive Plot ----
# Create hive plot
alt_wgcna_hive_plot <- ggraph::ggraph(alt_network_graph, 'hive', axis = module, sort.by = 'degreee') +
  # Add edges first
  geom_edge_hive(
    aes(
      colour = same.module,
      alpha = total.energy
    ),
    show.legend = TRUE
  ) +
  # Add nodes with module-based coloring
  geom_node_point(
    aes(
      colour = module
    ),
    alpha = 0.8,
    size = 2
  ) +
  # # 81 labels can't be drawn, so I will comment it out
  # geom_node_text(
  #   aes(label = name),
  #   size = 2,
  #   repel = TRUE,
  #   max.overlaps = 10
  # ) +
  scale_edge_color_manual(
    values = edge_colors,
    name = 'Module interaction'
  ) +
  scale_edge_alpha_continuous(
    name = 'Binding Energy (kcal/mol)',
    range = c(0.05, 1),
    trans = 'reverse'
  ) +
  scale_color_manual(
    values = module_colors,
    name = 'Node Module'
  ) +
  coord_fixed() +
  theme_graph() +
  theme(
    legend.position = 'bottom',
    legend.title.position = 'top',
    legend.title = element_text(hjust = 0.5)
  )
# theme_no_axes(base.theme = theme_linedraw())
alt_wgcna_hive_plot

# With facet
alt_wgcna_hive_plot_facet <- alt_wgcna_hive_plot + facet_edges(~same.module)
alt_wgcna_hive_plot_facet

