# Last Edited: 2025/02/02

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
library(ggvenn)
library(ggraph)
library(igraph)
library(tidygraph)
library(ggtree)
library(ape)
library(UpSetR)
library(ComplexUpset)
library(scales)
library(arrow)

## Read data ----
# Create variable for the fused dataset.
miRNA_mRNA_protein_data <- 'Data/Merged/mRNA_Protein_miRNA_Combined_Data_Venom_2025.01.22.parquet'
# Create a variable for the mRNA and protein data
mRNA_protein_data <- 'Data/Merged/mRNA_Protein_Combined_Data_2025.01.22.parquet'

# Read both in as data frames
miRNA_mRNA_protein_df <- read_parquet(file = miRNA_mRNA_protein_data)

# Create shorter df name and do some minor tweaks to it's structure for readability
mi_df <- miRNA_mRNA_protein_df %>% 
  filter(
    !sample.id == 'CV1082_viridis',
    str_detect(genes, 'Venom_'),
    !str_detect(genes, 'ADAM'),
    feature.type == 'three_prime_utr'
  ) %>% 
  distinct() %>% 
  # Change case of family names
  mutate(
    venom.family = case_when(
      grepl('SVMP', genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('VEGF', genes) ~ 'VEGF',
      grepl('ohanin', genes) ~ 'Ohanin',
      grepl('vQC', genes) ~ 'vQC',
      grepl('SVSP', genes) ~ 'SVSP',
      grepl('PLA2', genes) ~ 'PLA2',
      grepl('CRISP', genes) ~ 'CRISP',
      grepl('CTL', genes) ~ 'CTL',
      grepl('EXO', genes) ~ 'EXO',
      grepl('LAAO', genes) ~ 'LAAO',
      grepl('myotoxin', genes) ~ 'Myotoxin',
      grepl('BPP', genes) ~ 'bpp',
      TRUE ~ 'others'
    )
  )
glimpse(mi_df)
rm(miRNA_mRNA_protein_df)

# Read the protein/mRNA data
genes_df <- read_parquet(mRNA_protein_data) %>%
  filter(
    !sample.id == 'CV1082_viridis',
    str_detect(genes, 'Venom_'),
    !str_detect(genes, 'ADAM')
  ) %>%
  # Remove genes if they already exist in the mi_df
  filter(!genes %in% mi_df$genes) %>% 
  # Change case of family names
  mutate(
    venom.family = case_when(
      grepl('SVMP', genes) ~ 'SVMP', # Add a Venom.Families Column
      grepl('VEGF', genes) ~ 'VEGF',
      grepl('ohanin', genes) ~ 'Ohanin',
      grepl('vQC', genes) ~ 'vQC',
      grepl('SVSP', genes) ~ 'SVSP',
      grepl('PLA2', genes) ~ 'PLA2',
      grepl('CRISP', genes) ~ 'CRISP',
      grepl('CTL', genes) ~ 'CTL',
      grepl('EXO', genes) ~ 'EXO',
      grepl('LAAO', genes) ~ 'LAAO',
      grepl('myotoxin', genes) ~ 'Myotoxin',
      grepl('BPP', genes) ~ 'bpp',
      TRUE ~ 'others'
    )
  )



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
miRNA_color <- 'grey14'
novel_miRNA_color <- 'grey'
three_prime_color <- 'black'
# five_prime_color <- '#0072b2'
five_prime_color <- '#1B9E77'
# cds_color <- '#d55e00'
cds_color <- '#4A70B5'

# Create an order for the genes
venom_gene_order <- c(
  'BPP',
  'myotoxin',
  'ohanin',
  'PLA2B1', 'PLA2K', 'PLA2C1', 'PLA2A1',
  'SVSP1', 'SVSP2', 'SVSP3', 'SVSP10', 'SVSP6', 'SVSP11', 'SVSP7', 'SVSP8', 'SVSP9',
  'SVMP1', 'SVMP2', 'SVMP3', 'SVMP4', 'SVMP5', 'SVMP6', 'SVMP7', 'SVMP8', 'SVMP9', 'SVMP10', 'SVMP11', 'SVMP12',
  'CRISP1', 'CRISP2', 'CRISP3', 'CRISP4',
  'CTL1', 'CTL2', 'CTL3', 'CTL4', 'CTL5', 'CTL6',
  'EXO1', 'EXO2', 'EXO3',
  'LAAO1', 'LAAO2', 'LAAO3',
  'VEGF1', 'VEGF2',
  'vQC1', 'vQC2'
) 

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

# Create color scheme for the 
loci_colors <- c(
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
  bpp = BPP_color,
  others = other_color,
  novel_miRNA = novel_miRNA_color,
  miRNA = miRNA_color
)


# Hierachical edge modeling ----

# Bind rows
mi_df2 <- bind_rows(mi_df, genes_df)

# Remove the Venom_ prefix
mi_df2 <- mi_df2 %>% 
  mutate(genes = str_remove(genes, '^Venom_'))
glimpse(mi_df2)


# Create a data frame with the connections between miRNAs and genes
connect <- mi_df2 %>% 
  distinct(miRNA.cluster, genes) %>% 
  rename(
    from = miRNA.cluster, to = genes
  ) %>%
  select(from, to) # Rearrange
rownames(connect) <- NULL


# Get miRNA edges
mi_edges <- mi_df2 %>% 
  distinct(
    miRNA.cluster
  ) %>% 
  rename(to = miRNA.cluster) %>% 
  mutate(
    from = case_when(
      grepl('Cluster_', to) ~ 'novel_miRNA',
      grepl('cvi-', to) ~ 'miRNA', 
      TRUE ~ 'others'
    )
  ) %>% 
  filter(!is.na(to))

# Get gene edges
gene_edges <- mi_df2 %>% 
  distinct(genes, venom.family) %>% rename(to = genes, from = venom.family)


# Create a data frame containing all of the information need to create the hierachical edge modeling graph
edges <- rbind(mi_edges, gene_edges)

# Get group edges
group_edges <- edges %>% 
  distinct(from) %>% 
  rename(to = from) %>% 
  mutate(from = 'origin')

# Update edges df
edges <- rbind(group_edges, edges) %>% select(from, to)
rownames(edges) <- NULL

# Create a data frame that contains only miRNAs and their expression levels
mirna_df <- mi_df2 %>% 
  select(
    miRNA.cluster, miRNA.vst
  ) %>%
  distinct() %>% 
  group_by(miRNA.cluster) %>% 
  summarise(
    expression = mean(miRNA.vst)
  ) %>% 
  rename(name = miRNA.cluster)

# Create a data frame containing only genes and their expression levels
genes_df <- mi_df2 %>% 
  select(
    genes, mRNA.vst
  ) %>% 
  distinct() %>% 
  group_by(genes) %>% 
  summarise(
    expression = mean(mRNA.vst)
  ) %>% 
  rename(name = genes)

# Create a data frame containing expression in vst
expression_df <- rbind(mirna_df, genes_df) %>% filter(!is.na(name))
rownames(expression_df) <- NULL

# Create a vertices data frame
vertices <- data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to)))
) %>% 
  left_join(
    expression_df,
    by = 'name'
  )
vertices$group <- edges$from[match(vertices$name, edges$to)]

# Arrange vertices by group order
vertices <- vertices %>% arrange(group, name)


#Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
#calculate the ANGLE of the labels
vertices$id <- NA
myleaves1 <- which(is.na(match(vertices$name, edges$from)))
nleaves1 <- length(myleaves1)
vertices$id[myleaves1] <- seq(1:nleaves1)
vertices$angle <- 90 - 360 * vertices$id / nleaves1

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust <- ifelse(vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle + 180, vertices$angle)


# Check if all edge nodes exist in vertices
all_nodes_in_vertices <- all(unique(c(edges$to, edges$from)) %in% vertices$name)
print(all_nodes_in_vertices)  # Should return TRUE


# Create a graph object
graph_data <- igraph::graph_from_data_frame(edges, vertices = vertices)

# Convert to tidy graph
graph_tidy <- as_tbl_graph(graph_data)
glimpse(graph_tidy)

# Acces node data as a tibble
nodes_df <- graph_tidy %>% 
  activate(nodes) %>% 
  as_tibble()

# Access edge data as a tibble (data frame)
edges_df <- graph_tidy %>%
  activate(edges) %>%
  as_tibble()

# The connection object must refer to the ids of the leaves:
miRNA = match(connect$from, vertices$name)
gene = match(connect$to, vertices$name)

vertices

# Create the hierarchical network plot
hem_plot <- ggraph::ggraph(graph_data, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = miRNA, to = gene), alpha = 0.2, width = 0.5, aes(colour = after_stat(index))) +
  scale_edge_colour_distiller(palette = "RdPu") +
  
  geom_node_text(aes(x = x * 1.15, y = y * 1.15, filter = leaf, label = name, angle = angle, hjust = hjust, colour = group), size = 2, alpha = 1) +
  
  geom_node_point(aes(filter = leaf, x = x * 1.07, y = y * 1.07, colour = group, size = expression, alpha = 0.2)) +
  scale_color_manual(values = loci_colors) +
  scale_size_continuous(range = c(0.05, 7)) +
  
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
hem_plot
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Hierachical_Network/miRNA-Venom_gene_network_2025.02.02.pdf', plot = hem_plot, create.dir = T, width = 10, height = 10, dpi = 900)

# Shared miRNA dendrograms ----
## Create dendrogram of venom genes by shared miRNAs ----

# Create dendrogram data
dendrogram_df <- mi_df %>% 
  mutate(genes = str_remove(genes, '^Venom_')) %>% 
  distinct(genes, miRNA.cluster) %>% 
  mutate(Value = 1) %>% 
  pivot_wider(
    names_from = miRNA.cluster,
    values_from = Value,
    values_fill = 0
  )
glimpse(dendrogram_df)
class(dendrogram_df)

# Set the data frame to be a regular data frame
dendrogram_df <- as.data.frame(dendrogram_df)
class(dendrogram_df)

# Set rownames to the genes column
rownames(dendrogram_df) <- dendrogram_df$genes
# Remove the genes column
dendrogram_df <- dendrogram_df %>% select(-genes)

# Perform hierarchical clustering
hc <- hclust(dist(dendrogram_df), method = 'ward.D2')

# Convert the clustering to a tree structure
mi_tree <- as.phylo(hc)


# Create dendrogram with ggtree
gene_dendrogram <- ggtree(mi_tree, layout = 'rectangular', size = 0.8) +
  geom_tiplab(size = 4, angle = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))  # Shrink horizontal axis spacing
gene_dendrogram
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Dendrogram/miRNA-Venom_gene_dendrogram_2025.02.02.pdf', plot = gene_dendrogram, create.dir = T, width = 10, height = 10, dpi = 900)



## Create dendrogram of venom families by shared miRNAs ----

# Create dendrogram data
dendrogram_df2 <- mi_df %>% 
  mutate(genes = str_remove(genes, '^Venom_')) %>% 
  distinct(venom.family, miRNA.cluster) %>% 
  mutate(Value = 1) %>% 
  pivot_wider(
    names_from = miRNA.cluster,
    values_from = Value,
    values_fill = 0
  )
glimpse(dendrogram_df2)
class(dendrogram_df2)

# Set the data frame to be a regular data frame
dendrogram_df2 <- as.data.frame(dendrogram_df2)
class(dendrogram_df2)

# Set rownames to the genes column
rownames(dendrogram_df2) <- dendrogram_df2$venom.family
# Remove the genes column
dendrogram_df2 <- dendrogram_df2 %>% select(-venom.family)

# Perform hierarchical clustering
hc2 <- hclust(dist(dendrogram_df2), method = 'ward.D2')

# Convert the clustering to a tree structure
mi_tree2 <- as.phylo(hc2)


# Create dendrogram with ggtree
family_dendrogram <- ggtree(mi_tree2, layout = 'rectangular', size = 0.8) +
  geom_tiplab(size = 4, angle = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.2)))  # Shrink horizontal axis spacing
family_dendrogram
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Dendrogram/miRNA-Venom_family_dendrogram_2025.02.02.pdf', plot = family_dendrogram, create.dir = T, width = 10, height = 10, dpi = 900)



# Venn Diagrams ----
## Set up data for venn diagram ----

# Set up data frame to be able to use ggvenn
venn_df <- mi_df %>% 
  select(
    miRNA.cluster, venom.family
  ) %>% 
  distinct() %>% 
  pivot_wider(
    names_from = venom.family,
    values_from = venom.family,
    values_fn = list(venom.family = ~ TRUE),
    values_fill = FALSE
  ) 
rownames(venn_df) <- NULL
glimpse(venn_df)

## Venn 1 ----
# Create a test venn diagram
venn_1 <- ggplot(
  venn_df, aes(
    A = `Myotoxin`, B = `SVMP`, C = `SVSP`, D = `PLA2`
  )
) +
  geom_venn(fill_color = c(myotoxin_color, SVMP_color, SVSP_color, PLA2_color)) + theme_void() + coord_fixed() +
  labs(
    title = 'miRNAs shared between Myotoxin, PLA2s, SVMPs, and SVSPs'
  ) +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', size = 15, hjust = 0)
  )
venn_1
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Venn_Diagram/Myotoxin_PLA2_SVMP_SVSP_2025.02.02.pdf', plot = venn_1, create.dir = T, width = 10, height = 10, dpi = 900)

## Venn 2 ----
venn_2 <- ggplot(
  venn_df, aes(
    A = `CTL`, B = `EXO`, C = `CRISP`, D = `SVSP`
  )
) +
  geom_venn(fill_color = c(CTL_color, EXO_color, CRISP_color, SVSP_color)) + theme_void() + coord_fixed() +
  labs(
    title = 'miRNAs shared between CTLs, EXOs, CRISPs, and SVSPs'
  ) +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', size = 15, hjust = 1)
  )
venn_2
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Venn_Diagram/CTL_EXO_CRISP_SVSP_2025.02.02.pdf', plot = venn_2, create.dir = T, width = 10, height = 10, dpi = 900)

## Venn 3 ----
venn_3 <- ggplot(
  venn_df, aes(
    A = `EXO`, B = `CRISP`, C = `SVMP`, D = `SVSP`
  )
) +
  geom_venn(fill_color = c(EXO_color, CRISP_color, SVMP_color, SVSP_color)) + theme_void() + coord_fixed() +
  labs(
    title = 'miRNAs shared between EXOs, CRISPs, SVMPs, and SVSPs'
  ) +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', size = 15, hjust = 1)
  )
venn_3
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Venn_Diagram/EXO_CRISP_SVMP_SVSP_2025.02.02.pdf', plot = venn_3, create.dir = T, width = 10, height = 10, dpi = 900)

## Venn 4 ----
venn_4 <- ggplot(
  venn_df, aes(
    A = `Ohanin`, B = `vQC`, C = `Myotoxin`, D = `PLA2`
  )
) +
  geom_venn(fill_color = c(ohanin_color, vQC_color, myotoxin_color, PLA2_color)) + theme_void() + coord_fixed() +
  labs(
    title = 'miRNAs shared between Ohanin, vQCs, Myotoxin, and PLA2s'
  ) +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', size = 15, hjust = 1)
  )
venn_4
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Venn_Diagram/Ohanin_vQC_Myotoxin_PLA2_2025.02.02.pdf', plot = venn_4, create.dir = T, width = 10, height = 10, dpi = 900)

## Venn 5 ----
venn_5 <- ggplot(
  venn_df, aes(
    A = `VEGF`, B = `vQC`, C = `Myotoxin`, D = `PLA2`
  )
) +
  geom_venn(fill_color = c(VEGF_color, vQC_color, myotoxin_color, PLA2_color)) + theme_void() + coord_fixed() +
  labs(
    title = 'miRNAs shared between VEGFs, vQCs, Myotoxins, and PLA2s'
  ) +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', size = 15, hjust = 1)
  )
venn_5
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Venn_Diagram/VEGF_vQC_Myotoxin_PLA2_2025.02.02.pdf', plot = venn_5, create.dir = T, width = 10, height = 10, dpi = 900)

## Venn 6 ----
venn_6 <- ggplot(
  venn_df, aes(
    A = `EXO`, B = `SVSP`, C = `Myotoxin`, D = `VEGF`
  )
) +
  geom_venn(fill_color = c(EXO_color, SVSP_color, myotoxin_color, VEGF_color)) + theme_void() + coord_fixed() +
  labs(
    title = 'miRNAs shared between EXOs, SVSPs, Myotoxins, and VEGFs'
  ) +
  theme(
    plot.title = element_text(color = 'black', face = 'bold', size = 15, hjust = 1)
  )
venn_6
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Venn_Diagram/EXO_SVSP_Myotoxin_VEGF_2025.02.02.pdf', plot = venn_6, create.dir = T, width = 10, height = 10, dpi = 900)




# Up Set Plot ----

# Set up data
set_df <- mi_df %>% 
  select(miRNA.cluster, venom.family) %>% 
  distinct() %>% 
  pivot_wider(
    names_from = venom.family,
    values_from = venom.family,
    values_fn = ~ 1,  # Fill with 1 if the miRNA.cluster belongs to the venom.family
    values_fill = 0   # Fill with 0 otherwise
  ) %>% 
  column_to_rownames(var = "miRNA.cluster")  # Set miRNA.cluster as rownames

# Create upset graph
up_set_plot <- UpSetR::upset(
  set_df, 
  sets = colnames(set_df), 
  nintersects = NA, 
  order.by = c('freq')
)
up_set_plot



# Create the plot with ggplot2-style color mapping
up_set_plot2 <- ComplexUpset::upset(
  # Set data
  set_df,
  # Set what intersections to use
  intersect = colnames(set_df),
  # Set intersection size colors and text color
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = 'bars_color'),
      text_colors = c(
              on_background = 'black', on_bar = 'black'
            )
    ) + 
      scale_fill_manual(values = c('bars_color' = 'grey14'), guide = 'none') +
      labs(y = "Number of miRNAs") # Change Intersection size label 
  ),
  # Set width ratio
  width_ratio = 0.25,
  # Set stripe color for the dots
  stripes = c(
    'white'
  ),
  # Set name under dots
  name = 'Venom Family', 
  # # Sort intesections by degree
  # sort_intersections_by = 'degree',
  # # Group intersections by set
  # group_by = 'sets',
  # Set venom family colors
  queries = list(
    # Make the SVMPs red
    upset_query(
      intersect = c('SVMP'),
      color = 'coral',
      fill = 'coral',
      only_components = c('Intersection size')
    ),
    # Each of these queries gives the color for each gene family
    upset_query(
      set = 'SVMP',
      fill = SVMP_color
    ),
    upset_query(
      set = 'CTL',
      fill = CTL_color
    ),
    upset_query(
      set = 'SVSP',
      fill = SVSP_color
    ),
    upset_query(
      set = 'EXO',
      fill = EXO_color
    ),
    upset_query(
      set = 'CRISP',
      fill = CRISP_color
    ),
    upset_query(
      set = 'VEGF',
      fill = VEGF_color
    ),
    upset_query(
      set = 'Myotoxin',
      fill = myotoxin_color
    ),
    upset_query(
      set = 'vQC',
      fill = vQC_color
    ),
    upset_query(
      set = 'PLA2',
      fill = PLA2_color
    ),
    upset_query(
      set = 'Ohanin',
      fill = ohanin_color
    )
  ),
  # Display counts
  set_sizes = (
    upset_set_size() +
      geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = 'count') +
      theme(
        axis.text.x = element_text(angle = 0)
        # axis.title.x = element_text(size = 10),   # Modify size for axis title
        # axis.text = element_text(size = 8),      # Modify size for tick labels
        # plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
      )
  )
)
up_set_plot2
ggsave('Figures/miRNA-Gene_Relationships/3UTR/UpSetPlot/Venom_family_UpSetPlot_2025.06.20.pdf', plot = up_set_plot2, create.dir = T, width = 20, height = 10, dpi = 900)


# Create a second set df that contains the genes instead of families
gene_set_df <- mi_df %>% 
  select(miRNA.cluster, genes) %>% 
  distinct() %>% 
  # filter(!(genes %in% c('Venom_ohanin', 'Venom_vQC1', 'Venom_VEGF1'))) %>% 
  pivot_wider(
    names_from = genes,
    values_from = genes,
    values_fn = ~ 1,  # Fill with 1 if the miRNA.cluster belongs to the venom.family
    values_fill = 0   # Fill with 0 otherwise
  ) %>% 
  column_to_rownames(var = "miRNA.cluster")  # Set miRNA.cluster as rownames


# Gene upset plot
gene_up_set_plot <- ComplexUpset::upset(
  # Set data
  gene_set_df,
  # Set what intersections to use
  intersect = colnames(gene_set_df),
  # Set intersection size colors and text color
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = 'bars_color'),
      text_colors = c(
        on_background = 'black', on_bar = 'white'
      )
    ) + 
      scale_fill_manual(values = c('bars_color' = 'grey14'), guide = 'none') +
      labs(y = "Number of miRNA Clusters") # Change Intersection size label 
  ),
  # Set width ratio
  width_ratio = 0.25,
  # Set stripe color for the dots
  stripes = c(
    'white'
  ),
  # Set name under dots
  name = 'genes', 
  # Display counts
  set_sizes = (
    upset_set_size() +
      geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = 'count') +
      theme(
        axis.text.x = element_text(angle = 0)
        # axis.title.x = element_text(size = 10),   # Modify size for axis title
        # axis.text = element_text(size = 8),      # Modify size for tick labels
        # plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
      )
  )
)
gene_up_set_plot
ggsave('Figures/miRNA-Gene_Relationships/3UTR/UpSetPlot/Venom_Gene_UpSetPlot_2025.02.02.pdf', plot = gene_up_set_plot, create.dir = T, width = 20, height = 10, dpi = 900)


# Heat map of shared miRNAs per gene ----

# Format data so that the number of shared miRNAs can be used to make a heatmap
genes_per_mirna <- mi_df %>% 
  select(genes, miRNA.cluster) %>% 
  distinct() %>%
  group_by(miRNA.cluster) %>% 
  summarize(genes = list(genes), .groups = 'drop')


# Save list of genes
targeted_genes <- unique(mi_df$genes)
targeted_genes

# Create an emty matrix that will be filled with the number of miRNAs
gene_matrix <- matrix(
  0, nrow = length(targeted_genes),
  ncol = length(targeted_genes),
  dimnames = list(targeted_genes, targeted_genes)
)

# This one doesn't count miRNAs shared between genes themselves
# # For loop that will populate the matrix
# for (row in seq_len(nrow(genes_per_mirna))) {
#   genes <- genes_per_mirna$genes[[row]]
#   for (gene1 in genes) {
#     for (gene2 in genes) {
#       if (gene1 != gene2) {  # Skip self-pairing
#         gene_matrix[gene1, gene2] <- gene_matrix[gene1, gene2] + 1
#       }
#     }
#   }
# }

# This one does count them between themselves
# For loop that will populate the matrix
for (row in seq_len(nrow(genes_per_mirna))) {
  genes <- genes_per_mirna$genes[[row]]
  for (gene1 in genes) {
    for (gene2 in genes) {
      gene_matrix[gene1, gene2] <- gene_matrix[gene1, gene2] + 1
    }
  }
}

# load reshape
library(reshape2)

# Create data frame for ggplot to use
shared_mirnas_df <-  melt(gene_matrix, varnames = c('Gene1', 'Gene2'), value.name = 'Shared.miRNAs') %>% 
  mutate(
    Gene1 = str_remove(Gene1, '^Venom_'),
    Gene2 = str_remove(Gene2, '^Venom_')
  )

# Enforce venom gene order
shared_mirnas_df$Gene1 <- factor(shared_mirnas_df$Gene1, levels = venom_gene_order)
shared_mirnas_df$Gene2 <- factor(shared_mirnas_df$Gene2, levels = venom_gene_order)

# view(shared_mirnas_df)


# Create heatmap
shared_mirnas_heatmap <- ggplot(
  shared_mirnas_df,
  aes(x = Gene1, y = Gene2, fill = Shared.miRNAs)
) + 
  geom_tile() +
  scale_fill_viridis_c(option = 'magma', limits = c(0, 15), oob = scales::squish) +
  # scale_fill_gradient(
  #   low = 'grey',
  #   high = scales::muted('red')
  # ) +
  labs(
    title = 'miRNA shared by gene',
    fill = 'miRNAs shared'
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1),
    axis.title = element_blank(),
    plot.title = element_text(face = 'bold')
  )
shared_mirnas_heatmap
ggsave('Figures/miRNA-Gene_Relationships/3UTR/Heat_Maps/miRNAs_shared_by_gene_2025.02.02.pdf', plot = shared_mirnas_heatmap, create.dir = T, width = 12, height = 10, dpi = 900)
