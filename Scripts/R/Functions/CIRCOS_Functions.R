# This file contains all of the functions for drawing mRNA-miRNA-protein CIRCOS plots

#### Draw CIRCOS Points ####

# Function for drawing points on a track
draw_circos_tracks_points <- function(data, max_value, track_height, chromosome_col, loci_location_col, value_col, color) {
  # Required packages
  require(circlize)
  require(foreach)
  require(dplyr)
  
  # Draw the protein expression data as points
  circos.track(
    ylim = max_value,
    track.height = track_height,
    panel.fun = function(x, y) {
       
      # Filter data for current chromosome
      current_chrom <- data[data[[chromosome_col]] == CELL_META$sector.index, ]
      # data[[chromosome_col]]: Access the column named by the value of chromosome_col in the data frame data.
      # CELL_META$sector.index: Access the value of the sector index in the CELL_META object.
      
      # data[data[[chromosome_col]] == CELL_META$sector.index, ]:
      # Subset the data frame data by selecting rows where the values in the column named by chromosome_col are equal to the value of the sector index in CELL_META.
      
      # current_chrom <- data[data[[chromosome_col]] == CELL_META$sector.index, ]:
      # Store the subsetted data frame in the variable current_chrom.
      
      # Use foreach for iteration over rows
      foreach(loci = 1:nrow(current_chrom)) %do% {
        
        # Variables for x and y in the graph
        horizontal <- current_chrom[[loci_location_col]][loci]
        vertical <- log(current_chrom[[value_col]][loci] + 1)
        
        # Print the above for debugging.
        print(horizontal)
        print(vertical)
        
        # Draw points
        circos.points(
          horizontal,
          vertical,
          col = color,
          pch = 16,
          cex = 0.5
        )
      }
      
      # Add text to CIRCOS about chromosomes
      circos.text(
        CELL_META$xcenter, 
        CELL_META$cell.ylim[2] + mm_y(7),
        gsub('scaffold-', '', CELL_META$sector.index), 
        facing = 'downward',
        cex = 0.6)
      }
    )
}


# Function for drawing points on a track without any text.
draw_circos_tracks_points_without_text <- function(data, max_value, track_height, chromosome_col, loci_location_col, value_col, color) {
  # Required packages
  require(circlize)
  require(foreach)
  require(dplyr)
  
  # Draw the protein expression data as points
  circos.track(
    ylim = max_value,
    track.height = track_height,
    panel.fun = function(x, y) {
      # filter data for current chromosome
      # Create variable to hold the chomosome type (miRNA or miRNA target gene)
      
      # Filter data for current chromosome
      current_chrom <- data[data[[chromosome_col]] == CELL_META$sector.index, ]
      # data[[chromosome_col]]: Access the column named by the value of chromosome_col in the data frame data.
      # CELL_META$sector.index: Access the value of the sector index in the CELL_META object.
      
      # data[data[[chromosome_col]] == CELL_META$sector.index, ]:
      # Subset the data frame data by selecting rows where the values in the column named by chromosome_col are equal to the value of the sector index in CELL_META.
      
      # current_chrom <- data[data[[chromosome_col]] == CELL_META$sector.index, ]:
      # Store the subsetted data frame in the variable current_chrom.
      
      # Use foreach for iteration over rows
      foreach(loci = 1:nrow(current_chrom)) %do% {
        
        # Variables for x and y in the graph
        horizontal <- current_chrom[[loci_location_col]][loci]
        vertical <- log(current_chrom[[value_col]][loci] + 1)
        
        # Print the above for debugging.
        print(horizontal)
        print(vertical)
        
        # Draw points
        circos.points(
          horizontal,
          vertical,
          col = color,
          pch = 16,
          cex = 0.5
        )
      }
    }
  )
}

# Example Usage:
# Note: function must be called in this way:
# Call the function 
# draw_circos_tracks_points(RVG5_circos_df, 
#                   max_value = max_protein_value, 
#                   track_height = 0.1, 
#                   chromosome_col = 'miRNA.Target.Chrom',
#                   loci_location_col =  'miRNA.Target.Start',
#                   value_col = 'Protein.Value.RVG.5S.CV1087.viridis.North.F',
#                   color = 'purple'
#                   )

# draw_circos_tracks_points_without_text(RVG5_circos_df,
#                                 max_value = max_miRNA_count,
#                                 track_height = 0.1,
#                                 chromosome_col = 'miRNA.Sequence.Chrom',
#                                 loci_location_col = 'miRNA.Start',
#                                 value_col = 'miRNA.Counts.RVG.5S.CV1087.viridis.North.F',
#                                 color = 'blue')
# 
# draw_circos_tracks_points_without_text(RVG5_circos_df,
#                                 max_value = max_mRNA_expression,
#                                 track_height = 0.1,
#                                 chromosome_col = 'miRNA.Target.Chrom',
#                                 loci_location_col = 'miRNA.Target.Start',
#                                 value_col = 'RNA.Expression.RVG.5S.CV1087.viridis.North.F',
#                                 color = 'red')





#### Draw CIRCOS Points With Colored Dots ####

# Function for drawing points on a track
draw_circos_tracks_points2 <- function(data, max_value, track_height, chromosome_col, loci_location_col, value_col, color_ramp, chrom_labels = TRUE) {
  # Required packages
  require(circlize)
  require(foreach)
  require(dplyr)
  
  # Draw the protein expression data as points
  circos.track(
    ylim = max_value,
    track.height = track_height,
    panel.fun = function(x, y) {
      
      # Filter data for current chromosome
      current_chrom <- data[data[[chromosome_col]] == CELL_META$sector.index, ]
      # data[[chromosome_col]]: Access the column named by the value of chromosome_col in the data frame data.
      # CELL_META$sector.index: Access the value of the sector index in the CELL_META object.
      
      # data[data[[chromosome_col]] == CELL_META$sector.index, ]:
      # Subset the data frame data by selecting rows where the values in the column named by chromosome_col are equal to the value of the sector index in CELL_META.
      
      # current_chrom <- data[data[[chromosome_col]] == CELL_META$sector.index, ]:
      # Store the subsetted data frame in the variable current_chrom.
      
      # Use foreach for iteration over rows
      foreach(loci = 1:nrow(current_chrom)) %do% {
        
        # Variables for x and y in the graph
        horizontal <- current_chrom[[loci_location_col]][loci]
        vertical <- log(current_chrom[[value_col]][loci] + 1)
        point_color <- color_ramp(vertical)
        
        # Print the above for debugging.
        print(horizontal)
        print(vertical)
        
        # Draw points
        circos.points(
          horizontal,
          vertical,
          col = point_color,
          pch = 16,
          cex = 1
        )
      }
      
      # Add text to CIRCOS about chromosomes
      if (chrom_labels) {
        circos.text(
          CELL_META$xcenter, 
          CELL_META$cell.ylim[2] + mm_y(7),
          gsub('scaffold-', '', CELL_META$sector.index), 
          facing = 'downward',
          cex = 0.8
        )
      }
    }
  )
}










#### Draw CIRCOS arrows ####


# Function for drawing the CIRCOS links
draw_target_arrows <- function(data, sector1_col, position1_col, sector2_col, position2_col, line_color_col, line_type_col) {
  # Iterate through each row of the data frame
  for (cluster in 1:nrow(data)) {
    # Extract necessary information for plotting
    sector1 <- data[[sector1_col]][cluster]
    position1 <- data[[position1_col]][cluster]
    sector2 <- data[[sector2_col]][cluster]
    position2 <- data[[position2_col]][cluster]
    line_and_arrow_color <- data[[line_color_col]][cluster]
    line_type <- data[[line_type_col]][cluster]
    
    # Print for debugging
    print(data[cluster, ])
    
    # Draw links between sectors based on sector and position
    circos.link(sector1, 
                position1, 
                sector2, 
                position2,
                col = line_and_arrow_color,
                # col = add_transparency(line_and_arrow_color, 0.95), # Adjust the line color as needed
                directional = 1,
                lwd = 0.7,
                arr.width = 0.1,
                arr.length = 0.1,
                lty = line_type # Add the line type parameter
    )
  }
}





# Example Usage:
# Note: function must be called in this way:
# Call the function with RVG5_circos_df as argument
# draw_target_arrows(RVG5_circos_df,
#                         sector1_col = 'miRNA.Sequence.Chrom',
#                         position1_col = 'miRNA.Start',
#                         sector2_col = 'miRNA.Target.Chrom',
#                         position2_col = 'miRNA.Target.Start',
#                         line_color_col = 'Color')


# Example of how it used to work:

# # Load CIRCOS plot packages
# library(circlize)
# library(rlist)
# library(officer)
# library(rvg)
# library(foreach)
# library(viridis)
# 
# # Create smaller data frame for circos plot circos.link function to draw lines
# # This data frame also only has RVG_6S information in it and is color coded for venom gene families
# RVG6_venom_circos_df <- mi_df %>% 
#   select(miRNA.Cluster, miRNA.Sequence.Chrom, miRNA.Start, miRNA.End, Genes, miRNA.Target.Chrom, miRNA.Target.Start, miRNA.Target.End, Origin,
#          miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M,
#          RNA.Expression.RVG.6S.CV0987.lutosus.Other.M,
#          Protein.Value.RVG.6S.CV0987.lutosus.Other.M) %>%
#   filter(str_starts(Genes, 'Venom')) %>% # Filter out everything but venom
#   filter(!Origin %in% c('exon', 'five_prime_utr')) %>% # Filter out exons and five_prime_utr as needed
#   mutate(Color = case_when(
#     str_starts(Genes, 'Venom_PLA2') ~ 'purple',
#     str_starts(Genes, 'Venom_SVSP') ~ 'blue',
#     str_starts(Genes, 'Venom_SVMP') ~ 'red',
#     str_starts(Genes, 'Venom_vQC') ~ 'yellow',
#     str_starts(Genes, 'Venom_ohanin') ~ 'green',
#     str_starts(Genes, 'Venom_LAAO') ~ 'orange',
#     str_starts(Genes, 'Venom_CTL') ~ 'magenta',
#     str_starts(Genes, 'Venom_CRISP') ~ 'coral',
#     str_starts(Genes, 'Venom_BPP') ~ 'cyan',
#     str_starts(Genes, 'Venom_VEGF') ~ 'navy',
#     str_starts(Genes, 'Venom_myotoxin') ~ 'maroon',
#     str_starts(Genes, 'Venom_BPP') ~ 'gray',
#     str_starts(Genes, 'Venom_EXO') ~ 'pink',
#     str_starts(Genes, 'Venom_ADAM') ~ 'aquamarine',
#     TRUE ~ 'black'
#   )) %>% 
#   #filter(!Origin %in% c('exon', 'five_prime_utr')) %>% # Filter out exons as needed
#   filter(!miRNA.Target.Chrom %in% c('PE-reconstructed-10x-myo', 'scaffold-un187')) %>% 
#   filter(!miRNA.Sequence.Chrom %in% c('scaffold-un11', 'scaffold-un619', 'scaffold-un31', 'scaffold-un147')) %>% # These are weird unplaced scaffolds that don't correspond to anything and unfortunately myotoxin and Venom_BPP are on them
#   # Also note that BPP is on scaffold-un187, and there is a cluster (Cluster_2557) that target myotoxin, PLA2B1, and Venom_SVMP6 through 9 is on scaffold-un11. 
#   distinct()
# 
# # Load scaffold sizes
# # This will be used to initiate the sizes of each sector of the
# # circos plot
# scaffold_size = read.delim('Data/Misc/scaffold_sizes.txt', header = F) %>% 
#   select(-V4, -V5)
# 
# # rename columns in scaffold_size
# names(scaffold_size) = c('Chrom','size','genome_position')
# 
# # create a column of zeros, to indicate the starting
# # value of each sector
# scaffold_size = scaffold_size %>% mutate(start = 0)
# 
# # Begin saving the plot as a pdf
# pdf(file = 'Figures/CIRCOS/Venom_Genes_Colored_CIRCOS_3UTR_only_lutosus.pdf', width = 10, height = 10)
# 
# circos.clear()
# 
# # Set basic parameters for the circos plot
# circos.par('track.height' = 0.1,
#            cell.padding = c(0.01, 0, 0.01, 0),
#            start.degree = 90,
#            gap.degree = c(rep(1, 17), 20))
# 
# # Initiate the circos plot
# # xlim is used to specify the start and the end value of each chromosome (sectors)
# circos.initialize(scaffold_size$Chrom, xlim = cbind(scaffold_size$start, scaffold_size$size))
# 
# # Initialize max count values
# max_protein_value <- range(log(RVG6_venom_circos_df$Protein.Value.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)
# max_miRNA_count <- range(log(RVG6_venom_circos_df$miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)
# max_mRNA_expression <- range(log(RVG6_venom_circos_df$RNA.Expression.RVG.6S.CV0987.lutosus.Other.M + 1), na.rm = T)
# 
# # Draw the protein expression data as points
# circos.track(ylim = max_protein_value, track.height = 0.1,
#              panel.fun = function(x, y) {
#                # filter data for current chromosome
#                # Create variable to hold the chomosome type (miRNA or miRNA target gene)
#                miRNA_gene_chrom <- RVG6_venom_circos_df$miRNA.Target.Chrom
#                
#                current_chrom <-  RVG6_venom_circos_df %>% 
#                  filter(miRNA_gene_chrom == CELL_META$sector.index)
#                
#                # Use foreach for iterating over rows
#                foreach(protein_target = 1:nrow(current_chrom)) %do% {
#                  
#                  # Variables for x and y
#                  protein_start <- current_chrom$miRNA.Target.Start[protein_target]
#                  protein_value <- log(current_chrom$Protein.Value.RVG.6S.CV0987.lutosus.Other.M[protein_target] + 1) # Log scale to make things more visible
#                  
#                  print(protein_start)
#                  print(protein_value)
#                  
#                  circos.points(protein_start, 
#                                protein_value, 
#                                col = "purple", 
#                                pch = 16, 
#                                cex = 0.5)
#                }
#                
#                circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + mm_y(7), 
#                            gsub('scaffold-', '', CELL_META$sector.index), 
#                            facing = 'downward',
#                            cex = 0.6)
#              })
# 
# # Draw the miRNA expression data as points
# circos.track(ylim = max_miRNA_count, track.height = 0.1,
#              panel.fun = function(x, y) {
#                # filter data for current chromosome
#                # Create variable to hold the chomosome type (miRNA or miRNA target gene)
#                miRNA_gene_chrom <- RVG6_venom_circos_df$miRNA.Sequence.Chrom
#                
#                current_chrom <-  RVG6_venom_circos_df %>% 
#                  filter(miRNA_gene_chrom == CELL_META$sector.index)
#                
#                # Use foreach for iterating over rows
#                foreach(miRNA_cluster = 1:nrow(current_chrom)) %do% {
#                  
#                  # Variables for x and y
#                  miRNA_start <- current_chrom$miRNA.Start[miRNA_cluster]
#                  miRNA_count <- log(current_chrom$miRNA.Counts.RVG.6S.CV0987.lutosus.Other.M[miRNA_cluster] + 1)
#                  
#                  print(miRNA_start)
#                  print(miRNA_count)
#                  
#                  circos.points(miRNA_start, 
#                                miRNA_count, 
#                                col = "blue", 
#                                pch = 16, 
#                                cex = 0.5)
#                }
#              })
# 
# # Draw the miRNA expression data as points
# circos.track(ylim = max_mRNA_expression, track.height = 0.1,
#              panel.fun = function(x, y) {
#                # filter data for current chromosome
#                # Create variable to hold the chomosome type (miRNA or miRNA target gene)
#                miRNA_gene_chrom <- RVG6_venom_circos_df$miRNA.Target.Chrom
#                
#                current_chrom <-  RVG6_venom_circos_df %>% 
#                  filter(miRNA_gene_chrom == CELL_META$sector.index)
#                
#                # Use foreach for iterating over rows
#                foreach(mRNA_target = 1:nrow(current_chrom)) %do% {
#                  
#                  # Variables for x and y
#                  mRNA_start <- current_chrom$miRNA.Target.Start[mRNA_target]
#                  mRNA_expression <- log(current_chrom$RNA.Expression.RVG.6S.CV0987.lutosus.Other.M[mRNA_target] + 1)
#                  
#                  print(mRNA_start)
#                  print(mRNA_expression)
#                  
#                  circos.points(mRNA_start, 
#                                mRNA_expression, 
#                                col = "red", 
#                                pch = 16, 
#                                cex = 0.5)
#                }
#              })
# 
# # Plot each miRNA target line one by one with a for loop
# for(cluster in 1:nrow(RVG6_venom_circos_df)) {
#   
#   # Need to specify sector1, position1 and sector2, position2 to plot the miRNA targeting line
#   sector1 <- RVG6_venom_circos_df$miRNA.Sequence.Chrom[cluster]
#   position1 <- RVG6_venom_circos_df$miRNA.Start[cluster]
#   sector2 <- RVG6_venom_circos_df$miRNA.Target.Chrom[cluster]
#   position2 <- RVG6_venom_circos_df$miRNA.Target.Start[cluster]
#   line_and_arrow_color <- RVG6_venom_circos_df$Color[cluster]
#   
#   # Print for debugging
#   print(RVG6_venom_circos_df[cluster,])
#   
#   # Draw links between sectors based on sector and position.
#   circos.link(sector1, 
#               position1, 
#               sector2, 
#               position2,
#               col = line_and_arrow_color,
#               #col = add_transparency(line_and_arrow_color, 0.95), # Adjust the line color as needed
#               directional = 1,
#               lwd = 0.3,
#               arr.width = 0.1,
#               arr.length = 0.1
#   )
# }
# 
# # Finish pdf
# dev.off()





