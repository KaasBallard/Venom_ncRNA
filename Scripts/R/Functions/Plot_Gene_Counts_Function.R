# Last edited: 2025/02/07

# DESeq2 gene plots function ----
# Define the function
plot_gene_counts <- function(dds, gene, intgroup, saveDir, date) {
  # Require the following packages
  require(ggplot2)
  require(dplyr)
  require(DESeq2)
  
  # Extract normalized counts for the specified gene
  gene_df <- plotCounts(dds, gene = gene, intgroup = intgroup, returnData = TRUE)
  
  # Calculate mean and standard deviation for each group
  gene_stats <- gene_df %>%
    group_by(!!sym(intgroup)) %>%
    summarise(
      mean_count = mean(count),
      sd_count = sd(count),
      n = n(),
      se_count = sd_count / sqrt(n)  # Standard error for variance intervals
    ) %>% 
    mutate(
      color = case_when(
        grepl('Crotalus_oreganus', !!sym(intgroup)) ~ oreganus_color,
        grepl('Crotalus_viridis', !!sym(intgroup)) ~ viridis_color
      )
    ) %>% 
    mutate(
      species = case_when(
        grepl('Crotalus_oreganus', !!sym(intgroup)) ~ 'C. oreganus',
        grepl('Crotalus_viridis', !!sym(intgroup)) ~ 'C. viridis'
      )
    )
  
  # Set the colors 
  colors <- setNames(gene_stats$color, gene_stats$species)
  
  # Create the plot with mean and error bars for variance intervals
  gene_plot <- ggplot(gene_stats, aes(x = !!sym(intgroup), y = mean_count, fill = !!sym(intgroup))) +
    geom_bar(stat = 'identity', width = 0.25, color = 'black', linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.2) +
    scale_fill_manual(values = colors) +
    labs(y = 'Mean Normalized Counts', x = intgroup, title = gene) +
    theme_classic2() +
    theme(
      plot.title = element_text(hjust = 0.5, margin = margin(b = 5, t = 5), size = 10),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      legend.position = 'none'
    )
  
  # Create the filename
  file_name <- paste(gene, "DESeq_Mean_Variance_Plot", date, "pdf", sep = ".")
  
  # Save the plot as a PDF file
  ggsave(filename = file.path(saveDir, file_name), plot = gene_plot, device = "pdf", create.dir = TRUE)
  
  # Return the plot
  return(gene_plot)
}