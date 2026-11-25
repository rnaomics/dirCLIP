#!/usr/bin/env Rscript

# Figure: Read length distribution visualization
# Author: Mohamed Bahrudeen

# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

#--------------------------------------------------
# Draw violin plot of read lengths
#--------------------------------------------------
draw_violin_plot <- function(read_lengths_df,
                             group_col,
                             sample_names,
                             y_tick_interval,
                             Colors,
                             text_size,
                             x_lab,
                             y_lab,
                             Title) {
  
  # Generate sample-level color mapping when plotting by Condition
  if (group_col == "Condition") {
    
    condition_color_map <- unique(read_lengths_df[, c("Condition", "Pool_condition")])
    
    condition_color_map$Colors <-Colors[condition_color_map$Pool_condition]
    
    Colors <- condition_color_map$Colors
    names(Colors) <- condition_color_map$Condition
  }
  
  # Standardize grouping column name
  colnames(read_lengths_df)[colnames(read_lengths_df) == group_col] <- "group"
  
  read_lengths_df$group <- factor(read_lengths_df$group,unique(sample_names))
  Angle = 0 
  violin_plot <- ggplot(read_lengths_df, aes(x = group, y = x, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.1, fill = "white") +
    scale_fill_manual(values = Colors) +
    labs(
      x = x_lab,
      y = y_lab,
      title = Title
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      # General text size
      text = element_text(size = text_size),
      
      # Panel appearance
      panel.grid = element_blank(),
      panel.border = element_blank(),
      
      # Axis lines
      axis.line = element_line(color = "black", linewidth = 0.8),
      
      # X-axis text and title
      axis.text.x = element_text(size = text_size, color = "black", angle = Angle, hjust = 0.6),
      axis.title.x = element_text(size = text_size),
      
      # Y-axis text and title
      axis.text.y = element_text(size = text_size, color = "black"),
      axis.title.y = element_text(size = text_size),
      
      # Tick marks
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      
      # Legend
      legend.position = "none"
    )
  
  y_min <- floor(min(read_lengths_df$x) /
                   y_tick_interval) * y_tick_interval
  
  y_max <- ceiling(max(read_lengths_df$x) /
                     y_tick_interval) * y_tick_interval
  
  y_ticks = seq(y_min, y_max, by = y_tick_interval)
  
  violin_plot +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = expansion(mult = c(0, 0.02))
    )
}

###############################################################
# Generate violin plots of read-length distributions
###############################################################

make_violin_plot <- function(
    all_read_lengths_df,
    x_lab_all,
    y_lab_all,
    filtered_read_lengths_df,
    x_lab_filt,
    y_lab_filt,
    sample_names,
    Pool_sample_names,
    y_tick_interval,
    Title_filt,
    group_colors,
    output_file,
    output_width,
    output_height){
  
  text_size <- 16
  
  ###########################################################
  # Complete read-length distributions
  ###########################################################
  
  all_reads_pool_violin <- draw_violin_plot(
    read_lengths_df = all_read_lengths_df,
    group_col = "Pool_condition",
    sample_names = Pool_sample_names,
    y_tick_interval = y_tick_interval,
    Colors = group_colors,
    text_size = text_size,
    x_lab = x_lab_all,
    y_lab = y_lab_all,
    Title = ""
  )
  
  all_reads_sample_violin <- draw_violin_plot(
    read_lengths_df = all_read_lengths_df,
    group_col = "Condition",
    sample_names = sample_names,
    y_tick_interval = y_tick_interval,
    Colors = group_colors,
    text_size = text_size,
    x_lab = x_lab_all,
    y_lab = y_lab_all,
    Title = ""
  )
  
  ###########################################################
  # Filtered read-length distributions
  ###########################################################
  
  filtered_reads_pool_violin <- draw_violin_plot(
    read_lengths_df = filtered_read_lengths_df,
    group_col = "Pool_condition",
    sample_names = Pool_sample_names,
    y_tick_interval = y_tick_interval,
    Colors = group_colors,
    text_size = text_size,
    x_lab = x_lab_filt,
    y_lab = y_lab_filt,
    Title = Title_filt
  )
  
  filtered_reads_sample_violin <- draw_violin_plot(
    read_lengths_df = filtered_read_lengths_df,
    group_col = "Condition",
    sample_names = sample_names,
    y_tick_interval = y_tick_interval,
    Colors = group_colors,
    text_size = text_size,
    x_lab = x_lab_filt,
    y_lab = y_lab_filt,
    Title = ""
  )
  
  ###########################################################
  # Calculate relative panel widths
  ###########################################################
  
  pool_panel_width <- length(unique(filtered_read_lengths_df$Pool_condition))
  
  sample_panel_width <- length(unique(filtered_read_lengths_df$Condition))
  
  total_panel_width <- pool_panel_width + sample_panel_width
  
  pool_panel_width <- pool_panel_width / total_panel_width
  
  sample_panel_width <- sample_panel_width / total_panel_width
  
  ###########################################################
  # Assemble combined figure
  ###########################################################
  
  combined_plot <- ggdraw()
  
  combined_plot <- combined_plot +
    draw_plot(
      filtered_reads_pool_violin,
      x = 0,
      y = 0,
      width = pool_panel_width,
      height = 1
    )
  
  combined_plot <- combined_plot +
    draw_plot(
      filtered_reads_sample_violin,
      x = pool_panel_width,
      y = 0,
      width = sample_panel_width,
      height = 1
    )
  
  ###########################################################
  # Save figure
  ###########################################################
  
  ggsave(
    filename = output_file,
    plot = combined_plot,
    width = output_width,
    height = output_height,
    limitsize = FALSE
  )
  
  return(combined_plot)
}
#--------------------------------------------------
# 2D probability heatmap
#--------------------------------------------------
make_2d_prob_heatmap <- function(df, Breaks, Ranges,
                                 x_lab, y_lab,
                                 group_col,
                                 fig_name,
                                 fig_width,
                                 fig_height){

  colnames(df)[colnames(df) == group_col] <- 'group'

  df <- df %>%
    mutate(y_bin=cut(x, breaks=Breaks, labels=Ranges))

  df_prob <- df %>%
    count(group, y_bin) %>%
    group_by(group) %>%
    mutate(probability=n/sum(n)) %>% ungroup()

  blues_palette <- colorRampPalette(brewer.pal(8,'Blues'))(25)

  text_size = 16
  Angle = 0 # 45
  
  fig <- ggplot(df_prob, aes(group, y_bin, fill=probability)) +
    geom_tile() +
    scale_fill_gradientn(colors = blues_palette, name = "Probability") +  # Use probabilities for color scale
    labs(title = "2D Probability Distribution Heatmap", x = x_lab, y = y_lab) +
    theme_minimal() +
    theme(
      text = element_text(size = text_size),
      
      # Remove grid lines
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, size = text_size),  # Left-align title
      axis.text.x = element_text(size = text_size, angle = Angle, hjust = 0.6),
      axis.title.x = element_text(size = text_size),              # X-axis label size
      axis.text.y = element_text(size = text_size, color = "black"),  # Show x-axis labels
      axis.title.y = element_text(size = text_size)            # Y-axis label size
      
    )

  ggsave(fig_name, fig, width=fig_width, height=fig_height, limitsize = FALSE)

  fig
}

#--------------------------------------------------
# Main plotting function
#--------------------------------------------------
plot_read_length <- function(RL_paths,
                             Results_path,
                             sample_names,
                             Pool_sample_names,
                             read_min,
                             Read_max_lim,
                             y_tick_interval,
                             Colors,
                             violin_file_name,
                             violin_sz,
                             filename_2d_prob,
                             sz_2d_prob){

  all_read_lengths_df <- data.frame()
  filtered_read_lengths_df <- data.frame()
  
  for (sample_idx in seq_along(RL_paths)) {
    
    read_length_counts <- read.table(RL_paths[sample_idx])
    
    # Reconstruct individual read lengths from count table
    read_lengths <- unlist(
      mapply(
        function(count, read_length) replicate(count, read_length),
        read_length_counts$V1,
        read_length_counts$V2
      )
    )
    
    read_lengths <- read_lengths[read_lengths > read_min]
    
    all_read_lengths_df <- rbind(
      all_read_lengths_df,
      data.frame(
        x = read_lengths,
        Condition = sample_names[sample_idx],
        Pool_condition = Pool_sample_names[sample_idx]
      )
    )
    
    if (Read_max_lim == "use_99_quantile") {
      filtered_read_lengths <- read_lengths[
        read_lengths <= quantile(read_lengths, 0.99)
      ]
    } else {
      filtered_read_lengths <- read_lengths[
        read_lengths <= Read_max_lim
      ]
    }
    
    filtered_read_lengths_df <- rbind(
      filtered_read_lengths_df,
      data.frame(
        x = filtered_read_lengths,
        Condition = sample_names[sample_idx],
        Pool_condition = Pool_sample_names[sample_idx]
      )
    )
  }

  dir.create(Results_path, recursive=TRUE, showWarnings=FALSE)

  violin_plot <- make_violin_plot(
    all_read_lengths_df = all_read_lengths_df,
    x_lab_all = "Samples",
    y_lab_all = "Read length (bases)",
    
    filtered_read_lengths_df = filtered_read_lengths_df,
    x_lab_filt = "Samples",
    y_lab_filt = "Read length (bases)",
    
    sample_names = sample_names,
    Pool_sample_names = Pool_sample_names,
    
    y_tick_interval = y_tick_interval,
    Title_filt = "",
    
    group_colors = Colors,
    
    output_file = file.path(Results_path, violin_file_name),
    
    output_width = violin_sz[1],
    output_height = violin_sz[2]
  )
  
  length_distribution_heatmap <- make_2d_prob_heatmap(
    all_read_lengths_df,
    Breaks = c(0, 200, 300, 400, 500, 1000, Inf),
    Ranges = c("<200", "200-300", "300-400", "400-500", "500-1000", ">1000"),
    x_lab = "Samples",
    y_lab = "Read length ranges",
    group_col = "Pool_condition",
    fig_name = file.path(Results_path, filename_2d_prob),
    fig_width = sz_2d_prob[1],
    fig_height = sz_2d_prob[2]
  )

  return(list(violin=violin_plot, heatmap=length_distribution_heatmap))
}

######################################################################
# Project directory
######################################################################

Project_path <- "/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/"
setwd(Project_path)

######################################################################
# Direct RNA SRSF3 samples
######################################################################

RL_paths <- c(
  "Result/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/read_length.txt",
  "Result/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/read_length.txt",
  "Result/3_Genome_Mapping/SRSF3_dcDNA/rep3_1pct/read_length.txt",
  "Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/read_length.txt",
  "Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/read_length.txt",
  "Result/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep3_0.2pct/read_length.txt"
)

Conditions <- c(
  "CLIP rep1",
  "CLIP rep2",
  "CLIP rep3",
  "-UV rep1",
  "-UV rep2",
  "-UV rep3"
)

pooling_Conditions <- c(
  "CLIP",
  "CLIP",
  "CLIP",
  "-UV",
  "-UV",
  "-UV"
)

Results_path <- "Result/figure/6_Read_length/"

Colors <- c(
  "CLIP" = "#66c2a5",
  "-UV" = "#08519C"
)

######################################################################
# Plot parameters
######################################################################

read_min <- 20
y_tick_interval <- 1000

Read_max_lim <- "use_99_quantile"

######################################################################
# Generate figures
######################################################################

f10 <- plot_read_length(
  RL_paths = RL_paths,
  Results_path = Results_path,
  sample_names = Conditions,
  Pool_sample_names = pooling_Conditions,
  Colors = Colors,
  read_min = read_min,
  Read_max_lim = Read_max_lim,
  y_tick_interval = y_tick_interval,
  violin_file_name = "RL_violin_from_bam.pdf",
  violin_sz = c(10, 7),
  filename_2d_prob = "RL_2d_prob_heatmap_from_bam.pdf",
  sz_2d_prob = c(7, 7)
)