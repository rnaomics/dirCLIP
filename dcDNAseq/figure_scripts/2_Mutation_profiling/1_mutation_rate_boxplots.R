#!/usr/bin/env Rscript

# Figure: Mutation Rate Distribution Plots
# Author: Mohamed Bahrudeen
# Purpose: Generate mutation-rate boxplots used in manuscript figures.

# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(scales)
library(S4Vectors)
library(matrixStats)
library(GenomicRanges)
library(tibble)

############ Function to get box plot values ###################

get_boxplot_values = function(mutation_rates){
  
  q1 <- quantile(mutation_rates, 0.25)   # First quartile
  median_value <- quantile(mutation_rates, 0.5)
  q3 <- quantile(mutation_rates, 0.75)   # Third quartile
  
  IQR_value <- q3 - q1
  
  # Lower and Upper Whiskers based on 1.5 × IQR rule
  lower_whisker <- min(mutation_rates[mutation_rates >= (q1 - 1.5 * IQR_value)])
  upper_whisker <- max(mutation_rates[mutation_rates <= (q3 + 1.5 * IQR_value)])
  
  # Number of data points
  sample_count <- length(mutation_rates)
  
  box_plot_values <- data.frame(lower_whisker, q1, median_value, q3, upper_whisker, sample_count)
  
  return(box_plot_values)
}

############### Function to get the boxplot values for coverage thresholds ################################

get_boxplot_values_for_Coverage_thresholds = function(
    Paths,
    Samples,
    isdel,
    Mut_file,
    Coverage_thresholds
){
  
  box_plot_values_list = list()
  
  sample_idx = 1
  
  for (Path in Paths){
    
    mutation_rate_df = read.delim(paste0(Path, Mut_file), header = FALSE)
    
    colnames(mutation_rate_df) = c("Chr", "loc", "Coverage", "Mut_rate")
    
    #### If Common_locs length is greater than 0, then filter the common locations ####
    threshold_idx = 1
    
    for (cov_thresh in Coverage_thresholds){
      
      #### Mutation rates above a given coverage threshold ####
      
      mutation_rates = mutation_rate_df[mutation_rate_df$Coverage >= cov_thresh,]$Mut_rate
      
      if (threshold_idx == 1){
        
        box_plot_values = get_boxplot_values(mutation_rates)
        
      } else {
        
        box_plot_values = rbind(box_plot_values, get_boxplot_values(mutation_rates))
      }
      
      threshold_idx = threshold_idx + 1
    }
    
    rownames(box_plot_values) = Coverage_thresholds
    
    box_plot_values_list[[sample_idx]] = box_plot_values
    
    sample_idx = sample_idx + 1
  }
  
  box_plot_values_tibble = tibble(
    box_plot_values_list,
    Samples,
    isdel
  )
  
  return(box_plot_values_tibble)
}

############### Function to make_box_plot_type1 ################################
make_box_plot_type1 = function(
    boxplot_stats,
    Colors,
    x_lab,
    y_lab,
    Title,
    y_limit,
    text_size,
    text_size2
){
  
  boxplot_figure = ggplot(boxplot_stats, aes(x = Samples, fill = Samples)) +
    
    geom_boxplot(
      aes(
        ymin = lower_whisker, lower = q1, middle = median_value, upper = q3, ymax = upper_whisker),
      stat = "identity",
      width = 0.6,
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    ) +
    
    scale_fill_manual(values = Colors) +
    
    geom_text(
      data = boxplot_stats,
      aes(x = Samples, y = max(upper_whisker) * 1.05, label = paste0("n=", sample_count)),
      vjust = 0,
      size = text_size2,
      fontface = "bold",
      color = "black"
    ) +
    
    labs(
      y = y_lab,
      x = x_lab,
      fill = x_lab,
      title = Title
    )
  
  boxplot_figure = add_theme(boxplot_figure, text_size)
  
  return(boxplot_figure)
}

############### Function to make_box_plot_type2 ################################
make_box_plot_type2 = function(
    boxplot_stats,
    Colors,
    x_lab,
    y_lab,
    Title,
    legend_tile,
    y_limit,
    text_size,
    text_size2
){
  
  boxplot_figure = ggplot(boxplot_stats, aes(x = is_del, fill = Samples)) +
    
    geom_boxplot(
      aes(ymin = lower_whisker, lower = q1, middle = median_value, upper = q3, ymax = upper_whisker),
      stat = "identity",
      width = 0.6,
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    ) +
    
    scale_fill_manual(values = Colors) +
    
    labs(
      y = y_lab,
      x = x_lab,
      fill = legend_tile,
      title = Title
    )
  
  boxplot_figure = add_theme(boxplot_figure, text_size)
  
  return(boxplot_figure)
}

############### Function to add_theme ################################
add_theme = function(plot_object, text_size){
  
  plot_object = plot_object + theme_minimal(text_size) +
    
    theme(
      text = element_text(size = text_size),
      
      panel.grid = element_blank(),
      
      panel.border = element_blank(),
      
      axis.line = element_line(color = "black", linewidth = 0.8),
      
      axis.text.x = element_text(size = text_size, color = "black", angle = 0, hjust = 1),
      
      axis.title.x = element_text(size = text_size),
      
      axis.text.y = element_text(size = text_size, color = "black"),
      
      axis.title.y = element_text(size = text_size),
      
      legend.title = element_text(size = text_size),
      
      legend.text = element_text(size = text_size),
      
      axis.ticks = element_line(color = "black", size = 0.8),
      
      axis.ticks.length = unit(0.25, "cm"),
      
      axis.ticks.x = element_line(color = "black"),
      
      axis.ticks.y = element_line(color = "black")
    )
  
  return(plot_object)
}

################# Function to make box plot #####################################

make_box_plot = function(
    box_plot_values_list,
    cov_thresh,
    Samples,
    isdel,
    Sample_order,
    Sample_Colors,
    del_Colors,
    x_lab,
    y_lab,
    Title,
    y_lim,
    y_tick_interval
){
  
  # isdel = "No"
  
  selected_rows = box_plot_values_list$isdel %in% isdel
  
  box_plot_values_list = box_plot_values_list[selected_rows,]
  
  selected_samples = box_plot_values_list$Samples
  
  is_del = box_plot_values_list$isdel
  
  for (sample_idx in 1:length(selected_samples)) {
    
    current_sample = selected_samples[sample_idx]
    
    print(current_sample)
    
    ###### Get the box_plot values ##################
    
    if (sample_idx == 1){
      
      boxplot_stats = box_plot_values_list[sample_idx,]$box_plot_values_list[[1]][as.character(cov_thresh),]
      
    } else {
      
      boxplot_stats = rbind(boxplot_stats, box_plot_values_list[sample_idx,]$box_plot_values_list[[1]][as.character(cov_thresh),])
    }
  }
  
  boxplot_stats$Samples = factor(selected_samples, levels = Sample_order)
  
  boxplot_stats$is_del = is_del
  
  # Plot using ggplot2
  
  text_size = 16
  text_size2 = 6
  
  if (length(y_lim) == 0){
    y_lim = c(0, max(boxplot_stats$upper_whisker) * 1.1)
  }
  
  ###########################################################
  
  if (length(unique(is_del)) == 1){
    
    boxplot_figure =
      make_box_plot_type1(
        boxplot_stats,
        Sample_Colors,
        x_lab,
        y_lab,
        Title,
        y_lim,
        text_size,
        text_size2
      )
    
  } else {
    
    boxplot_figure =
      make_box_plot_type2(
        boxplot_stats,
        Sample_Colors,
        x_lab,
        y_lab,
        Title,
        legend_tile = "Sample",
        y_lim,
        text_size,
        text_size2
      )
  }
  
  ########## setting the ticks ############################
  
  y_min = floor(y_lim[1] / y_tick_interval) * y_tick_interval
  
  y_max = ceiling(y_lim[2] / y_tick_interval) * y_tick_interval
  
  y_ticks = seq(y_min, y_max, by = y_tick_interval)
  
  boxplot_figure = boxplot_figure +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = expansion(mult = c(0, 0.02))
    )
  
  return(boxplot_figure)
}

############### Merging the box plots ########################################

make_joined_box_plots = function(
    box_plot_values_list,
    Coverage_thresholds,
    Samples,
    isdel,
    Sample_order,
    Sample_Colors,
    del_Colors,
    x_lab,
    y_lab,
    y_lim,
    y_tick_interval,
    fig_name,
    fig_sz,
    sub_plots
){
  
  combined_figure = ggdraw()
  
  sub_fig_sz = 1 / sub_plots
  
  combinations <- expand.grid(x = 0:(sub_plots[1] - 1), y = 0:(sub_plots[2] - 1))
  
  plot_idx = 1
  
  for (cov_thresh in Coverage_thresholds){
    
    print(cov_thresh)
    
    subplot_figure =
      make_box_plot(
        box_plot_values_list,
        cov_thresh,
        Samples = Samples,
        isdel,
        Sample_order,
        Sample_Colors,
        del_Colors,
        x_lab = x_lab,
        y_lab = y_lab,
        Title = paste0("Cov thresh = ", cov_thresh),
        y_lim = y_lim,
        y_tick_interval
      )
    
    combined_figure = combined_figure +
      draw_plot(
        subplot_figure,
        x = combinations$x[plot_idx] * sub_fig_sz[1],
        y = combinations$y[plot_idx] * sub_fig_sz[2],
        width = sub_fig_sz[1],
        height = sub_fig_sz[2]
      )
    
    plot_idx = plot_idx + 1
  }
  
  ####### saving the plot #################################################################
  
  ggsave(
    fig_name,
    plot = combined_figure,
    width = fig_sz[1] * sub_plots[1],
    height = fig_sz[2] * sub_plots[2],
    limitsize = FALSE,
    device = cairo_pdf
  )
}
########################################################################################

Project_path <- "Tools/"
setwd(Project_path)

# Sample definitions
mutation_rate_paths <- c(
  "Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/q_8_Q_8/1_No_del/",
  "Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/q_8_Q_8/1_No_del/",
  
  "Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/q_8_Q_8/1_No_del/",
  "Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/q_8_Q_8/1_No_del/",
  
  "Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep1_1pct/q_8_Q_8/2_Inc_del/",
  "Result/5_Mutation_profiling/3_Genome_Mapping/SRSF3_dcDNA/rep2_1pct/q_8_Q_8/2_Inc_del/",
  
  "Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep1_0.2pct/q_8_Q_8/2_Inc_del/",
  "Result/5_Mutation_profiling/3_Genome_Mapping/TotalRNA_dcDNA_no_UV/rep2_0.2pct/q_8_Q_8/2_Inc_del/"
)

sample_names = c(
  "SRSF3 dcDNA 1",
  "SRSF3 dcDNA 2",
  "noUV dcDNA 1",
  "noUV dcDNA 2",
  "SRSF3 dcDNA 1",
  "SRSF3 dcDNA 2",
  "noUV dcDNA 1",
  "noUV dcDNA 2"
)

sample_order = c(
  "SRSF3 dcDNA 1",
  "SRSF3 dcDNA 2",
  "noUV dcDNA 1",
  "noUV dcDNA 2"
)

include_deletions = c(
  "No", "No", "No", "No",
  "Yes", "Yes", "Yes", "Yes"
)

mutation_rate_file = "Mutation_rate_cov_10_mut_0.tsv.gz"

######## Colors
sample_colors = c(
  "SRSF3 dcDNA 1" = "#66c2a5",
  "SRSF3 dcDNA 2" = "#66c2a5",
  
  "noUV dcDNA 1" = "#08519C",
  "noUV dcDNA 2" = "#08519C"
)

deletion_colors = c(
  "Yes" = "gray40",
  "No" = "gray70"
)

########## Output settings
result_path_no_del = "Result/figure/1_Mutation_profiling/"

result_path_with_del = "Result/figure/1_Mutation_profiling/"

single_plot_size = c(6, 7)
comparison_plot_size = c(12, 7)

coverage_thresholds = c(50)

######### Generate boxplot statistics
boxplot_values_list =
  get_boxplot_values_for_Coverage_thresholds(
    mutation_rate_paths,
    sample_names,
    include_deletions,
    mutation_rate_file,
    coverage_thresholds
  )

######### 1. Plot generation: Without deletions ######
y_lim = c()
y_tick_interval = 0.01

make_joined_box_plots(
  boxplot_values_list,
  coverage_thresholds,
  Samples = sample_names,
  isdel = "No",
  sample_order,
  sample_colors,
  deletion_colors,
  x_lab = "Sample",
  y_lab = "Mut rate",
  y_lim = y_lim,
  y_tick_interval = y_tick_interval,
  fig_name = paste0(result_path_no_del, "Mut_Rate_nodel.pdf"),
  fig_sz = single_plot_size,
  sub_plots = c(1,1)
)

######### 2. Plot generation: With deletions ######
y_lim = c()
y_tick_interval = 0.025

make_joined_box_plots(
  boxplot_values_list,
  coverage_thresholds,
  Samples = sample_names,
  isdel = "Yes",
  sample_order,
  sample_colors,
  deletion_colors,
  x_lab = "Sample",
  y_lab = "Mut rate",
  y_lim = y_lim,
  y_tick_interval = y_tick_interval,
  fig_name = paste0(result_path_no_del, "Mut_Rate_del.pdf"),
  fig_sz = single_plot_size,
  sub_plots = c(1,1)
)

######### 3. Plot generation: With/without deletions ######
y_lim = c()
y_tick_interval = 0.025

make_joined_box_plots(
  boxplot_values_list,
  coverage_thresholds,
  Samples = sample_names,
  isdel = c("No", "Yes"),
  sample_order,
  sample_colors,
  deletion_colors,
  x_lab = "Sample",
  y_lab = "Mut rate",
  y_lim = y_lim,
  y_tick_interval = y_tick_interval,
  fig_name = paste0(result_path_no_del, "Mut_Rate_del_nodel.pdf"),
  fig_sz = single_plot_size,
  sub_plots = c(1,1)
)

############################################