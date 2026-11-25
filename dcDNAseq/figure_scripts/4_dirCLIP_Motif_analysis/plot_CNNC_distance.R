#!/usr/bin/env Rscript
# CNNC Distance Plotting

# Generates frequency and cumulative-frequency plots of CNNC motif distances relative to mutation sites.
# Preserve original function implementations from the manuscript script.

# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

## Libraries ######
library(ggplot2)
library(tidyverse)
library(writexl)
library(dplyr)
library("cowplot")

# ---------------------------------------------------------------------------------------------------
# 1. Function to plot CNNC distance versus cumulative fraction
# ---------------------------------------------------------------------------------------------------

plot_CNNC_dist_vs_cum_frac = function(
    type,
    motif_path,
    CNNC_dist_file,
    CNNC_dist_randomers_file,
    do_abs_cum_freq,
    do_diff,
    Title,
    y_lim,
    y_tick_interval,
    x_tick_interval
){
  
  CNNC_dist <- read.delim(
    paste0(
      motif_path,
      CNNC_dist_file
    )
  )
  
  colnames(CNNC_dist)[-1] =
    "Actual_seq"
  
  CNNC_dist_randomers <- read.delim(
    paste0(
      motif_path,
      CNNC_dist_randomers_file
    )
  )
  
  colnames(CNNC_dist_randomers)[-1] =
    paste0(
      type,
      "_",
      1:(ncol(CNNC_dist_randomers) - 1)
    )
  
  CNNC_dist =
    merge(
      CNNC_dist,
      CNNC_dist_randomers,
      by = "CNNC_locs"
    )
  
  #################################################
  # Normalize frequencies
  #################################################
  
  column_sums <- colSums(
    CNNC_dist[-1]
  )
  
  CNNC_dist <- CNNC_dist %>%
    mutate(
      across(
        -CNNC_locs,
        ~ . / column_sums[
          as.character(
            cur_column()
          )
        ]
      )
    )
  
  #################################################
  # Cumulative frequencies
  #################################################
  
  if (do_abs_cum_freq){
    
    CNNC_dist$CNNC_locs =
      abs(CNNC_dist$CNNC_locs)
    
    CNNC_dist <- CNNC_dist %>%
      group_by(CNNC_locs) %>%
      summarise(
        across(
          everything(),
          sum,
          na.rm = TRUE
        )
      )
    
    CNNC_dist <- CNNC_dist %>%
      mutate(
        across(
          -CNNC_locs,
          cumsum
        )
      )
    
    if (do_diff){
      
      CNNC_dist <- CNNC_dist %>%
        mutate(
          across(
            -CNNC_locs,
            ~ Actual_seq - .
          )
        )
      
      y_lab =
        "Diff Cum Freq (Act-Shuf)"
      
    } else {
      
      y_lab =
        "Cum Freq"
    }
    
  } else {
    
    y_lab =
      "Freq"
  }
  
  #################################################
  # Mean and SD of random controls
  #################################################
  
  CNNC_dist$random_mean =
    rowMeans(
      CNNC_dist[
        ,
        paste0(
          type,
          "_",
          1:(ncol(CNNC_dist_randomers) - 1)
        )
      ]
    )
  
  CNNC_dist$random_sd =
    apply(
      CNNC_dist[
        ,
        paste0(
          type,
          "_",
          1:(ncol(CNNC_dist_randomers) - 1)
        )
      ],
      1,
      sd
    )
  
  CNNC_dist =
    CNNC_dist[
      ,
      c(
        "CNNC_locs",
        "Actual_seq",
        "random_mean",
        "random_sd"
      )
    ]
  
  #################################################
  # Remove infinite distances
  #################################################
  
  selected_idx =
    CNNC_dist$CNNC_locs == Inf
  
  CNNC_dist =
    CNNC_dist[
      !selected_idx,
    ]
  
  #################################################
  # Plot settings
  #################################################
  
  line_colors = c(
    "Actual_seq"  = "#2ca25f",
    "random_mean" = "#99d8c9"
  )
  
  error_colors = c(
    "random_sd" = "gray60"
  )
  
  alpha_value = 0.8
  
  text_size = 16
  
  #################################################
  # Create plot
  #################################################
  
  cnnc_plot =
    ggplot(
      CNNC_dist,
      aes(x = CNNC_locs)
    ) +
    
    geom_ribbon(
      aes(
        ymin = random_mean - random_sd,
        ymax = random_mean + random_sd,
        fill = "random_sd"
      ),
      alpha = alpha_value,
      color = NA
    ) +
    
    geom_line(
      aes(
        y = Actual_seq,
        color = "Actual_seq"
      ),
      alpha = alpha_value,
      size = 1
    ) +
    
    geom_line(
      aes(
        y = random_mean,
        color = "random_mean"
      ),
      alpha = alpha_value,
      size = 1
    ) +
    
    geom_point(
      aes(
        y = Actual_seq,
        color = "Actual_seq"
      ),
      alpha = alpha_value,
      size = 1.5
    ) +
    
    geom_point(
      aes(
        y = random_mean,
        color = "random_mean"
      ),
      alpha = alpha_value,
      size = 1.5
    ) +
    
    scale_color_manual(
      name = "Line",
      values = line_colors
    ) +
    
    scale_fill_manual(
      name = "Error",
      values = error_colors
    )
  
  cnnc_plot = cnnc_plot +
    labs(
      x = "CNNC distance from Mut sites",
      y = y_lab,
      title = Title
    ) +
    theme_minimal(base_size = text_size) +
    theme(
      text = element_text(size = text_size),
      
      # Remove grid lines
      panel.grid = element_blank(),
      
      # Add border only on left and bottom
      panel.border = element_blank(),
      
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      
      axis.text = element_text(
        size = text_size,
        hjust = 1
      ),
      
      axis.title = element_text(
        size = text_size
      ),
      
      axis.text.x = element_text(
        hjust = 0.6
      ),
      
      legend.title = element_text(
        size = text_size
      ),
      
      legend.text = element_text(
        size = text_size
      ),
      
      # draw tick lines
      axis.ticks = element_line(
        color = "black",
        size = 0.8
      ),
      
      # tick mark length
      axis.ticks.length = unit(
        0.25,
        "cm"
      ),
      
      axis.ticks.x = element_line(
        color = "black"
      ),
      
      axis.ticks.y = element_line(
        color = "black"
      )
    )
  
  #################################################
  # Y-axis
  #################################################
  
  if (length(y_lim) > 0){
    
    y_min = y_lim[1]
    y_max = y_lim[2]
    
  } else {
    
    y_min =
      floor(
        min(
          CNNC_dist[
            ,
            c(
              "Actual_seq",
              "random_mean"
            )
          ]
        ) /
          y_tick_interval
      ) *
      y_tick_interval
    
    y_max =
      ceiling(
        max(
          CNNC_dist[
            ,
            c(
              "Actual_seq",
              "random_mean"
            )
          ]
        ) /
          y_tick_interval
      ) *
      y_tick_interval
  }
  
  y_ticks =
    seq(
      y_min,
      y_max,
      by = y_tick_interval
    )
  
  cnnc_plot =
    cnnc_plot +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = c(0, 0)
    )
  
  #################################################
  # X-axis
  #################################################
  
  x_min =
    floor(
      min(CNNC_dist$CNNC_locs) /
        x_tick_interval
    ) *
    x_tick_interval
  
  x_max =
    ceiling(
      max(CNNC_dist$CNNC_locs) /
        x_tick_interval
    ) *
    x_tick_interval
  
  x_ticks =
    seq(
      x_min,
      x_max,
      by = x_tick_interval
    )
  
  cnnc_plot =
    cnnc_plot +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = x_ticks,
      expand = c(0, 0)
    )
  
  return(cnnc_plot)
}

# ---------------------------------------------------------------------------------------------------
# 2. Function to make a merged CNNC-distance plot
# ---------------------------------------------------------------------------------------------------

plot_CNNC_dist_vs_cum_frac_joined_plots = function(
    Path,
    Coverage_thresholds,
    Diff_thresholds,
    window_length,
    mut_type,
    kmer,
    do_abs_cum_freq,
    do_diff,
    y_lim,
    y_tick_interval,
    x_tick_interval,
    fig_name,
    fig_sz,
    sub_plots
){
  
  combined_plot = ggdraw()
  
  subplot_size = 1 / sub_plots
  
  combinations <- expand.grid(
    x = 0:(sub_plots[1] - 1),
    y = 0:(sub_plots[2] - 1)
  )
  
  plot_idx = 1
  
  motif_path = paste0(
    Path,
    "cov_thresh_",
    Coverage_thresholds,
    "/",
    mut_type, "/win_", window_length,
    "/kmer_",
    kmer,
    "/"
  )
  
  CNNC_dist_file = paste0(
    "CNNC_dist_thresh_",
    Diff_thresholds,
    ".tsv"
  )
  
  CNNC_dist_randomers_file = paste0(
    "CNNC_dist_randomers_thresh_",
    Diff_thresholds,
    ".tsv"
  )
  
  CNNC_dist_rnd_bg_file = paste0(
    "CNNC_dist_rnd_bg_thresh_",
    Diff_thresholds,
    ".tsv"
  )
  
  Title = paste0(
    "Cov_thresh = ",
    Coverage_thresholds,
    ", Diff_thresh = ",
    Diff_thresholds
  )
  
  ###########################################################
  # Shuffled-sequence control
  ###########################################################
  
  shuffled_control_plot =
    plot_CNNC_dist_vs_cum_frac(
      type = "shuffling",
      motif_path,
      CNNC_dist_file,
      CNNC_dist_randomers_file,
      do_abs_cum_freq,
      do_diff,
      Title = paste0(
        "shuffling: ",
        Title
      ),
      y_lim,
      y_tick_interval,
      x_tick_interval
    )
  
  ###########################################################
  # Random-background control
  ###########################################################
  
  random_background_plot =
    plot_CNNC_dist_vs_cum_frac(
      type = "rnd_bg",
      motif_path,
      CNNC_dist_file,
      CNNC_dist_rnd_bg_file,
      do_abs_cum_freq,
      do_diff,
      Title = paste0(
        "rnd_bg: ",
        Title
      ),
      y_lim,
      y_tick_interval,
      x_tick_interval
    )
  
  ###########################################################
  # Combine plots
  ###########################################################
  
  combined_plot =
    combined_plot +
    draw_plot(
      shuffled_control_plot,
      x = 0,
      y = 0,
      width = 0.5,
      height = 1
    )
  
  combined_plot =
    combined_plot +
    draw_plot(
      random_background_plot,
      x = 0.5,
      y = 0,
      width = 0.5,
      height = 1
    )
  
  ###########################################################
  # Save figure
  ###########################################################
  
  ggsave(
    fig_name,
    plot = combined_plot,
    width = fig_sz[1] * 2,
    height = fig_sz[2] * 1,
    limitsize = FALSE,
    device = cairo_pdf
  )
}

# ---------------------------------------------------------------------------------------------------
# 3. Function to Create folder
# ---------------------------------------------------------------------------------------------------

create_folder = function(Path){
  ifelse(!dir.exists(Path), dir.create(Path, recursive = TRUE), FALSE)
}


################################ Paths ##########################################################################################################

motif_analysis_paths = c(
  "Result/figure/2_Motif_analysis/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/"
)

###################################################################################################################################

################## Setting the Coverage thresholds and difference thresholds ###############

Coverage_thresholds = c(50)
Diff_thresholds = c(0.03)

sub_plots = c(1,1)

############################################################################################

kmer = 5

window_length = 100

################### Choose sites or region centre #############################################

mut_types = c(
  "sites"
)

##############################################################################################

for (mut_type in mut_types){
  
  for (path_idx in seq_along(motif_analysis_paths)){
    
    motif_path = motif_analysis_paths[path_idx]
    
    print(motif_path)
    
    output_path = paste0(
      motif_path,
      "plot_CNNC_dist/",
      mut_type,
      "/win_", window_length, "/"
    )
    
    create_folder(output_path)
    
    ####################################################
    # 1. CNNC distance vs cumulative frequency
    ####################################################
    
    y_tick_interval = 0.15
    x_tick_interval = 10
    
    plot_CNNC_dist_vs_cum_frac_joined_plots(
      Path = motif_path,
      Coverage_thresholds,
      Diff_thresholds,
      window_length,
      mut_type,
      kmer,
      do_abs_cum_freq = TRUE,
      do_diff = FALSE,
      y_lim = c(),
      y_tick_interval,
      x_tick_interval,
      fig_name = paste0(
        output_path,
        "CNNC_dist_cumfreq_L_",
        window_length,
        "_sites.pdf"
      ),
      fig_sz = c(7, 7),
      sub_plots = sub_plots
    )
    
    ####################################################
    # 2. CNNC distance vs frequency
    ####################################################
    
    y_tick_interval = 0.005
    x_tick_interval = 25
    
    plot_CNNC_dist_vs_cum_frac_joined_plots(
      Path = motif_path,
      Coverage_thresholds,
      Diff_thresholds,
      window_length,
      mut_type,
      kmer,
      do_abs_cum_freq = FALSE,
      do_diff = FALSE,
      y_lim = c(),
      y_tick_interval,
      x_tick_interval,
      fig_name = paste0(
        output_path,
        "CNNC_dist_freq_L_",
        window_length,
        "_sites.pdf"
      ),
      fig_sz = c(7, 7),
      sub_plots = sub_plots
    )
  }
}