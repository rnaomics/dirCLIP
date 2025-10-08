# Optional alternative library path
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

# Core libraries for sequence and motif visualization
library(ShortRead)
library(stringr)
library(stringi)
library(ggplot2)
library(ggseqlogo)      # Installed via: devtools::install_github("omarwagih/ggseqlogo")
library(Biostrings)
library(tidyverse)
library(writexl)
library(dplyr)
library(cowplot)

##### 🎯 Function: Plot CNNC Motif Distance vs. Cumulative Frequency
plot_CNNC_dist_vs_cum_frac <- function(motif_path, CNNC_dist_file, CNNC_dist_randomers_file, 
                                       do_abs_cum_freq, do_diff, Title, y_lim, y_tick_interval, x_tick_interval) {
  
  ### 1. 🗂️ Load CNNC motif position data
  type <- "shuffling"
  CNNC_dist <- read.delim(file.path(motif_path, CNNC_dist_file))
  colnames(CNNC_dist)[-1] <- "Actual_seq"
  
  CNNC_dist_randomers <- read.delim(file.path(motif_path, CNNC_dist_randomers_file))
  colnames(CNNC_dist_randomers)[-1] <- paste0(type, "_", seq_len(ncol(CNNC_dist_randomers) - 1))
  
  CNNC_dist <- merge(CNNC_dist, CNNC_dist_randomers, by = "CNNC_locs")
  
  ### 2. 📊 Normalize motif counts
  col_sums <- colSums(CNNC_dist[,-1])
  CNNC_dist <- CNNC_dist %>%
    mutate(across(-CNNC_locs, ~ . / col_sums[cur_column()]))
  
  ### 3. 🔢 Compute absolute positions and cumulative frequencies if requested
  if (do_abs_cum_freq) {
    CNNC_dist$CNNC_locs <- abs(CNNC_dist$CNNC_locs)
    
    CNNC_dist <- CNNC_dist %>%
      group_by(CNNC_locs) %>%
      summarise(across(everything(), sum, na.rm = TRUE)) %>%
      mutate(across(-CNNC_locs, cumsum))
    
    if (do_diff) {
      CNNC_dist <- CNNC_dist %>% mutate(across(-CNNC_locs, ~ Actual_seq - .))
      y_lab <- "Diff Cum Freq (Act − Shuf)"
    } else {
      y_lab <- "Cum Freq"
    }
    
  } else {
    y_lab <- "Freq"
  }
  
  ### 4. 📉 Compute mean and standard deviation of shuffled motif distributions
  CNNC_dist$random_mean <- rowMeans(select(CNNC_dist, starts_with(type)))
  CNNC_dist$random_sd   <- apply(select(CNNC_dist, starts_with(type)), 1, sd)
  
  CNNC_dist <- CNNC_dist[, c("CNNC_locs", "Actual_seq", "random_mean", "random_sd")]
  CNNC_dist <- CNNC_dist[is.finite(CNNC_dist$CNNC_locs), ]  # Remove Inf positions
  
  ### 5. 🎨 Set plotting aesthetics
  line_colors <- c("Actual_seq" = "#2ca25f", "random_mean" = "#99d8c9")
  error_colors <- c("random_sd" = "gray60")
  text_size <- 16
  alpha <- 0.8
  
  ### 6. 📈 Build cumulative motif plot
  fig <- ggplot(CNNC_dist, aes(x = CNNC_locs)) +
    geom_ribbon(aes(ymin = random_mean - random_sd,
                    ymax = random_mean + random_sd,
                    fill = "random_sd"), alpha = alpha, color = NA) +
    geom_line(aes(y = Actual_seq, color = "Actual_seq"), alpha = alpha, size = 1) +
    geom_line(aes(y = random_mean, color = "random_mean"), alpha = alpha, size = 1) +
    geom_point(aes(y = Actual_seq, color = "Actual_seq"), alpha = alpha, size = 1.5) +
    geom_point(aes(y = random_mean, color = "random_mean"), alpha = alpha, size = 1.5) +
    scale_color_manual(name = "Line", values = line_colors) +
    scale_fill_manual(name = "Error", values = error_colors) +
    labs(title = Title, x = "CNNC distance from mutation sites", y = y_lab) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.text = element_text(size = text_size, color = "black"),
      axis.text.x = element_text(hjust = 0.6),
      axis.title = element_text(size = text_size),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_text(size = text_size),
      legend.text  = element_text(size = text_size)
    )
  
  ### 7. 🧮 Define y-axis breaks
  if (length(y_lim) > 0) {
    y_min <- y_lim[1]
    y_max <- y_lim[2]
  } else {
    y_min <- floor(min(CNNC_dist[, c("Actual_seq", "random_mean")]) / y_tick_interval) * y_tick_interval
    y_max <- ceiling(max(CNNC_dist[, c("Actual_seq", "random_mean")]) / y_tick_interval) * y_tick_interval
  }
  
  fig <- fig + scale_y_continuous(limits = c(y_min, y_max),
                                  breaks = seq(y_min, y_max, y_tick_interval),
                                  expand = c(0, 0))
  
  ### 8. 🧮 Define x-axis breaks
  x_min <- floor(min(CNNC_dist$CNNC_locs) / x_tick_interval) * x_tick_interval
  x_max <- ceiling(max(CNNC_dist$CNNC_locs) / x_tick_interval) * x_tick_interval
  
  fig <- fig + scale_x_continuous(limits = c(x_min, x_max),
                                  breaks = seq(x_min, x_max, x_tick_interval),
                                  expand = c(0, 0))
  
  return(fig)
}


##### 🎨 Function: Generate and Merge CNNC Distribution Subplots
plot_CNNC_dist_vs_cum_frac_joined_plots <- function(Path, CNNC_motif_group, Coverage_thresholds, Diff_thresholds, 
                                                    L, mut_type, kmer, do_abs_cum_freq, do_diff, 
                                                    y_lim, y_tick_interval, x_tick_interval, 
                                                    fig_name, fig_sz, sub_plots) {
  fig <- ggdraw()
  sub_fig_sz <- 1 / sub_plots
  combinations <- expand.grid(x = 0:(sub_plots[1]-1), y = 0:(sub_plots[2]-1))
  
  i <- 1
  
  for (j in seq_along(Coverage_thresholds)) {
    for (k in seq_along(Diff_thresholds)) {
      
      motif_path <- file.path(Path, paste0("cov_thresh_", Coverage_thresholds[j]),
                              paste0("win_", L), mut_type,
                              paste0("kmer_", kmer),
                              paste0("CNNC_", CNNC_motif_group))
      
      CNNC_dist_file <- paste0("CNNC_dist_thresh_", Diff_thresholds[k], ".tsv")
      CNNC_dist_randomers_file <- paste0("CNNC_dist_randomers_thresh_", Diff_thresholds[k], ".tsv")
      
      Title <- paste0("Cov_thresh = ", Coverage_thresholds[j], ", Diff_thresh = ", Diff_thresholds[k])
      
      sub_plot <- plot_CNNC_dist_vs_cum_frac(
        motif_path, CNNC_dist_file, CNNC_dist_randomers_file,
        do_abs_cum_freq, do_diff, Title,
        y_lim, y_tick_interval, x_tick_interval
      )
      
      fig <- fig + draw_plot(sub_plot,
                             x = (j - 1) * sub_fig_sz[1],
                             y = (k - 1) * sub_fig_sz[2],
                             width = sub_fig_sz[1],
                             height = sub_fig_sz[2])
      
      i <- i + 1
    }
  }
  
  ### Save composite figure
  ggsave(fig_name,
         plot = fig,
         width = fig_sz[1] * sub_plots[1],
         height = fig_sz[2] * sub_plots[2],
         limitsize = FALSE,
         device = cairo_pdf)
}

######################## Create folder ######################################################
create_folder <- function(Path) {
  ifelse(!dir.exists(Path), dir.create(Path, recursive = TRUE), FALSE)
}


##################################################################################################################################
### set the project path as the work directory ###
Project_path = "/data/groups/anko/mohamed/oxford_nanopore/"
Tools_path = "Tools/DIRCLIP_analysis/5_Motif_analysis/4_Find_motifs_from_mutation_sites/"
setwd(Project_path)

################################ Paths ##########################################################################################################
# Set root project directory
Project_path <- "/data/groups/anko/mohamed/oxford_nanopore/"
setwd(Project_path)

# Experimental variables
paths <- c("Results/Nanopore/3_SRSF3_DIRCLIP/Motif_analysis/minimap/mm9/Mauro/CLIP_cDNA_merge/")
ref_subfolder <- c("ref_TotalRNA_cDNA_UV", "ref_TotalRNA_cDNA_no_UV")
CNNC_motif_groups <- c("group14", "group15")
subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

Coverage_thresholds <- c(50)
Diff_thresholds <- c(0.03)
sub_plots <- c(1, 1)

kmer <- 5
Lengths <- c(100)
mut_types <- c("regions_centre", "sites")

# Iterate through groups and parameters
for (CNNC_motif_group in CNNC_motif_groups) {
  for (subfolder1 in subfolders1) {
    for (subfolder2 in subfolders2) {
      Motifs_path <- file.path(paths, subfolder1, subfolder2, ref_subfolder)
      
      for (L in Lengths) {
        for (mut_type in mut_types) {
          for (i in seq_along(Motifs_path)) {
            
            motif_path <- Motifs_path[i]
            print(motif_path)
            
            Res_path <- file.path(motif_path, "plot_CNNC_dist", mut_type, paste0("CNNC_", CNNC_motif_group))
            create_folder(Res_path)
            
            ### 1. Cumulative frequency plot
            plot_CNNC_dist_vs_cum_frac_joined_plots(
              motif_path, CNNC_motif_group, Coverage_thresholds, Diff_thresholds,
              L, mut_type, kmer, TRUE, FALSE, c(), 0.15, 10,
              fig_name = file.path(Res_path, paste0("CNNC_dist_cumfreq_L_", L, "_sites.pdf")),
              fig_sz = c(7, 7), sub_plots
            )
            
            ### 2. Cumulative difference plot
            plot_CNNC_dist_vs_cum_frac_joined_plots(
              motif_path, CNNC_motif_group, Coverage_thresholds, Diff_thresholds,
              L, mut_type, kmer, TRUE, TRUE, c(), 0.15, 10,
              fig_name = file.path(Res_path, paste0("CNNC_dist_cumfreq_diff_L_", L, "_sites.pdf")),
              fig_sz = c(7, 7), sub_plots
            )
            
            ### 3. Raw frequency plot
            plot_CNNC_dist_vs_cum_frac_joined_plots(
              motif_path, CNNC_motif_group, Coverage_thresholds, Diff_thresholds,
              L, mut_type, kmer, FALSE, FALSE, c(), 0.005, 25,
              fig_name = file.path(Res_path, paste0("CNNC_dist_freq_L_", L, "_sites.pdf")),
              fig_sz = c(7, 7), sub_plots
            )
          }
        }
      }
    }
  }
}


#####################################################################################################################################


