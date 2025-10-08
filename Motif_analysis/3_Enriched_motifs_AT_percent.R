# 🧭 Set R package library path
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

# 📚 Load required libraries
library(ShortRead)
library(stringr)
library(stringi)
library(ggplot2)
library(Biostrings)
library(tidyverse)
library(writexl)
library(dplyr)
library(cowplot)
require(ggseqlogo)  # install via devtools::install_github("omarwagih/ggseqlogo")

# 🔗 Reference resources
# Motif plotting: motifStack vignette
# Motif discovery tutorials: 
#   - http://cbdm-01.zdv.uni-mainz.de/~jibnsale/teaching/motif_analysis_exercises.html
#   - http://rsat.france-bioinformatique.fr/teaching/peak-motifs.cgi


#### 🔤 Function: Calculate AT Percentage of Motif
calculate_AT_percentage <- function(motif) {
  at_count <- str_count(motif, "[AT]")
  total_length <- nchar(motif)
  return((at_count / total_length) * 100)
}

###🌟 Function: Get Top Enriched Motifs and AT Composition
get_enriched_motifs_AT_percentage <- function(motif_path_tgt, Sample_tgt, L, mut_type, Top_N,
                                              cov_thresh, diff_thesh, LFC_thresh, Count_thresh, kmers) {
  
  motif_path1 <- paste0(motif_path_tgt, 
                        "cov_thresh_", cov_thresh, 
                        "/win_", L, "/", 
                        mut_type, 
                        "/kmer_", kmers, "/")
  
  print(motif_path1)
  print(diff_thesh)
  
  # 📂 Read motif file
  motif_file1 <- paste0(motif_path1, "Pos_zscore_thresh_", diff_thesh, ".tsv")
  Motifs1 <- read.delim(motif_file1)
  
  # 🏅 Select top motifs
  Top_motifs <- Motifs1[1:Top_N, ]
  
  # 🔬 Annotate with AT%
  Top_motifs$AT_percentage <- sapply(Top_motifs$motifs, calculate_AT_percentage)
  
  # 📊 Construct output data frame
  df1 <- data.frame(AT_precent = Top_motifs$AT_percentage,
                    Sample = Sample_tgt)
  
  return(df1)
}


#### 📈 Function: Generate Box and Violin Plot of AT Content

make_box_plot <- function(df, Colors, Samples, Top_N, Coverage_thresholds, Diff_thresholds, fig_name, fig_sz) {
  
  df$Sample <- factor(df$Sample, levels = Samples)
  
  # ➕ Add sample size annotations
  n_labels <- df %>%
    group_by(Sample) %>%
    summarise(n = n(), max_value = max(AT_precent))
  
  # 🖋️ Styling
  Title <- paste0("cov thresh = ", Coverage_thresholds, 
                  ", diff thresh = ", Diff_thresholds, 
                  ", n = ", Top_N)
  text_size <- 16
  
  set.seed(123)
  
  fig <- ggplot(df, aes(x = Sample, y = AT_precent, fill = Sample)) +
    geom_violin(width = 0.8, linewidth = 0.6, color = "black", alpha = 0.85) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = Colors) +
    labs(title = Title, x = "Sample", y = "AT % in enriched motifs") +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.text.x = element_text(size = text_size, color = "black"),
      axis.title.x = element_text(size = text_size),
      axis.text.y = element_text(size = text_size, color = "black"),
      axis.title.y = element_text(size = text_size),
      legend.position = "none",
      plot.title = element_text(size = text_size, hjust = 0),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    )
  
  # 📏 Set Y-axis range
  y_tick_intervals <- 20
  y_max <- ceiling(max(df$AT_precent) / y_tick_intervals) * y_tick_intervals
  y_ticks <- seq(0, y_max, by = y_tick_intervals)
  
  fig <- fig + scale_y_continuous(limits = c(0, y_max), breaks = y_ticks,
                                  expand = expansion(mult = c(0, 0.02)))
  
  # 🖨️ Save to file
  ggsave(fig_name, plot = fig,
         width = fig_sz[1], height = fig_sz[2],
         limitsize = FALSE, device = cairo_pdf)
}

################################################################################################

# Set the project root directory and analysis tool path
Project_path <- "/data/groups/anko/mohamed/oxford_nanopore/"
Tools_path   <- "Tools/DIRCLIP_analysis/5_Motif_analysis/4_Find_motifs_from_mutation_sites/"
setwd(Project_path)

# Paths to CLIP-cDNA motif folders for each condition
Paths <- list(
  c("Results/Nanopore/3_SRSF3_DIRCLIP/Motif_analysis/minimap/mm9/CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Hmga1B_CLIP_cDNA_merge/"),
  
  c("Results/Nanopore/3_SRSF3_DIRCLIP/Motif_analysis/minimap/mm9/CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/CLIP_cDNA_merge/",
    "Results/Nanopore/2_Hmga1_DIRCLIP/Motif_analysis/minimap/mm9/Hmga1B_CLIP_cDNA_merge/")
)

# Matching reference subfolders (±UV)
ref_subfolders <- list(
  c("ref_TotalRNA_cDNA_no_UV", "ref_P19_WT_TotalRNA_cDNA_no_UV", "ref_P19_WT_TotalRNA_cDNA_no_UV"),
  c("ref_TotalRNA_cDNA_UV",    "ref_TotalRNA_cDNA_UV",           "ref_TotalRNA_cDNA_UV")
)

# Sample labels per condition
Samples <- list(
  c("SRSF3", "Hmga1a", "Hmga1b"),
  c("SRSF3", "Hmga1a", "Hmga1b")
)

# Output result paths
Result_paths <- c(
  "Results/Nanopore/3_SRSF3_DIRCLIP/Motif_analysis/minimap/mm9/CLIP_cDNA_merge/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Motif_analysis/minimap/mm9/CLIP_cDNA_merge/"
)

# Corresponding result subfolders
Result_subfolders <- c("ref_TotalRNA_cDNA_no_UV", "ref_TotalRNA_cDNA_UV")

# Sample color palette
Colors <- c("SRSF3" = "#66c2a5", "Hmga1a" = "#E26A6D", "Hmga1b" = "#FC9272")

##### ⚙️ Analysis Parameters
# Folder filters
subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")

# Analysis thresholds and parameters
Coverage_thresholds <- c(50)
Diff_thresholds     <- c(0.03)
kmers               <- c(5)
Lengths             <- c(100)
mut_types           <- c("regions")  # Mutation annotation types

##### 🔁 Main Loop: Extract Motif AT Content and Visualize

for (i in seq_along(Paths)) {
  for (subfolder1 in subfolders1) {
    for (subfolder2 in subfolders2) {
      
      # 🧭 Construct full motif paths for each sample
      motif_path1 <- paste0(Paths[[i]][1], subfolder1, "/", subfolder2, "/", ref_subfolders[[i]][1], "/")
      motif_path2 <- paste0(Paths[[i]][2], subfolder1, "/", subfolder2, "/", ref_subfolders[[i]][2], "/")
      motif_path3 <- paste0(Paths[[i]][3], subfolder1, "/", subfolder2, "/", ref_subfolders[[i]][3], "/")
      
      Sample1 <- Samples[[i]][1]
      Sample2 <- Samples[[i]][2]
      Sample3 <- Samples[[i]][3]
      
      Res_path <- paste0(Result_paths[i], subfolder1, "/", subfolder2, "/", Result_subfolders[i], "/")
      
      # 🚀 Iterate over sequence lengths and mutation annotation types
      for (L in Lengths) {
        for (mut_type in mut_types) {
          
          # STEP 1: Compute AT% for enriched motifs
          Enriched_AT1 <- get_enriched_motifs_AT_percentage(motif_path1, Sample1,
                                                            L, mut_type, Top_N = 50,
                                                            Coverage_thresholds, Diff_thresholds,
                                                            LFC_thresh = 1, Count_thresh = 20, kmers)
          
          Enriched_AT2 <- get_enriched_motifs_AT_percentage(motif_path2, Sample2,
                                                            L, mut_type, Top_N = 50,
                                                            Coverage_thresholds, Diff_thresholds,
                                                            LFC_thresh = 1, Count_thresh = 20, kmers)
          
          Enriched_AT3 <- get_enriched_motifs_AT_percentage(motif_path3, Sample3,
                                                            L, mut_type, Top_N = 50,
                                                            Coverage_thresholds, Diff_thresholds,
                                                            LFC_thresh = 1, Count_thresh = 20, kmers)
          
          # 🧬 Combine across all samples
          Enriched_AT_all <- rbind(Enriched_AT1, Enriched_AT2, Enriched_AT3)
          
          # STEP 2: Plot AT% distributions
          make_box_plot(
            df        = Enriched_AT_all,
            Colors    = Colors,
            Samples   = Samples[[i]],
            Top_N     = 50,
            Coverage_thresholds = Coverage_thresholds,
            Diff_thresholds     = Diff_thresholds,
            fig_name  = paste0(Res_path, "Enriched_motifs_AT_percent_L_", 
                               L, "_", mut_type, "_", Sample1, "_", Sample2, ".pdf"),
            fig_sz    = c(5, 7)
          )
        }
      }
    }
  }
}

################################################################################################