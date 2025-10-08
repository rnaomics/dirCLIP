################################################################################
# Mutation Site Distribution Analysis
# Author: Mohamed B
# Purpose: Load, normalize, and visualize mutation sites across gene bodies.
################################################################################

### 📁 Set Working Directory and Library Paths
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### 📚 Load Required Libraries
library(Rsubread)
library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggcorrplot)
library(scales)
library(GenomicRanges)
library(S4Vectors)
library(ComplexUpset)
library(patchwork)
library(rtracklayer)
library(tidyverse, warn.conflicts = FALSE)
library(ggupset)

### ⚙ Load Helper Functions
# source("Tools/DIRCLIP_analysis/8_Mutation_profiling/1_mm9_genome/count_genewise_sites_func.R")

### 🧬 Function to Load Mutation Sites as GRanges
load_sites <- function(path, left_offset, right_offset) {
  sites <- read.delim(path, header = FALSE)
  sites$strand <- "*"
  sites <- sites[!is.na(sites$V2), ]
  sites <- sites[!grepl("random", sites$V1), ]
  
  return(GRanges(seqnames = sites$V1,
                 ranges = IRanges(sites$V2 - left_offset, sites$V2 + right_offset),
                 strand = sites$strand))
}

### 🧬 Function to Load Mutation Regions
load_regions <- function(path, left_offset, right_offset) {
  regions <- read.delim(path, header = TRUE)
  regions$strand <- "*"
  regions <- regions[!is.na(regions$start), ]
  regions <- regions[!grepl("random", regions$chr), ]
  
  return(GRanges(seqnames = regions$chr,
                 ranges = IRanges(regions$start - left_offset, regions$end + right_offset),
                 strand = regions$strand))
}

### 📍 Function to Center Genomic Ranges
get_regions_centre <- function(gr) {
  centered <- gr
  is_range <- width(centered) > 1
  midpoints <- round((start(centered)[is_range] + end(centered)[is_range]) / 2)
  start(centered)[is_range] <- midpoints
  end(centered)[is_range] <- midpoints
  return(centered)
}

### 📐 Normalize Mutation Sites by Gene Length
normalise_by_gene_length <- function(mut_gr, gene_gr) {
  hits <- findOverlaps(mut_gr, gene_gr)
  
  mut_df <- data.frame(mut_gr[queryHits(hits)])[, c("seqnames", "start")]
  colnames(mut_df) <- c("Mut_chr", "Mut_loc")
  
  gene_df <- data.frame(gene_gr[subjectHits(hits)])[, c("seqnames", "start", "end", "width", "strand", "gene_id", "gene_type", "gene_name")]
  combined <- cbind(mut_df, gene_df)
  
  combined <- combined %>%
    mutate(rel_mut_loc = ifelse(strand == "+",
                                Mut_loc - start,
                                end - Mut_loc),
           norm_loc = rel_mut_loc / width)
  
  return(combined)
}

### 📊 Main Function: Generate Patched Density Plots
patched_mutsites_distribution_plot <- function(Mut_path, filt_by, type, gene_GR,
                                               Coverage_thresholds, Diff_thresholds,
                                               fig_name, fig_sz, sub_plots) {
  
  figure_list <- list()
  text_size <- 16
  Diff_thresholds <- sort(Diff_thresholds, decreasing = TRUE)
  
  for (diff in Diff_thresholds) {
    for (cov in Coverage_thresholds) {
      
      mut_file <- if (type == "sites") {
        paste0("cov_thresh_", cov, "/", filt_by, "/Pos_mut_thresh_", diff, ".tsv")
      } else {
        paste0("cov_thresh_", cov, "/", filt_by, "/Pos_mut_thresh_", diff, "_region.tsv")
      }
      
      mut_GR <- if (type == "sites") {
        load_sites(file.path(Mut_path, mut_file), 0, 0)
      } else {
        get_regions_centre(load_regions(file.path(Mut_path, mut_file), 0, 0))
      }
      
      mut_df <- normalise_by_gene_length(mut_GR, gene_GR)
      title <- paste0("Cov_thresh = ", cov, ", Diff_thresh = ", diff)
      plot <- make_density_plot(mut_df, title, text_size)
      
      figure_list[[length(figure_list) + 1]] <- plot
    }
  }
  
  final_fig <- wrap_plots(figure_list, ncol = sub_plots[1], nrow = sub_plots[2], byrow = TRUE)
  
  ggsave(fig_name,
         plot = final_fig,
         width = fig_sz[1] * sub_plots[1],
         height = fig_sz[2] * sub_plots[2],
         limitsize = FALSE,
         device = cairo_pdf)
}


################################################################################

### 📉 Function to Draw Normalized Density Plot
make_density_plot <- function(df, Title, text_size) {
  
  fig <- ggplot(df, aes(x = norm_loc)) +
    geom_density(alpha = 0.5, size = 1, color = "black", fill = "skyblue") +
    labs(title = Title,
         x = "Mut sites locs normalised by gene length",
         y = "Density") +
    theme_minimal(base_size = text_size) +
    theme(
      text = element_text(size = text_size),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.text.x = element_text(size = text_size, color = "black"),
      axis.text.y = element_text(size = text_size, color = "black"),
      axis.title.x = element_text(size = text_size),
      axis.title.y = element_text(size = text_size),
      plot.title  = element_text(size = text_size, hjust = 0.5),
      axis.ticks  = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    )
  
  ### Calculate Tick Breaks ###
  y_max <- max(density(df$norm_loc)$y)
  y_ticks <- pretty(c(0, y_max), n = 5)
  x_ticks <- seq(0, 1, by = 0.2)
  
  fig <- fig + scale_y_continuous(limits = c(min(y_ticks), max(y_ticks)),
                                  breaks = y_ticks,
                                  expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(limits = c(0, 1),
                       breaks = x_ticks,
                       expand = expansion(mult = c(0, 0.02)))
  
  return(fig)
}

  
###############################################################################################################################

### 🧬 Paths and Parameters
Mut_paths <- c(
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/"
)

ref_subfolder <- c("ref_TotalRNA_cDNA_no_UV", "ref_TotalRNA_cDNA_UV",
                   "ref_P19_WT_TotalRNA_cDNA_no_UV", "ref_TotalRNA_cDNA_UV",
                   "ref_P19_WT_TotalRNA_cDNA_no_UV", "ref_TotalRNA_cDNA_UV")

subfolders1 <- c("q_8_Q_8")
subfolders2 <- c("2_Inc_del", "1_No_del")
filt_by <- "FDR"

gtf_file <- "genome/mm9/Annotation/gencode.vM1.annotation.gtf"
gtf_GR <- import(gtf_file)
gene_GR <- gtf_GR[gtf_GR$type == "gene"]

Coverage_thresholds <- c(50)
Diff_thresholds <- c(0.03)
sub_plots <- c(1, 1)
fig_sz <- c(6, 7)

################################################################################################################################
for (i in seq_along(Mut_paths)) {
  for (subfolder1 in subfolders1) {
    for (subfolder2 in subfolders2) {
      
      Mut_path <- file.path(Mut_paths[i], subfolder1, subfolder2, ref_subfolder[i])
      Res_path <- file.path(Mut_path, filt_by)
      if (!dir.exists(Res_path)) dir.create(Res_path, recursive = TRUE)
      
      # 📊 Mutation Site Distribution
      patched_mutsites_distribution_plot(Mut_path = Mut_path, filt_by = filt_by, type = "sites", gene_GR,
                                         Coverage_thresholds, Diff_thresholds,
                                         fig_name = file.path(Res_path, "Mutsites_distribution_across_gene_length.pdf"),
                                         fig_sz = fig_sz, sub_plots = sub_plots)
      
      # 📊 Mutation Region Distribution (Centered)
      patched_mutsites_distribution_plot(Mut_path = Mut_path, filt_by = filt_by, type = "regions", gene_GR,
                                         Coverage_thresholds, Diff_thresholds,
                                         fig_name = file.path(Res_path, "Mutregions_distribution_across_gene_length.pdf"),
                                         fig_sz = fig_sz, sub_plots = sub_plots)
    }
  }
}

#############################################################################################################################
