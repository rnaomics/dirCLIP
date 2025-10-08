################################################################################
# mutation_overlap_analysis.R
# Author: Mohamed B
# Description:
# Loads site/region mutations and iCLIP tags, applies positional offsets,
# and computes overlap matrices across transcriptome samples.
################################################################################

### 📁 Set working environment ###
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### 📚 Load required libraries ###
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
# library(ComplexUpset)
# library(UpSetR)
library(patchwork)
library(ggbreak)
library(tidyverse, warn.conflicts = FALSE)
library(ggupset)  # Your custom package?

### 🔍 Load iCLIP tags with directional offsets ###
load_CLIP_tags <- function(iclip_path, Peak_caller, Left_offset, Right_offset) {
  if (Peak_caller == "Michael_unknown") {
    iclip_bed <- read.table(iclip_path, skip = 1)
    iclip_bed$Strand <- ifelse(iclip_bed$V4 < 0, "-", "+")
    
    # Apply strand-specific offsets
    iclip_bed$V2[iclip_bed$Strand == "+"] <- iclip_bed$V2[iclip_bed$Strand == "+"] - Left_offset
    iclip_bed$V3[iclip_bed$Strand == "+"] <- iclip_bed$V3[iclip_bed$Strand == "+"] + Right_offset
    iclip_bed$V2[iclip_bed$Strand == "-"] <- iclip_bed$V2[iclip_bed$Strand == "-"] - Right_offset
    iclip_bed$V3[iclip_bed$Strand == "-"] <- iclip_bed$V3[iclip_bed$Strand == "-"] + Left_offset
    
    return(GRanges(iclip_bed$V1, IRanges(iclip_bed$V2, iclip_bed$V3), strand = iclip_bed$Strand))
  } else {
    iclip_bed <- read.delim(iclip_path, header = FALSE)
    
    # Apply strand-specific offsets
    iclip_bed$V2[iclip_bed$V6 == "+"] <- iclip_bed$V2[iclip_bed$V6 == "+"] - Left_offset
    iclip_bed$V3[iclip_bed$V6 == "+"] <- iclip_bed$V3[iclip_bed$V6 == "+"] + Right_offset
    iclip_bed$V2[iclip_bed$V6 == "-"] <- iclip_bed$V2[iclip_bed$V6 == "-"] - Right_offset
    iclip_bed$V3[iclip_bed$V6 == "-"] <- iclip_bed$V3[iclip_bed$V6 == "-"] + Left_offset
    
    return(GRanges(iclip_bed$V1, IRanges(iclip_bed$V2, iclip_bed$V3), strand = iclip_bed$V6))
  }
}

### 🧬 Load mutation sites as GRanges ###
load_sites <- function(Path, Left_offset, Right_offset) {
  Sites <- read.delim(Path, header = FALSE)
  Sites$strand <- "*"
  Sites <- Sites[!is.na(Sites$V2),]
  
  # Exclude entries from random chromosomes
  id <- grep("random", Sites$V1)
  if (!isEmpty(id)) {
    Sites <- Sites[-id, ]
  }
  
  return(GRanges(
    seqnames = Sites[,"V1"],
    ranges = IRanges(Sites[,"V2"] - Left_offset, Sites[,"V2"] + Right_offset),
    strand = Sites$strand
  ))
}

### 🧬 Load mutation regions as GRanges ###
load_regions <- function(Path, Left_offset, Right_offset) {
  Sites <- read.delim(Path, header = TRUE)
  Sites$strand <- "*"
  Sites <- Sites[!is.na(Sites$start),]
  
  id <- grep("random", Sites$chr)
  if (!isEmpty(id)) {
    Sites <- Sites[-id, ]
  }
  
  return(GRanges(
    seqnames = Sites[,"chr"],
    ranges = IRanges(Sites[,"start"] - Left_offset, Sites[,"end"] + Right_offset),
    strand = Sites$strand
  ))
}

### 🧮 Compute overlap matrix across mutation samples & iCLIP ###
find_overlaps <- function(Mut_sites_paths, mut_file, type, Left_offset, Right_offset, Transcritome_refs, iCLIP_GR) {
  Mut_GR_list <- list()
  
  for (i in seq_along(Mut_sites_paths)) {
    path <- paste0(Mut_sites_paths[i], mut_file)
    
    Mut_GR <- if (type == "sites") {
      load_sites(path, Left_offset, Right_offset)
    } else {
      load_regions(path, Left_offset, Right_offset)
    }
    
    Mut_GR_list[[i]] <- Mut_GR
  }
  
  # Combine all GRanges and remove duplicates
  all_gr <- unique(do.call(c, Mut_GR_list))
  
  # Optionally include iCLIP tags
  if (length(iCLIP_GR) > 0) {
    Mut_GR_list[[length(Mut_GR_list) + 1]] <- iCLIP_GR
  }
  
  # Build overlap data frame
  df <- data.frame(ID = seq_along(all_gr))
  
  for (i in seq_along(Mut_GR_list)) {
    hits <- findOverlaps(all_gr, Mut_GR_list[[i]])
    overlap_logical <- seq_along(all_gr) %in% queryHits(hits)
    
    df[[Transcritome_refs[i]]] <- overlap_logical
  }
  
  return(df)
}

########################## Patching the upset plot #############################

patched_upset_plot <- function(Mut_sites_paths, filt_by, type, iCLIP_GR, 
                               Left_offset, Right_offset, Transcritome_refs, 
                               Sample_interest, Coverage_thresholds, Diff_thresholds, 
                               fig_name, fig_sz, sub_plots) {
  
  figure_list <- list()
  text_size <- 16
  text_size2 <- 6
  Diff_thresholds <- sort(Diff_thresholds, decreasing = TRUE)
  
  for (diff_thresh in Diff_thresholds) {
    for (cov_thresh in Coverage_thresholds) {
      
      # 🔗 Construct filename for mutation input
      mut_file <- if (type == "sites") {
        paste0("cov_thresh_", cov_thresh, "/", filt_by, "/Pos_mut_thresh_", diff_thresh, ".tsv")
      } else {
        paste0("cov_thresh_", cov_thresh, "/", filt_by, "/Pos_mut_thresh_", diff_thresh, "_region.tsv")
      }
      
      # 📎 Add iCLIP reference if provided
      refs <- Transcritome_refs
      if (length(iCLIP_GR) > 0) {
        refs <- c(refs, "iCLIP")
      }
      
      # 🧬 Get binary matrix of overlaps
      df <- find_overlaps(Mut_sites_paths, mut_file, type, Left_offset, Right_offset, refs, iCLIP_GR)
      df[, -1] <- lapply(df[, -1], function(x) ifelse(x, 1, 0))  # Convert to binary
      
      # 🔖 Optional query annotation setup
      annot_interaction_list <- list(refs)
      for (ref in refs) {
        annot_interaction_list[[length(annot_interaction_list) + 1]] <- upset_query(
          intersect = c(Sample_interest, ref),
          color = 'blue', fill = 'blue',
          only_components = c('intersections_matrix', 'Intersection size')
        )
      }
      
      # 🎨 Generate UpSet figure
      Title <- paste0("Cov_thresh = ", cov_thresh, ", Diff_thresh = ", diff_thresh)
      figure <- make_ggupset(df, Title, text_size, text_size2, y_lim = c(), y_tick_interval = 6000)
      
      figure_list[[length(figure_list) + 1]] <- figure
    }
  }
  
  # 📦 Combine all plots in patchwork grid
  fig <- wrap_plots(figure_list, ncol = sub_plots[1], nrow = sub_plots[2], byrow = TRUE)
  
  # 💾 Save figure
  ggsave(fig_name,
         plot = fig,
         width = fig_sz[1] * sub_plots[1],
         height = fig_sz[2] * sub_plots[2],
         limitsize = FALSE,
         device = cairo_pdf)
}


########## Function to make ggupset plot #####################################
make_ggupset <- function(df, Title, text_size, text_size2, y_lim, y_tick_interval) {
  
  condition_list <- apply(df[, -1], 1, function(x) names(x)[x == 1])
  df_plot <- tibble(ID = df$ID, group = condition_list)
  
  fig <- ggplot(df_plot, aes(x = group)) +
    geom_bar() +
    scale_x_upset() +
    geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = text_size2) +
    labs(x = "group", y = "Intersection size", title = Title) +
    theme_minimal(base_size = text_size) +
    theme(
      text = element_text(size = text_size),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.text = element_text(size = text_size, color = "black"),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      plot.title = element_text(hjust = 0.5)
    )
  
  # 📏 Y-axis scaling
  if (length(y_lim) == 0) {
    y_ticks <- pretty(c(0, max(table(df[, -1])) * 1.1), n = 5)
  } else {
    y_min <- floor(y_lim[1] / y_tick_interval) * y_tick_interval
    y_max <- ceiling(y_lim[2] / y_tick_interval) * y_tick_interval
    y_ticks <- seq(y_min, y_max, by = y_tick_interval)
  }
  
  fig <- fig + scale_y_continuous(
    limits = range(y_ticks),
    breaks = y_ticks,
    expand = expansion(mult = c(0, 0.02))
  )
  
  return(fig)
}

################################################################################

# 🗂 Mutation paths grouped by condition/control
Mut_paths <- list(
  c("Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/1_No_del/ref_TotalRNA_cDNA_no_UV/",
    "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/1_No_del/ref_TotalRNA_cDNA_UV/"),
  
  c("Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_TotalRNA_cDNA_no_UV/",
    "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_TotalRNA_cDNA_UV/")
)

# 📌 Reference condition labels and sample of interest
Trans_refs <- list(c("-UV", "+UV"), c("-UV", "+UV"))
Sample_interest <- "+UV"

# 📁 Output base paths for figures
Res_paths <- c(
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/1_No_del/ref_TotalRNA_cDNA_no_UV/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_TotalRNA_cDNA_no_UV/"
)

# ⚙️ Filtering criteria
filt_by <- "FDR"  # or "pvalue"

# 🔬 iCLIP peak annotation setup
Peak_caller <- "iCount"
iclip_path <- c("Results/Illumina/iCLIP_SRSF3/Nextflow/SRSF3_ctrl_iCLIP/results/icount/SRSF3_merged.3nt_3nt.peaks.bed.gz")
Left_offset <- 0
Right_offset <- 0
iCLIP_GR <- load_CLIP_tags(iclip_path, Peak_caller, Left_offset, Right_offset)

# 📊 Plot parameters
Coverage_thresholds <- c(50)
Diff_thresholds <- c(0.03)
sub_plots <- c(1, 1)
fig_sz <- c(5, 7)

# 🔁 Iterate over mutation groups
for (i in seq_along(Mut_paths)) {
  
  # 📁 Construct output path
  Res_path <- file.path(Res_paths[i], filt_by)
  Mutation_paths <- Mut_paths[[i]]
  Transcritome_refs <- Trans_refs[[i]]
  
  Left_offset <- 5
  Right_offset <- 5
  
  # 🔹 Generate site-based UpSet plot
  patched_upset_plot(
    Mut_sites_paths = Mutation_paths,
    filt_by = filt_by,
    type = "sites",
    iCLIP_GR = c(),  # iCLIP_GR optionally included
    Left_offset = Left_offset,
    Right_offset = Right_offset,
    Transcritome_refs = Transcritome_refs,
    Sample_interest = Sample_interest,
    Coverage_thresholds = Coverage_thresholds,
    Diff_thresholds = Diff_thresholds,
    fig_name = file.path(Res_path, "Mutsites_olap_between_controls.pdf"),
    fig_sz = fig_sz,
    sub_plots = sub_plots
  )
  
  # 🔹 Generate region-based UpSet plot
  patched_upset_plot(
    Mut_sites_paths = Mutation_paths,
    filt_by = filt_by,
    type = "regions",
    iCLIP_GR = c(),
    Left_offset = Left_offset,
    Right_offset = Right_offset,
    Transcritome_refs = Transcritome_refs,
    Sample_interest = Sample_interest,
    Coverage_thresholds = Coverage_thresholds,
    Diff_thresholds = Diff_thresholds,
    fig_name = file.path(Res_path, "Mutregions_olap_between_controls.pdf"),
    fig_sz = fig_sz,
    sub_plots = sub_plots
  )
}

######################################################################################################

#  Suggested Notes for a Supplementary Methods Section

# Mutation profiles for SRSF3 DIRCLIP experiments were processed across two conditions: UV-treated and untreated transcriptomes, with and without deletion-based filtering. UpSet plots were generated using patched_upset_plot() to visualize intersecting mutation sites and regions. iCLIP peak data (iCount method) was optionally included but excluded from comparative plots shown. Coverage and mutation difference thresholds were set at 50 and 0.03 respectively.
