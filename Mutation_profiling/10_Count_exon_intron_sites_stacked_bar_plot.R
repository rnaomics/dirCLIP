################################################################################
# annotation_utilities.R
# Author: Mohamed B
# Description:
# Utility functions for processing GTF annotations, extracting genomic features
# (exons, introns, UTRs), loading CLIP sites/regions, and annotating them.
# Includes bar plot helper for visualizing annotations.
################################################################################

### 📂 Set working directory and library paths ###
setwd("/data/groups/anko/mohamed/oxford_nanopore/")
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### 📦 Load required libraries ###
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(GGally)
library(ggcorrplot)
library(scales)
library(S4Vectors)
library(matrixStats)
library(reshape2)
library(RColorBrewer)

################################################################################
# 🧬 1. Extract exons and introns from GTF annotation
################################################################################
get_exons_introns_from_gtf <- function(gtf_file, out_file) {
  gtf <- import(gtf_file)
  gtf <- gtf[!grepl("_random", seqnames(gtf))]  # Remove random chromosomes
  
  transcripts <- gtf[gtf$type == "transcript"]
  transcripts_reduced <- reduce(split(transcripts, transcripts$transcript_id))
  
  exons <- gtf[gtf$type == "exon"]
  exons_by_transcript <- split(exons, exons$transcript_id)
  exons_reduced <- reduce(exons)
  
  # Derive introns by subtracting exons from transcripts
  introns_list <- lapply(seq_along(transcripts_reduced), function(i) {
    GenomicRanges::setdiff(transcripts_reduced[[i]], reduce(exons_by_transcript[[i]]))
  })
  
  introns <- do.call(c, introns_list)
  introns <- introns[!overlapsAny(introns, exons_reduced)]
  introns_reduced <- reduce(introns)
  
  # Convert to data frames and save
  save(data.frame(exons_reduced), data.frame(introns_reduced), file = out_file)
}

################################################################################
# 🧬 2. Extract and classify UTRs (5' and 3') from GTF
################################################################################
get_UTR_from_gtf <- function(gtf_file, out_file) {
  gtf <- import(gtf_file)
  gtf <- gtf[!grepl("_random", seqnames(gtf))]
  utrs <- gtf[gtf$type == "UTR"]
  
  start_codons <- gtf[gtf$type == "start_codon"]
  stop_codons  <- gtf[gtf$type == "stop_codon"]
  utrs$UTR_type <- NA
  
  for (tx in unique(utrs$transcript_id)) {
    tx_utrs <- utrs[utrs$transcript_id == tx]
    tx_strand <- unique(as.character(strand(tx_utrs)))
    start_codon <- start_codons[start_codons$transcript_id == tx]
    stop_codon  <- stop_codons[stop_codons$transcript_id == tx]
    if (length(start_codon) == 0 || length(stop_codon) == 0) next
    
    if (tx_strand == "+") {
      start_pos <- max(start(start_codon))
      stop_pos  <- min(end(stop_codon))
      tx_utrs$UTR_type[start(tx_utrs) < start_pos] <- "5UTR"
      tx_utrs$UTR_type[end(tx_utrs) > stop_pos]   <- "3UTR"
    } else {
      start_pos <- min(end(start_codon))
      stop_pos  <- max(start(stop_codon))
      tx_utrs$UTR_type[end(tx_utrs) > start_pos] <- "5UTR"
      tx_utrs$UTR_type[start(tx_utrs) < stop_pos] <- "3UTR"
    }
    utrs[utrs$transcript_id == tx]$UTR_type <- tx_utrs$UTR_type
  }
  
  save(utrs, utrs, file = out_file)
}

################################################################################
# 📥 3. Load CLIP sites as GenomicRanges
################################################################################
load_sites <- function(path) {
  sites <- read.delim(path, header = FALSE)
  gr <- GRanges(seqnames = sites$V1, ranges = IRanges(start = sites$V2, end = sites$V2))
  gr <- gr[!grepl("_random", seqnames(gr))]
  return(gr)
}

################################################################################
# 📥 4. Load CLIP regions as GenomicRanges
################################################################################
load_regions <- function(path) {
  regions <- read.delim(path, header = TRUE)
  gr <- GRanges(seqnames = regions$chr, ranges = IRanges(start = regions$start, end = regions$end))
  gr <- gr[!grepl("_random", seqnames(gr))]
  return(gr)
}

################################################################################
# 🧠 5. Annotate genomic positions (exon/intron/UTR/intergenic)
################################################################################
annotate_exons_introns_UTRs <- function(gr, exons, introns, utrs, gtf) {
  gr$annotation <- "Intergenic"
  
  gr$annotation[queryHits(findOverlaps(gr, gtf))]   <- "Genic"
  gr$annotation[queryHits(findOverlaps(gr, exons))] <- "Exon"
  gr[gr$annotation == "Genic"]$annotation <- "Intron"
  
  gr$annotation[queryHits(findOverlaps(gr, utrs[utrs$UTR_type == "3UTR"]))] <- "3UTR"
  gr$annotation[queryHits(findOverlaps(gr, utrs[utrs$UTR_type == "5UTR"]))] <- "5UTR"
  
  return(gr)
}

################################################################################
# 🔬 6. Annotate gene types from GTF
################################################################################
annotate_gene_type <- function(gr, gtf) {
  hits <- findOverlaps(gr, gtf)
  gene_types <- mcols(gtf)$gene_type[subjectHits(hits)]
  mcols(gr)$gene_type <- "Not_annotated"
  mcols(gr)$gene_type[queryHits(hits)] <- gene_types
  return(gr)
}

################################################################################
# 📊 7. Customizable grouped bar plot for annotation results
################################################################################
make_bar_plot <- function(df, colors, title, y_lim, y_interval, x_lab, y_lab, legend_title) {
  text_size <- 16
  fig <- ggplot(df, aes(x = Group, y = y, fill = sub_group)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.6) +
    scale_fill_manual(values = colors) +
    labs(x = x_lab, y = y_lab, fill = legend_title, title = title) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.title = element_text(size = text_size),
      axis.text  = element_text(size = text_size, color = "black"),
      legend.title = element_text(size = text_size),
      legend.text  = element_text(size = text_size),
      axis.ticks = element_line(color = "black", size = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    )
  
  # Adjust Y-axis ticks
  if (length(y_lim) > 0) {
    y_min <- y_lim[1]
    y_max <- y_lim[2]
  } else {
    y_vals <- df %>% group_by(Group) %>% summarise(total_y = sum(y), .groups = "drop")
    y_max <- ceiling(max(y_vals$total_y) / y_interval) * y_interval
    y_min <- 0
  }
  
  fig <- fig + scale_y_continuous(
    limits = c(y_min, y_max),
    breaks = seq(y_min, y_max, by = y_interval),
    expand = expansion(mult = c(0, 0.02))
  )
  
  return(fig)
}

################################################################################

### 🧬 Annotate mutation sites/regions across coverage/difference thresholds ###
annotate_sites <- function(Mut_paths, filt_by, type, Coverage_thresholds, Diff_thresholds,
                           exons_GR, introns_GR, utrs_GR, gtf_GR) {
  ii <- 1
  annot_sites_list <- list()
  annot_genes_list <- list()
  Cov_thresh_list <- c()
  diff_thresh_list <- c()
  
  for (Cov_thresh in Coverage_thresholds) {
    for (diff_thresh in Diff_thresholds) {
      Cov_thresh_list[ii] <- Cov_thresh
      diff_thresh_list[ii] <- diff_thresh
      
      for (i in seq_along(Mut_paths)) {
        base_path <- paste0(Mut_paths[i], "cov_thresh_", Cov_thresh, "/", filt_by, "/")
        file_suffix <- if (type == "sites") paste0("Pos_mut_thresh_", diff_thresh, ".tsv") else
          paste0("Pos_mut_thresh_", diff_thresh, "_region.tsv")
        full_path <- paste0(base_path, file_suffix)
        
        sites_GR <- if (type == "sites") load_sites(full_path) else load_regions(full_path)
        sites_GR <- annotate_exons_introns_UTRs(sites_GR, exons_GR, introns_GR, utrs_GR, gtf_GR)
        sites_GR <- annotate_gene_type(sites_GR, gtf_GR)
        
        annot_counts <- table(sites_GR$annotation)
        gene_type_counts <- table(sites_GR$gene_type)
        
        # Fill missing categories with 0s
        annot_counts <- ifelse(annot_type %in% names(annot_counts), annot_counts[annot_type], 0)
        gene_type_counts <- ifelse(gene_type %in% names(gene_type_counts), gene_type_counts[gene_type], 0)
        
        # Accumulate across samples
        if (i == 1) {
          n_sites_annot_type <- annot_counts
          n_sites_gene_type  <- gene_type_counts
        } else {
          n_sites_annot_type <- rbind(n_sites_annot_type, annot_counts)
          n_sites_gene_type  <- rbind(n_sites_gene_type, gene_type_counts)
        }
      }
      
      rownames(n_sites_annot_type) <- Samples
      rownames(n_sites_gene_type)  <- Samples
      colnames(n_sites_annot_type) <- annot_type
      colnames(n_sites_gene_type)  <- gene_type
      
      annot_sites_list[[ii]] <- n_sites_annot_type
      annot_genes_list[[ii]] <- n_sites_gene_type
      ii <- ii + 1
    }
  }
  
  return(tibble(Cov_thresh_list, diff_thresh_list, annot_sites_list, annot_genes_list))
}

################################################################################
# 📊 Generate joined stacked bar plots per threshold combination
################################################################################
plot_genes_with_mutsites_overlap_with_iclip_joined_plots <- function(sites_annot, Coverage_thresholds, Diff_thresholds,
                                                                     Colors, y_lim, y_tick_interval, fig_name,
                                                                     fig_sz, sub_plots) {
  fig <- ggdraw()
  sub_fig_sz <- 1 / sub_plots
  
  for (k in seq_along(Coverage_thresholds)) {
    for (m in seq_along(Diff_thresholds)) {
      Cov_thresh <- Coverage_thresholds[k]
      Diff_thresh <- Diff_thresholds[m]
      
      id <- sites_annot$Cov_thresh_list == Cov_thresh & sites_annot$diff_thresh_list == Diff_thresh
      annot_sites <- sites_annot[id, ]$annot_sites_list[[1]]
      annot_sites <- annot_sites / rowSums(annot_sites)
      
      df <- melt(cbind(annot_sites, Group = rownames(annot_sites)),
                 id.vars = "Group", variable.name = "sub_group", value.name = "y")
      
      Title <- paste0("Cov_thresh = ", Cov_thresh, ", Diff_thresh = ", Diff_thresh)
      bar <- make_bar_plot(df, Colors, Title, y_lim, y_tick_interval,
                           x_lab = "Sample", y_lab = "Fraction of sites", legend_title = "Annotation")
      
      fig <- fig + draw_plot(bar,
                             x = (k - 1) * sub_fig_sz[1],
                             y = (m - 1) * sub_fig_sz[2],
                             width = sub_fig_sz[1],
                             height = sub_fig_sz[2])
    }
  }
  
  ggsave(fig_name, plot = fig,
         width  = fig_sz[1] * sub_plots[1],
         height = fig_sz[2] * sub_plots[2],
         limitsize = FALSE,
         device = cairo_pdf)
}


######### 🧬 Sample Paths and Setup #######################################################################
Samples <- c("Hmga1a", "Hmga1b")
Mutation_paths <- c("Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/",
                    "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/")

Res_path <- "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/q_8_Q_8/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/"
filt_by <- "FDR"
Coverage_thresholds <- c(50)
Diff_thresholds <- c(0.03)
sub_plots <- c(1, 1)


gtf_file  = "genome/mm9/Annotation/gencode.vM1.annotation.gtf" # Gencode annotation
exons_introns_RData = "genome/mm9/Annotation/gencode_exons_and_intron.RData"
UTR_RData = "genome/mm9/Annotation/gencode_UTR.RData"

##################### 1. Convert gtf_file to exons introns RData #################################################
get_exons_introns_flag = FALSE # TRUE
if (get_exons_introns_flag){
  get_exons_introns_from_gtf(gtf_file = gtf_file, out_file = exons_introns_RData)
}

#################### 2. Getting the 3' UTR and 5' UTR RData ######################################################
get_UTR_flag = FALSE # TRUE

if (get_UTR_flag){
  get_UTR_from_gtf(gtf_file = gtf_file, out_file = UTR_RData)
}


### 📁 Load Annotation Inputs ###
load(exons_introns_RData)
load(UTR_RData)
exons_GR   <- GRanges(exons_reduced)
introns_GR <- GRanges(introns_reduced)
utrs_GR    <- utrs[!is.na(utrs$UTR_type)]

gtf_GR <- import(gtf_file)
gtf_GR <- gtf_GR[gtf_GR$type == "gene" & !grepl("_random", seqnames(gtf_GR))]

gene_type  <- c(unique(gtf_GR$gene_type), "Not_annotated")
annot_type <- c("3UTR", "5UTR", "Exon", "Intergenic", "Intron")
Colors <- c("3UTR" = "#66c2a5", "5UTR" = "#fc8d62", "Exon" = "#8da0cb",
            "Intergenic" = "#e78ac3", "Intron" = "#a6d854")

### 🧠 Annotation Pipeline ###
Mut_sites_annot   <- annotate_sites(Mut_paths = Mutation_paths, filt_by, "sites", Coverage_thresholds, Diff_thresholds, exons_GR, introns_GR, utrs_GR, gtf_GR)
Mut_regions_annot <- annotate_sites(Mut_paths = Mutation_paths, filt_by, "regions", Coverage_thresholds, Diff_thresholds, exons_GR, introns_GR, utrs_GR, gtf_GR)

### 📊 Plot Site Annotations ###
plot_genes_with_mutsites_overlap_with_iclip_joined_plots(Mut_sites_annot,
                                                         Coverage_thresholds, Diff_thresholds, Colors, y_lim = c(), y_tick_interval = 0.2,
                                                         fig_name = paste0(Res_path, filt_by, "/Mut_sites_annotation_stacked_bar.pdf"),
                                                         fig_sz = c(6, 7), sub_plots = sub_plots)

### 📊 Plot Region Annotations ###
plot_genes_with_mutsites_overlap_with_iclip_joined_plots(Mut_regions_annot,
                                                         Coverage_thresholds, Diff_thresholds, Colors, y_lim = c(), y_tick_interval = 0.2,
                                                         fig_name = paste0(Res_path, filt_by, "/Mut_regions_annotation_stacked_bar.pdf"),
                                                         fig_sz = c(6, 7), sub_plots = sub_plots)

###########################################################################################################

#  Suggestions for Documentation or Supplementary Methods

# Genomic features were extracted using GTF annotations from mm9, where transcripts and exons were parsed to define introns. UTRs were classified based on codon positions and strand orientation. CLIP datasets were loaded as GenomicRanges and annotated using overlapping gene features. Site classification into exons, introns, UTRs, and intergenic regions was visualized using custom bar plots generated with ggplot2.


