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
# Generates horizontal bar plots for annotation categories from mutation data,
# visualizing either absolute counts or fractions across annotation types.
################################################################################

make_bar_plot <- function(df, type, Color, Cov_thresh, Diff_thresh,
                          x_lim, x_tick_interval, x_lab, y_lab,
                          fig_name, fig_sz) {
  
  ### 🔧 Format data column based on plot type ###
  if (type == "fraction") {
    colnames(df)[colnames(df) == "Frac_sites"] <- "Value"
  } else if (type == "abs") {
    colnames(df)[colnames(df) == "N"] <- "Value"
  }
  
  ### 🎯 Filter dataframe for current thresholds ###
  df <- df[df$Cov_thresh == Cov_thresh & df$Diff_thresh == Diff_thresh, ]
  
  ### 🖼️ Plot title ###
  Title <- paste0("Cov thresh = ", Cov_thresh, ", Diff thresh = ", Diff_thresh)
  text_size <- 16
  
  ### 🧱 Construct horizontal bar plot ###
  fig <- ggplot(df, aes(y = Annot_type, x = Value, fill = Annot_type)) +
    geom_bar(stat = "identity", width = 0.7, color = "black", fill = Color, linewidth = 0.6) +
    labs(title = Title, x = x_lab, y = y_lab) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid         = element_blank(),
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "black", linewidth = 0.8),
      axis.line.y        = element_line(color = "black", linewidth = 0.8),
      axis.text          = element_text(size = text_size),
      axis.title         = element_text(size = text_size),
      plot.title         = element_text(hjust = 0.7),
      axis.ticks         = element_line(color = "black", size = 0.8),
      axis.ticks.length  = unit(0.25, "cm"),
      axis.ticks.x       = element_line(color = "black"),
      axis.ticks.y       = element_line(color = "black")
    )
  
  ### 📏 Set X-axis ticks and limits ###
  if (length(x_lim) > 0) {
    x_min <- x_lim[1]
    x_max <- x_lim[2]
  } else {
    x_min <- floor(min(df$Value) / x_tick_interval) * x_tick_interval
    x_max <- ceiling(max(df$Value) / x_tick_interval) * x_tick_interval
  }
  
  x_ticks <- seq(x_min, x_max, by = x_tick_interval)
  
  fig <- fig + scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = x_ticks,
    expand = expansion(mult = c(0, 0.02))  # Add 2% padding
  )
  
  ### 💾 Save to file ###
  ggsave(fig_name,
         plot = fig,
         width = fig_sz[1],
         height = fig_sz[2],
         limitsize = FALSE)
}


################################################################################
# Processes mutation data from Nanopore DIRCLIP experiments, annotates with
# genomic features (exons, introns, UTRs, gene types), and generates visual
# summaries (bar plots) of mutation site distribution across annotations.
################################################################################

# ⛳ Directory Setup & Metadata
Mutation_paths <- c(
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/CLIP_cDNA_merge/",
  "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/minimap/mm9/Mauro/Hmga1B_CLIP_cDNA_merge/"
)
Samples <- c("SRSF3", "SRSF3", "Hmga1a", "Hmga1b")
Colors <- c("#66c2a5", "#66c2a5", "#E26A6D", "#FC9272")

ref_subfolder <- c("ref_TotalRNA_cDNA_no_UV", "ref_TotalRNA_cDNA_UV", 
                   "ref_P19_WT_TotalRNA_cDNA_no_UV", "ref_P19_WT_TotalRNA_cDNA_no_UV")
subfolder1 <- "q_8_Q_8"
subfolders2 <- c("2_Inc_del", "1_No_del")
filt_by <- "FDR"

Coverage_thresholds <- c(50, 100)
Diff_thresholds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)

# 🧬 Load and Prepare Annotations
gtf_file <- "genome/mm9/Annotation/gencode.vM1.annotation.gtf"
exons_introns_RData <- "genome/mm9/Annotation/gencode_exons_and_intron.RData"
UTR_RData <- "genome/mm9/Annotation/gencode_UTR.RData"

get_exons_introns_from_gtf(gtf_file, out_file = exons_introns_RData)
get_UTR_from_gtf(gtf_file, out_file = UTR_RData)

load(exons_introns_RData)
exons_GR <- GRanges(exons_reduced)
introns_GR <- GRanges(introns_reduced)

load(UTR_RData)
utrs_GR <- utrs[!is.na(utrs$UTR_type)]

gtf_GR <- import(gtf_file)
gtf_GR <- gtf_GR[!grepl("_random", seqnames(gtf_GR))]
gtf_GR <- gtf_GR[gtf_GR$type == "gene"]

# 🔁 Loop Over Conditions and Annotate Mutations
for (subfolder2 in subfolders2) {
  Mut_paths <- paste0(Mutation_paths, subfolder1, "/", subfolder2, "/", ref_subfolder, "/")
  
  for (i in seq_along(Mut_paths)) {
    Color <- Colors[i]
    print(Mut_paths[i])
    
    # 🧠 Initialize Matrices for Counts
    init_array <- function(rows, cols) {
      array(dim = c(length(rows), length(cols)), 
            dimnames = list(rows, paste0("T_", cols)))
    }
    
    gene_types <- c(unique(gtf_GR$gene_type), "Not_annotated")
    Conditions <- c()
    ii <- 1
    
    for (k in seq_along(Coverage_thresholds)) {
      for (m in seq_along(Diff_thresholds)) {
        Cov_thresh <- Coverage_thresholds[k]
        Diff_thresh <- Diff_thresholds[m]
        Mut_path <- paste0(Mut_paths[i], "cov_thresh_", Cov_thresh, "/", filt_by, "/")
        
        pos_mut_file <- paste0("Pos_mut_thresh_", Diff_thresh, ".tsv")
        pos_mut_region <- paste0("Pos_mut_thresh_", Diff_thresh, "_region.tsv")
        
        # 📥 Load Mutation Data
        Pos_Mut_sites_GR <- load_sites(paste0(Mut_path, pos_mut_file))
        Pos_Mut_regions_GR <- load_regions(paste0(Mut_path, pos_mut_region))
        
        # 📎 Annotate Sites with Genomic Features
        Pos_Mut_sites_GR <- annotate_exons_introns_UTRs(Pos_Mut_sites_GR, exons_GR, introns_GR, utrs_GR, gtf_GR)
        Pos_Mut_regions_GR <- annotate_exons_introns_UTRs(Pos_Mut_regions_GR, exons_GR, introns_GR, utrs_GR, gtf_GR)
        
        Pos_Mut_sites_GR <- annotate_gene_type(Pos_Mut_sites_GR, gtf_GR)
        Pos_Mut_regions_GR <- annotate_gene_type(Pos_Mut_regions_GR, gtf_GR)
        
        # 📊 Aggregate Counts
        sites_annot_type_counts <- table(Pos_Mut_sites_GR$annotation)
        regions_annot_type_counts <- table(Pos_Mut_regions_GR$annotation)
        
        sites_gene_type_counts <- table(Pos_Mut_sites_GR$gene_type)
        regions_gene_type_counts <- table(Pos_Mut_regions_GR$gene_type)
        
        Conditions[ii] <- paste0("cov_", Cov_thresh, "_diff_", Diff_thresh)
        
        # 🧾 Build Annotated Data Frames
        df_sites_annot_type <- data.frame(
          Annot_type = names(sites_annot_type_counts),
          N = as.integer(sites_annot_type_counts),
          Frac_sites = prop.table(sites_annot_type_counts),
          Cov_thresh = Cov_thresh,
          Diff_thresh = Diff_thresh
        )
        
        df_regions_annot_type <- data.frame(
          Annot_type = names(regions_annot_type_counts),
          N = as.integer(regions_annot_type_counts),
          Frac_sites = prop.table(regions_annot_type_counts),
          Cov_thresh = Cov_thresh,
          Diff_thresh = Diff_thresh
        )
        
        # 📚 Combine for Bar Plotting
        if (ii == 1) {
          DF_sites_annot_type <- df_sites_annot_type
          DF_regions_annot_type <- df_regions_annot_type
        } else {
          DF_sites_annot_type <- rbind(DF_sites_annot_type, df_sites_annot_type)
          DF_regions_annot_type <- rbind(DF_regions_annot_type, df_regions_annot_type)
        }
        
        ii <- ii + 1
      }
    }
    
    # 🎨 Generate Bar Plots
    Res_path <- paste0(Mut_paths[i], filt_by, "/")
    make_bar_plot(df = DF_sites_annot_type, type = "fraction", Color,
                  Cov_thresh = 50, Diff_thresh = 0.03, x_lim = c(), x_tick_interval = 0.1,
                  x_lab = "Fraction of sites", y_lab = "Annotation Type",
                  fig_name = paste0(Res_path, "Pos_Mutsites_annotation_barplot.pdf"), 
                  fig_sz = c(6, 7))
    
    make_bar_plot(df = DF_regions_annot_type, type = "fraction", Color,
                  Cov_thresh = 50, Diff_thresh = 0.03, x_lim = c(), x_tick_interval = 0.1,
                  x_lab = "Fraction of regions", y_lab = "Annotation Type",
                  fig_name = paste0(Res_path, "Pos_Mutregions_annotation_barplot.pdf"), 
                  fig_sz = c(6, 7))
  }
}

###########################################################################################################

#  Suggested Caption for Supplementary Figure

# Bar plots were generated using threshold-filtered mutation annotations (Cov ≥ 50, Δ ≥ 0.03). Values represent either absolute counts or fractions per genomic category (e.g., exon, intron, UTR). Custom color mappings and tick intervals were applied for visual clarity.


# 🧾 Suggested Figure Legend

# Fractional distribution of mutation sites and mutation regions across annotated transcript features (Exon, Intron, 3'UTR, 5'UTR, Intergenic). Each bar summarizes mutation enrichment for a given coverage and difference threshold (Cov ≥ 50, Δ ≥ 0.03).
