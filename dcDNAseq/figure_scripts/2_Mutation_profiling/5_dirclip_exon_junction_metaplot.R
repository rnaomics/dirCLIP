#!/usr/bin/env Rscript
# Distance of dirCLIP Mutation Sites from Exon Junctions
# Manuscript-specific figure-generation script.
#
# Purpose:
# Calculate strand-aware distances between significant dirCLIP
# mutation sites and exon donor/acceptor junctions and generate
# metaplots of their positional distributions.

############################################################
# Distance of dirCLIP mutation sites from exon donor
# and acceptor sites
############################################################

##############################
# Set working directory
##############################
# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")

##############################
# Set library path
##############################
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

##############################
# Load required packages
##############################
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(patchwork)

############################################################
# Function: Calculate distances from exon donor or
# acceptor sites
############################################################
get_site_distance_from_exon_ends <- function(
    dirclip_gr,
    exon_gr,
    distance_type) {
  
  ##########################################################
  # Define donor or acceptor positions in a strand-aware
  # manner
  ##########################################################
  if (distance_type == "dist_from_acceptor") {
    
    exon_boundary_pos <- ifelse(
      strand(exon_gr) == "+",
      start(exon_gr),
      end(exon_gr)
    )
    
  } else if (distance_type == "dist_from_donor") {
    
    exon_boundary_pos <- ifelse(
      strand(exon_gr) == "+",
      end(exon_gr),
      start(exon_gr)
    )
  }
  
  ##########################################################
  # Convert exon boundaries into single-base GRanges
  ##########################################################
  exon_boundary_gr <- GRanges(
    seqnames = seqnames(exon_gr),
    ranges = IRanges(
      start = exon_boundary_pos,
      end = exon_boundary_pos
    ),
    strand = strand(exon_gr)
  )
  
  exon_boundary_gr <- unique(exon_boundary_gr)
  
  ##########################################################
  # Identify nearest exon boundary for each dirCLIP site
  ##########################################################
  nearest_hits <- distanceToNearest(
    dirclip_gr,
    exon_boundary_gr
  )
  
  ##########################################################
  # Extract matching ranges
  ##########################################################
  query_idx <- queryHits(nearest_hits)
  subject_idx <- subjectHits(nearest_hits)
  
  dirclip_hit_gr <- dirclip_gr[query_idx]
  exon_hit_gr <- exon_boundary_gr[subject_idx]
  
  ##########################################################
  # Extract coordinates and strand information
  ##########################################################
  dirclip_pos <- start(dirclip_hit_gr)
  
  exon_boundary_pos <- start(exon_hit_gr)
  
  exon_strand <- as.character(
    strand(exon_hit_gr)
  )
  
  ##########################################################
  # Compute strand-aware distance
  ##########################################################
  distance <- ifelse(
    exon_strand == "+",
    dirclip_pos - exon_boundary_pos,
    exon_boundary_pos - dirclip_pos
  )
  
  ##########################################################
  # Build results table
  ##########################################################
  distance_df <- data.frame(
    dirclip_chr = as.character(
      seqnames(dirclip_hit_gr)
    ),
    
    dirclip_pos = dirclip_pos,
    
    exon_chr = as.character(
      seqnames(exon_hit_gr)
    ),
    
    exon_start = exon_boundary_pos,
    
    exon_strand = exon_strand,
    
    distance = distance,
    
    score = mcols(dirclip_hit_gr)$score
  )
  
  ##########################################################
  # Summary statistics
  ##########################################################
  summary(distance_df$distance)
  table(distance_df$exon_strand)
  
  return(distance_df)
}

############################################################
# Function: Apply plotting style
############################################################
add_plot_theme <- function(
    figure,
    text_size) {
  
  figure <- figure +
    theme_minimal(text_size) +
    theme(
      text = element_text(size = text_size),
      
      panel.grid = element_blank(),
      
      panel.border = element_blank(),
      
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      
      axis.text.x = element_text(
        size = text_size,
        color = "black",
        angle = 0,
        hjust = 0.5
      ),
      
      axis.title.x = element_text(
        size = text_size
      ),
      
      axis.text.y = element_text(
        size = text_size,
        color = "black"
      ),
      
      axis.title.y = element_text(
        size = text_size
      ),
      
      plot.title = element_text(
        size = text_size,
        hjust = 0.5
      ),
      
      axis.ticks = element_line(
        color = "black",
        size = 0.8
      ),
      
      axis.ticks.length =
        unit(0.25, "cm"),
      
      axis.ticks.x =
        element_line(color = "black"),
      
      axis.ticks.y =
        element_line(color = "black")
    )
  
  return(figure)
}

############################################################
# Function: Set y-axis limits and tick marks
############################################################
set_y_axis_ticks <- function(
    figure,
    y_tick_interval) {
  
  plot_build <- ggplot_build(figure)
  
  density_df <- plot_build$data[[1]]
  
  max_density <- max(
    density_df$density
  )
  
  y_min <- 0
  
  y_max <- ceiling(
    max_density / y_tick_interval
  ) * y_tick_interval
  
  y_ticks <- seq(
    y_min,
    y_max,
    by = y_tick_interval
  )
  
  figure <- figure +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = expansion(
        mult = c(0, 0.02)
      )
    )
  
  return(figure)
}

############################################################
# Function: Generate density plot
############################################################
make_density_plot <- function(
    distance_df,
    fill_color,
    x_label,
    y_label) {
  
  text_size <- 16
  
  density_plot <- ggplot(
    distance_df,
    aes(x = distance)
  ) +
    geom_density(
      linewidth = 1,
      adjust = 1,
      fill = fill_color,
      alpha = 0.5
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "grey40"
    ) +
    labs(
      x = x_label,
      y = y_label,
      title = ""
    )
  
  density_plot <- add_plot_theme(
    density_plot,
    text_size
  )
  
  density_plot <- set_y_axis_ticks(
    density_plot,
    y_tick_interval = 0.001
  )
  
  return(density_plot)
}

############################################################
# Function: Plot distances from donor and acceptor sites
############################################################
plot_distance_from_exon_ends <- function(
    donor_distance_df,
    acceptor_distance_df,
    distance_threshold,
    fill_color,
    figure_size,
    figure_name) {
  
  ##########################################################
  # Restrict distances to plotting window
  ##########################################################
  donor_distance_filtered <-
    donor_distance_df[
      abs(donor_distance_df$distance) <=
        distance_threshold,
    ]
  
  acceptor_distance_filtered <-
    acceptor_distance_df[
      abs(acceptor_distance_df$distance) <=
        distance_threshold,
    ]
  
  ##########################################################
  # Generate donor and acceptor density plots
  ##########################################################
  donor_plot <- make_density_plot(
    distance_df =
      donor_distance_filtered,
    fill_color = fill_color,
    x_label =
      "Distance from exon donor site (bp)",
    y_label = "Density"
  )
  
  acceptor_plot <- make_density_plot(
    distance_df =
      acceptor_distance_filtered,
    fill_color = fill_color,
    x_label =
      "Distance from exon acceptor site (bp)",
    y_label = "Density"
  )
  
  ##########################################################
  # Combine panels
  ##########################################################
  combined_figure <- donor_plot + acceptor_plot
  
  ##########################################################
  # Save figure
  ##########################################################
  ggsave(
    figure_name,
    plot = combined_figure,
    width = figure_size[1] * 2,
    height = figure_size[2],
    limitsize = FALSE,
    device = cairo_pdf
  )
}

############################################################
# Input files: SRSF3
############################################################
dirclip_file <-
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/cov_thresh_50/Pos_mut_thresh_0.03.tsv"

annotation_gtf <-
  "genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz"

figure_name <-
  "Result/figure/1_Mutation_profiling/SRSF3_mutsites_distance_from_exon_ends_metaplot.pdf"

fill_color <- "#66c2a5"

############################################################
# Alternative datasets
############################################################
# Hmga1a
# dirclip_file <- "..."
# figure_name <- "..."
# fill_color <- "#E26A6D"

# Hmga1b
# dirclip_file <- "..."
# figure_name <- "..."
# fill_color <- "#FC9272"

############################################################
# Load gene annotation
############################################################
gtf_gr <- import(annotation_gtf)

exon_gr <- gtf_gr[
  gtf_gr$type == "exon"
]

############################################################
# Load dirCLIP mutation sites
############################################################
dirclip_df <- read.table(
  dirclip_file,
  col.names = c(
    "chr",
    "pos",
    "score"
  )
)

dirclip_gr <- GRanges(
  seqnames = dirclip_df$chr,
  ranges = IRanges(
    start = dirclip_df$pos,
    end = dirclip_df$pos
  ),
  score = dirclip_df$score
)

############################################################
# Calculate strand-aware distances
############################################################
donor_distance_df <-
  get_site_distance_from_exon_ends(
    dirclip_gr,
    exon_gr,
    distance_type =
      "dist_from_donor"
  )

acceptor_distance_df <-
  get_site_distance_from_exon_ends(
    dirclip_gr,
    exon_gr,
    distance_type =
      "dist_from_acceptor"
  )

############################################################
# Generate metaplots
############################################################
plot_distance_from_exon_ends(
  donor_distance_df,
  acceptor_distance_df,
  distance_threshold = 200,
  fill_color = fill_color,
  figure_size = c(6, 7),
  figure_name = figure_name
)

############################################################
# End of script
############################################################