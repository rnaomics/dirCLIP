#!/usr/bin/env Rscript
# Mutation Sites Per Read
# Manuscript-specific figure-generation script.
# Quantifies significant mutation sites per read and per isoform.

############################################################
# Quantify dirCLIP mutation sites per read and per isoform
############################################################

# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

##############################
# Load required packages
##############################
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(grid)
library(patchwork)

############################################################
# Function: Count mutation sites per read
############################################################
count_mutation_sites_per_read <- function(
    read_mutation_locations,
    mutation_sites) {
  
  ##########################################################
  # Retain only reads mapped to transcripts that contain
  # identified dirCLIP mutation sites
  ##########################################################
  keep_idx <- read_mutation_locations$ref_name %in%
    mutation_sites$Transcript
  
  read_mutation_locations <-
    read_mutation_locations[keep_idx, ]
  
  ##########################################################
  # Count mutation sites per transcript
  ##########################################################
  transcript_mutation_counts <- mutation_sites %>%
    group_by(Transcript) %>%
    summarise(
      N_mut = n(),
      .groups = "drop"
    )
  
  ##########################################################
  # Expand comma-separated mutation positions
  ##########################################################
  read_positions_long <- read_mutation_locations %>%
    mutate(
      ref_mod_pos_clean =
        str_replace_all(
          ref_mod_pos,
          "\\s+",
          ""
        )
    ) %>%
    separate_rows(
      ref_mod_pos_clean,
      sep = ",",
      convert = TRUE
    ) %>%
    rename(pos = ref_mod_pos_clean)
  
  ##########################################################
  # Retain only mutation positions matching the
  # mutation catalog for the corresponding transcript
  ##########################################################
  matching_positions <- read_positions_long %>%
    inner_join(
      mutation_sites %>%
        select(
          Transcript,
          Loc
        ),
      by = c(
        "ref_name" = "Transcript",
        "pos" = "Loc"
      )
    )
  
  ##########################################################
  # Summarize mutation sites per read
  ##########################################################
  read_mutation_summary <- matching_positions %>%
    group_by(
      read_id,
      ref_name
    ) %>%
    summarise(
      ref_mut_pos =
        paste(
          sort(unique(pos)),
          collapse = ","
        ),
      
      N_mut_per_read =
        n_distinct(pos),
      
      .groups = "drop"
    )
  
  ##########################################################
  # Add transcript-level mutation counts
  ##########################################################
  read_mutation_summary <-
    read_mutation_summary %>%
    left_join(
      transcript_mutation_counts,
      by = c(
        "ref_name" = "Transcript"
      )
    )
  
  return(read_mutation_summary)
}

############################################################
# Function: Generate mutation-site summaries per read
############################################################
generate_mutation_sites_per_read <- function(
    read_mutation_paths,
    mutation_site_file,
    mutation_location_file) {
  
  for (read_mutation_path in read_mutation_paths) {
    
    print(read_mutation_path)
    
    ########################################################
    # Load read-level mutation locations
    ########################################################
    read_mutation_locations <- read.delim2(
      paste0(
        read_mutation_path,
        mutation_location_file
      )
    )
    
    ########################################################
    # Load mutation catalog
    ########################################################
    mutation_sites <- read.delim(
      mutation_site_file,
      header = FALSE
    )
    
    colnames(mutation_sites) <- c(
      "Transcript",
      "Loc",
      "mut_rate_diff"
    )
    
    ########################################################
    # Count mutation sites per read
    ########################################################
    mutation_sites_per_read <-
      count_mutation_sites_per_read(
        read_mutation_locations,
        mutation_sites
      )
    
    ########################################################
    # Save output as gzipped TSV
    ########################################################
    output_connection <- gzfile(
      paste0(
        read_mutation_path,
        "df_sites_per_read.tsv.gz"
      ),
      open = "wt"
    )
    
    write.table(
      mutation_sites_per_read,
      file = output_connection,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE,
      col.names = TRUE
    )
    
    close(output_connection)
  }
}

############################################################
# Function: Apply publication-ready plotting theme
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
# Function: Mutation sites per read and isoform
############################################################
make_mutation_sites_per_read_plot <- function(
    read_mutation_path,
    figure_name,
    figure_size) {
  
  ##########################################################
  # Load mutation summary table
  ##########################################################
  mutation_sites_per_read <- read.delim(
    paste0(
      read_mutation_path,
      "df_sites_per_read.tsv.gz"
    )
  )
  
  ##########################################################
  # Categorize reads based on the number of dirCLIP sites
  ##########################################################
  mutation_sites_per_read <-
    mutation_sites_per_read %>%
    mutate(
      Nmut_bin = case_when(
        N_mut_per_read >= 1 &
          N_mut_per_read <= 2 ~ "1–2",
        
        N_mut_per_read >= 3 &
          N_mut_per_read <= 4 ~ "3–4",
        
        N_mut_per_read >= 5 &
          N_mut_per_read <= 6 ~ "5–6",
        
        N_mut_per_read >= 7 &
          N_mut_per_read <= 8 ~ "7–8",
        
        N_mut_per_read >= 9 &
          N_mut_per_read <= 10 ~ "9–10",
        
        N_mut_per_read > 10 ~ ">10",
        
        TRUE ~ NA_character_
      ),
      
      Nmut_bin = factor(
        Nmut_bin,
        levels = c(
          "1–2",
          "3–4",
          "5–6",
          "7–8",
          "9–10",
          ">10"
        )
      )
    )
  
  ##########################################################
  # Plot 1: Mutation sites per read
  ##########################################################
  text_size <- 16
  
  read_plot <- ggplot(
    mutation_sites_per_read,
    aes(x = Nmut_bin)
  ) +
    geom_bar(
      fill = "#2C7FB8",
      color = "white",
      linewidth = 0.3
    ) +
    labs(
      title = "",
      x = "Number of dirCLIP sites per read",
      y = "Number of reads"
    )
  
  read_plot <- add_plot_theme(
    figure = read_plot,
    text_size = text_size
  )
  
  ##########################################################
  # Plot 2: Mutation sites per isoform
  ##########################################################
  mutation_sites_per_isoform <- unique(
    mutation_sites_per_read[
      ,
      c(
        "ref_name",
        "N_mut"
      )
    ]
  )
  
  isoform_plot <- ggplot(
    mutation_sites_per_isoform,
    aes(x = N_mut)
  ) +
    geom_histogram(
      binwidth = 10,
      color = "white",
      fill = "#2C7FB8"
    ) +
    labs(
      title = "",
      x = "Number of dirCLIP sites per isoform",
      y = "Number of isoforms"
    )
  
  isoform_plot <- add_plot_theme(
    figure = isoform_plot,
    text_size = text_size
  )
  
  ##########################################################
  # Combine both panels
  ##########################################################
  combined_figure <- ggdraw()
  
  combined_figure <- combined_figure +
    draw_plot(
      read_plot,
      x = 0,
      y = 0,
      width = 0.5,
      height = 1
    )
  
  combined_figure <- combined_figure +
    draw_plot(
      isoform_plot,
      x = 0.5,
      y = 0,
      width = 0.5,
      height = 1
    )
  
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
# Function: Set y-axis ticks
############################################################
set_yticks <- function(
    fig,
    y_ticks) {
  
  fig <- fig +
    scale_y_continuous(
      limits = c(
        min(y_ticks),
        max(y_ticks)
      ),
      breaks = y_ticks,
      expand = expansion(
        mult = c(0, 0.02)
      )
    )
  
  return(fig)
}

############################################################
# Function: Distance bin summaries
############################################################
get_bin_counts <- function(
    df,
    col_name,
    breaks,
    labels) {
  
  id <- colnames(df) %in% col_name
  colnames(df)[id] <- "distance"
  
  ##########################################################
  # Raw counts
  ##########################################################
  all_dist <- unlist(df$distance)
  
  dist_df <- data.frame(
    distance = all_dist
  ) %>%
    mutate(
      bin = cut(
        distance,
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE
      )
    )
  
  df_counts <- dist_df %>%
    count(bin)
  
  ##########################################################
  # Weighted counts
  #
  # Each read contributes total weight = 1
  ##########################################################
  weighted_df <- df %>%
    mutate(
      n_dist = lengths(distance),
      weight = 1 / n_dist
    ) %>%
    select(
      read_id,
      distance,
      weight
    ) %>%
    tidyr::unnest_longer(distance)
  
  df_weighted_counts <- weighted_df %>%
    mutate(
      bin = cut(
        distance,
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE
      )
    ) %>%
    group_by(bin) %>%
    summarise(
      n = sum(weight),
      .groups = "drop"
    )
  
  result <- list(
    df_counts,
    df_weighted_counts
  )
  
  names(result) <- c(
    "df_counts",
    "df_weighted_counts"
  )
  
  return(result)
}

############################################################
# Function: Bar plot
############################################################
get_bar_plot <- function(
    df,
    x_lab,
    y_lab) {
  
  text_size <- 16
  
  fig <- ggplot(
    df,
    aes(
      x = bin,
      y = n
    )
  ) +
    geom_col(
      fill = "#2C7FB8",
      color = "white",
      linewidth = 0.3
    ) +
    labs(
      title = "",
      x = x_lab,
      y = y_lab
    )
  
  fig <- add_plot_theme(
    fig,
    text_size
  )
  
  ymax <- max(fig$data$n)
  
  y_ticks <- pretty(
    c(0, ymax),
    n = 6
  )
  
  fig <- set_yticks(
    fig,
    y_ticks
  )
  
  return(fig)
}

############################################################
# Function: Distances between mutation sites per read
############################################################
plot_distance_between_mut_sites_per_read <- function(
    read_mutation_path,
    fig_name,
    fig_sz) {
  
  df_sites_per_read <- read.delim(
    paste0(
      read_mutation_path,
      "df_sites_per_read.tsv.gz"
    )
  )
  
  ##########################################################
  # Calculate mutation-site distances
  ##########################################################
  df <- df_sites_per_read %>%
    filter(
      N_mut_per_read > 1
    ) %>%
    mutate(
      pos = map(
        ref_mut_pos,
        ~ as.numeric(
          strsplit(.x, ",")[[1]]
        )
      ),
      
      # distance from first site
      dist_from_first =
        map(
          pos,
          ~ .x[-1] - .x[1]
        ),
      
      # distance from last site
      dist_from_last =
        map(
          pos,
          ~ .x[length(.x)] -
            .x[-length(.x)]
        ),
      
      # distance between consecutive sites
      dist_consecutive =
        map(
          pos,
          diff
        )
    ) %>%
    select(-pos)
  
  breaks <- c(
    0,
    25,
    100,
    200,
    500,
    Inf
  )
  
  labels <- c(
    "1-25",
    "26-100",
    "101-200",
    "201-500",
    ">500"
  )
  
  ##########################################################
  # First site
  ##########################################################
  first_df <- get_bin_counts(
    df,
    "dist_from_first",
    breaks,
    labels
  )
  
  f1_counts <- get_bar_plot(
    first_df[["df_counts"]],
    "Distance from the first site",
    "Counts"
  )
  
  f1_weighted <- get_bar_plot(
    first_df[["df_weighted_counts"]],
    "Distance from the first site",
    "Weighted Counts"
  )
  
  ##########################################################
  # Last site
  ##########################################################
  last_df <- get_bin_counts(
    df,
    "dist_from_last",
    breaks,
    labels
  )
  
  f2_counts <- get_bar_plot(
    last_df[["df_counts"]],
    "Distance from the last site",
    "Counts"
  )
  
  f2_weighted <- get_bar_plot(
    last_df[["df_weighted_counts"]],
    "Distance from the last site",
    "Weighted Counts"
  )
  
  ##########################################################
  # Consecutive sites
  ##########################################################
  consecutive_df <- get_bin_counts(
    df,
    "dist_consecutive",
    breaks,
    labels
  )
  
  f3_counts <- get_bar_plot(
    consecutive_df[["df_counts"]],
    "Distance between consecutive sites",
    "Counts"
  )
  
  f3_weighted <- get_bar_plot(
    consecutive_df[["df_weighted_counts"]],
    "Distance between consecutive sites",
    "Weighted Counts"
  )
  
  ##########################################################
  # Combine panels
  ##########################################################
  fig <- (f1_counts | f2_counts | f3_counts) /
    (f1_weighted | f2_weighted | f3_weighted)
  
  ggsave(
    fig_name,
    plot = fig,
    width = fig_sz[1] * 3,
    height = fig_sz[2] * 2,
    limitsize = FALSE,
    device = cairo_pdf
  )
}

############################################################
# Input files
############################################################
read_mutation_paths <- c(
  "figure_scripts/figure_data/4_Co-occurence_sites/SRSF3_DIRCLIP/SRSF3_dcDNA/"
)

mutation_location_file <- "filt_q8_mut_loc.tsv.gz"

mutation_site_file <- paste0(
  "figure_scripts/figure_data/1_Mutation_profiling/",
  "2_Isoform_Mapping/",
  "1_SRSF3_cDNA_Mg/",
  "ref_TotalRNA_cDNA_no_UV/",
  "2_Inc_del/",
  "cov_thresh_50/",
  "Pos_mut_thresh_0.03.tsv"
)

figure_directory <- paste0(
  "Result/figure/",
  "5_Co-occurence_sites/",
  "3_SRSF3_DIRCLIP/",
  "SRSF3_dcDNA/",
  "Figures/"
)

############################################################
# Generate mutation-site summaries per read
############################################################
generate_mutation_sites_per_read(
  read_mutation_paths,
  mutation_site_file,
  mutation_location_file
)

############################################################
# Create figure directory
############################################################
ifelse(
  !dir.exists(figure_directory),
  dir.create(
    figure_directory,
    recursive = TRUE
  ),
  FALSE
)

############################################################
# Figure 1: Mutation sites per read and isoform
############################################################
make_mutation_sites_per_read_plot(
  read_mutation_path =
    read_mutation_paths[1],
  
  figure_name = paste0(
    figure_directory,
    "mut_sites_per_read.pdf"
  ),
  
  figure_size = c(7, 7)
)

############################################################
# Figure 2: Distance between mutation sites per read
############################################################
plot_distance_between_mut_sites_per_read(
  read_mutation_path =
    read_mutation_paths[1],
  
  fig_name = paste0(
    figure_directory,
    "dist_mut_sites_per_read.pdf"
  ),
  
  fig_sz = c(7, 7)
)

############################################################
# End of script
############################################################
