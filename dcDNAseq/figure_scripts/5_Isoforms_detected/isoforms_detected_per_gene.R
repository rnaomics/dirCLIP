#!/usr/bin/env Rscript
# Isoforms Detected per Gene
# Manuscript-specific figure-generation script.
#
# Purpose:
# Quantify the number of detected isoforms per gene
# and compare isoform diversity between samples.

# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

### Loading library #######
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library("cowplot")
library(scales)
library(S4Vectors)
library(matrixStats)
library(dplyr)
library(tidyr)
library(ggbreak)
library(RColorBrewer)
library(patchwork)
library(rtracklayer)

################################## Function to get unique transcript counts #########################################

get_transcript_unique_counts = function(Paths, Samples){
  
  for (sample_idx in seq_along(Paths)) {
    
    sample_file_paths = Paths[[sample_idx]]
    
    transcript_counts_df =
      get_counts(sample_file_paths)
    
    if (sample_idx == 1){
      
      merged_counts_df =
        transcript_counts_df
      
    } else {
      
      merged_counts_df =
        merge(
          merged_counts_df,
          transcript_counts_df,
          by = c("TXNAME", "GENEID"),
          all = TRUE
        )
      
      merged_counts_df[
        is.na(merged_counts_df)
      ] <- 0
    }
  }
  
  colnames(
    merged_counts_df
  )[-(1:2)] = Samples
  
  return(merged_counts_df)
}


get_counts = function(Sample_paths){
  
  for (file_idx in seq_along(Sample_paths)) {
    
    count_file =
      Sample_paths[file_idx]
    
    transcript_counts_df <-
      read.delim(count_file)
    
    colnames(transcript_counts_df) =
      c(
        "TXNAME",
        "GENEID",
        "Counts"
      )
    
    if (file_idx == 1){
      
      merged_counts_df =
        transcript_counts_df
      
    } else {
      
      merged_counts_df <-
        merge(
          merged_counts_df,
          transcript_counts_df,
          by = c("TXNAME", "GENEID"),
          all = TRUE
        )
    }
  }
  
  ####################################################
  # Sum counts across replicates
  ####################################################
  
  if (file_idx > 1){
    
    merged_counts_df[
      is.na(merged_counts_df)
    ] <- 0
    
    merged_counts_df$Counts =
      rowSums(
        merged_counts_df[
          ,
          -(1:2)
        ]
      )
    
    merged_counts_df =
      merged_counts_df[
        ,
        c(
          "TXNAME",
          "GENEID",
          "Counts"
        )
      ]
  }
  
  return(merged_counts_df)
}

#################### Function to get the number of isoforms expressed per gene ######################

get_isoform_expressed_per_gene = function(
    df,
    subsample_percent,
    Read_count_thresh
){
  
  ####################################################
  # Convert counts into passing/failing threshold
  ####################################################
  
  count_columns =
    colnames(df)[-(1:2)]
  
  df[, count_columns] =
    (df[, count_columns] *
       subsample_percent / 100) >=
    Read_count_thresh
  
  ####################################################
  # Count expressed isoforms per gene
  ####################################################
  
  isoform_counts_df <- df %>%
    group_by(GENEID) %>%
    summarise(
      across(
        where(is.logical),
        ~ sum(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  ####################################################
  # Remove genes with zero expressed isoforms
  ####################################################
  
  expressed_gene_idx =
    rowSums(
      isoform_counts_df[, -1]
    ) > 0
  
  isoform_counts_df =
    isoform_counts_df[
      expressed_gene_idx,
    ]
  
  return(isoform_counts_df)
}

############ Function to count the number of isoforms per gene #####################

count_nisoforms_per_gene = function(df){
  
  method_id = 2
  
  #### Option 1 ####################################
  
  if (method_id == 1){
    
    isoform_count_df <- df %>%
      pivot_longer(
        -1,
        names_to = "Condition",
        values_to = "Value"
      ) %>%
      count(Condition, Value) %>%
      pivot_wider(
        names_from = Condition,
        values_from = n,
        values_fill = 0
      ) %>%
      rename(
        N_isoform = Value
      )
    
  } else if (method_id == 2){
    
    #### Option 2 ####################################
    
    sample_names =
      colnames(df)[-1]
    
    n_isoforms =
      sort(
        unique(
          unlist(df[, -1])
        )
      )
    
    isoform_matrix =
      df[, -1]
    
    isoform_count_df =
      data.frame(n_isoforms)
    
    for (sample_idx in seq_along(sample_names)) {
      
      isoform_frequency_df =
        data.frame(
          table(
            isoform_matrix[, sample_idx]
          )[as.character(n_isoforms)]
        )
      
      isoform_count_df =
        cbind(
          isoform_count_df,
          isoform_frequency_df$Freq
        )
    }
    
    isoform_count_df[
      is.na(isoform_count_df)
    ] = 0
    
    colnames(isoform_count_df) =
      c(
        "N_isoform",
        sample_names
      )
  }
  
  #### Remove N_isoform = 0 #########################
  
  isoform_count_df =
    isoform_count_df[
      isoform_count_df$N_isoform > 0,
    ]
  
  return(isoform_count_df)
}

############# Function to make bar plot of number of isoforms expressed per gene ################

make_bar_plot = function(
    df,
    Colors,
    x_lab,
    y_lab,
    legend_title,
    Title,
    fig_name,
    fig_sz
){
  
  plot_df <- df %>%
    pivot_longer(
      cols = -N_isoform,
      names_to = "Condition",
      values_to = "Count"
    )
  
  text_size = 16
  
  bar_plot =
    ggplot(
      plot_df,
      aes(
        x = as.factor(N_isoform),
        y = Count,
        fill = Condition
      )
    ) +
    
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    ) +
    
    scale_fill_manual(
      values = Colors
    ) +
    
    labs(
      x = x_lab,
      y = y_lab,
      fill = legend_title,
      title = Title
    )
  
  bar_plot =
    add_theme(
      bar_plot,
      text_size
    )
  
  bar_plot =
    add_log_scale(
      bar_plot,
      df
    )
  
  ggsave(
    fig_name,
    plot = bar_plot,
    width = fig_sz[1] * 2,
    height = fig_sz[2],
    limitsize = FALSE
  )
  
  return(bar_plot)
}

######### Function to add theme ####################################

add_theme = function(
    plot_object,
    text_size
){
  
  plot_object =
    plot_object +
    theme_minimal(
      base_size = text_size
    ) +
    theme(
      panel.border = element_blank(),
      
      axis.line.x = element_line(
        color = "black",
        linewidth = 0.8
      ),
      
      axis.line.y = element_line(
        color = "black",
        linewidth = 0.8
      ),
      
      panel.grid = element_blank(),
      
      axis.text = element_text(
        size = text_size,
        color = "black"
      ),
      
      axis.title = element_text(
        size = text_size
      ),
      
      plot.title = element_text(
        hjust = 0.7
      ),
      
      axis.ticks = element_line(
        color = "black",
        size = 0.8
      ),
      
      axis.ticks.length = unit(
        0.25,
        "cm"
      ),
      
      axis.ticks.x = element_line(
        color = "black"
      ),
      
      axis.ticks.y = element_line(
        color = "black"
      ),
      
      axis.ticks.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top")
    )
  
  return(plot_object)
}

##### Function to create figure with log scale ####################################

add_log_scale = function(
    plot_object,
    df
){
  
  log_scale_plot =
    plot_object +
    scale_y_log10(
      limits = c(1, 10000),
      breaks = c(
        1,
        10,
        100,
        1000,
        10000
      ),
      labels = scales::comma
    ) +
    annotation_logticks(
      sides = "l"
    )
  
  combined_plot =
    plot_object | log_scale_plot
  
  return(combined_plot)
}

############# Function to make bar plot of number of isoforms expressed per gene ################

make_bar_plot2 = function(
    df,
    Colors,
    x_lab,
    y_lab,
    Title,
    legend_title,
    fig_name,
    fig_sz
){
  
  # Reshape to long format
  # Create the bins you want as an ordered factor
  breaks <- c(0, 1, 2, 3, 4, Inf)
  labels <- c("1", "2", "3", "4", "≥5")
  
  df <- df %>%
    mutate(
      Nmut_bin = cut(
        N_isoform,
        breaks = breaks,
        labels = labels,
        include.lowest = TRUE
      )
    )
  
  df$N_isoform = df$Nmut_bin
  df <- df %>%
    dplyr::select(-Nmut_bin)
  
  # Reshape to long format
  plot_df <- df %>%
    pivot_longer(
      cols = -N_isoform,
      names_to = "Condition",
      values_to = "Count"
    )
  
  plot_df <- plot_df %>%
    group_by(N_isoform, Condition) %>%
    summarise(
      Count = sum(Count),
      .groups = "drop"
    )
  
  # Ensure order of N_isoform
  
  plot_df$N_isoform <- factor(
    plot_df$N_isoform,
    levels = labels
  )
  
  # Plot
  
  text_size = 16
  
  bar_plot = ggplot(
    plot_df,
    aes(
      x = as.factor(N_isoform),
      y = Count,
      fill = Condition
    )
  ) +
    
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      linewidth = 0.6,
      alpha = 0.8
    ) +
    
    scale_fill_manual(
      values = c(Colors)
    ) +
    
    labs(
      x = x_lab,
      y = y_lab,
      fill = legend_title,
      title = Title
    )
  
  bar_plot = add_theme(
    bar_plot,
    text_size
  )
  
  bar_plot = add_log_scale(
    bar_plot,
    df
  )
  
  ################################################################
  
  ggsave(
    fig_name,
    plot = bar_plot,
    width = fig_sz[1] * 2,
    height = fig_sz[2],
    limitsize = FALSE
  )
  
  ##############################################
  
  return(bar_plot)
}

################ Getting the number of genes/isoforms above certain thresholds ###################

gene_above_threshold = function(raw_counts, Count_thresholds){
  
  counts_above_threshold = c()
  
  for (threshold_idx in seq_along(Count_thresholds)) {
    
    count_threshold =
      Count_thresholds[threshold_idx]
    
    count_matrix =
      as.matrix(
        raw_counts[, Samples]
      )
    
    counts_above_threshold =
      rbind(
        counts_above_threshold,
        colSums(
          count_matrix >= count_threshold
        )
      )
  }
  
  rownames(counts_above_threshold) =
    Count_thresholds
  
  return(counts_above_threshold)
}

################### Equate library sizes ################

equating_lib_size = function(df){
  
  library_sizes =
    colSums(
      df[, -(1:2)]
    )
  
  scaling_factors =
    library_sizes /
    min(library_sizes)
  
  normalized_counts_df = df
  
  for (sample_name in names(scaling_factors)) {
    
    normalized_counts_df[, sample_name] =
      normalized_counts_df[, sample_name] /
      scaling_factors[sample_name]
  }
  
  return(normalized_counts_df)
}

################### Get number of expressed transcripts ################

get_expressed_transcripts = function(
    df,
    subsample_percent,
    Read_count_thresholds
){
  
  subsampled_counts_df = df
  
  subsampled_counts_df[, -(1:2)] =
    subsampled_counts_df[, -(1:2)] *
    subsample_percent / 100
  
  expressed_transcript_counts =
    gene_above_threshold(
      raw_counts = subsampled_counts_df,
      Read_count_thresholds
    )
  
  return(expressed_transcript_counts)
}

#################### Get transcripts with mutation sites ####################

get_transcripts_with_mut_sites = function(Mut_sites_path){
  
  mutation_sites_df <- read.delim(
    Mut_sites_path,
    header = FALSE
  )
  
  transcript_names =
    unique(mutation_sites_df$V1)
  
  transcript_names <-
    sub(
      "^Chr_",
      "",
      transcript_names
    )
  
  return(transcript_names)
}

################## Reset transcripts without mutation sites ##################

reset_tx_without_mut_sites = function(
    df,
    Sample,
    transcript_names
){
  
  missing_transcript_idx =
    !df$TXNAME %in%
    transcript_names
  
  df[, Sample][missing_transcript_idx] = 0
  
  return(df)
}

################# Function to get transcript lengths #####################################
get_transcript_length_from_gtf = function(gtf_GR){
  
  # Keep only exon entries
  
  exon_granges <- gtf_GR[
    gtf_GR$type == "exon"
  ]
  
  exon_df =
    data.frame(exon_granges)
  
  transcript_lengths_df <- exon_df %>%
    group_by(transcript_id) %>%
    summarise(
      transcript_length =
        sum(width, na.rm = TRUE)
    ) %>%
    ungroup()
  
  return(transcript_lengths_df)
}

################### 2. equating the library size based on the sample with the least library size ################
get_tpm = function(df, transcript_lengths_df){
  
  # Scale_size = lib_size/min(lib_size)
  count_cols <- colnames(df)[-(1:2)]
  
  colnames(transcript_lengths_df)[1] = "TXNAME"
  df_norm = merge(df, transcript_lengths_df, by = c("TXNAME"))
  
  # lib_size = colSums(df[,-(1:2)])
  for (Col in count_cols){
    rpk <- df_norm[[Col]] / (df_norm$transcript_length / 1000)
    tpm <- rpk / sum(rpk, na.rm = TRUE) * 1e6
    
    df_norm[,Col] = tpm 
  }
  
  df_norm = df_norm[, !colnames(df_norm) %in% "transcript_length"]
  return(df_norm)
}
###################### Paths #######################################################################

sample_paths = list(
  c(
    "figure_scripts/figure_data/3_Isoform_counts/SRSF3_dirCLIP/rep1/gencode/filt_q8_read_counts_per_transcript.tsv",
    "figure_scripts/figure_data/3_Isoform_counts/SRSF3_dirCLIP/rep2/gencode/filt_q8_read_counts_per_transcript.tsv"
  ),
  c(
    "figure_scripts/figure_data/3_Isoform_counts/noUV_transcriptome/rep1/gencode/filt_q8_read_counts_per_transcript.tsv"
  )
)

Samples = c(
  "SRSF3 CLIP",
  "TotalRNA cDNA no UV"
)

Colors = c(
  "SRSF3 CLIP" = "#66c2a5",
  "TotalRNA cDNA no UV" = "#08519C"
)

Mut_sites_path =
  "figure_scripts/figure_data/1_Mutation_profiling/2_Isoform_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/cov_thresh_50/Pos_mut_thresh_0.03.tsv"

Result_path =
  "Result/figure/8_Isoforms_detected_per_gene/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/gencode/"

annotation_gtf = "genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz"
########################################################################################################

# Scale_percents = c(70, 100)
Scale_percents = c(100)

# Read_count_thresholds = c(5, 50)
Read_count_thresholds = c(5)
tpm_thresh = 0.5
################## Load mutation sites ###########################################################

transcript_names =
  get_transcripts_with_mut_sites(
    Mut_sites_path
  )

################### 1. Get transcript unique read counts ###############################################

transcript_counts_df =
  get_transcript_unique_counts(
    sample_paths,
    Samples
  )

# ################### 2. Equate library size ##############################################################
# 
# normalized_counts_df =
#   equating_lib_size(
#     transcript_counts_df
#   )

################### 2. calculate cpm ###############
gtf_GR = import(annotation_gtf)

transcript_lengths_df =
  get_transcript_length_from_gtf(
    gtf_GR
  )
# df_norm = get_cpm(df_counts)

normalized_counts_df = get_tpm(transcript_counts_df, transcript_lengths_df)

############################################################
# Reset transcripts with no mutation sites in the CLIP sample
############################################################

filtered_counts_df =
  reset_tx_without_mut_sites(
    df = transcript_counts_df,
    Sample = "SRSF3 CLIP",
    transcript_names
  )

filtered_normalized_counts_df =
  reset_tx_without_mut_sites(
    df = normalized_counts_df,
    Sample = "SRSF3 CLIP",
    transcript_names
  )

############################################################
# Retain only genes containing mutation-site-associated
# transcripts in the SRSF3 CLIP sample
############################################################
mut_genes =
  filtered_counts_df[
    filtered_counts_df$`SRSF3 CLIP` > 0,
  ]$GENEID

filtered_counts_df =
  filtered_counts_df[
    filtered_counts_df$GENEID %in% mut_genes,
  ]

filtered_normalized_counts_df =
  filtered_normalized_counts_df[
    filtered_normalized_counts_df$GENEID %in% mut_genes,
  ]

########################################################################################################

for (subsample_percent in Scale_percents){
  
  ################### Number of transcripts above threshold ###########################################
  
  expressed_isoform_counts =
    get_expressed_transcripts(
      transcript_counts_df,
      subsample_percent,
      Read_count_thresholds
    )
  
  normalized_expressed_isoform_counts =
    get_expressed_transcripts(
      normalized_counts_df,
      subsample_percent,
      tpm_thresh
    )
  
  ####################################################################################################
  
  for (Read_count_thresh in Read_count_thresholds){
    
    ################### Number of isoforms expressed per gene ########################################
    
    isoforms_per_gene_df =
      get_isoform_expressed_per_gene(
        filtered_counts_df,
        subsample_percent,
        Read_count_thresh
      )
    
    normalized_isoforms_per_gene_df =
      get_isoform_expressed_per_gene(
        filtered_normalized_counts_df,
        subsample_percent,
        tpm_thresh
      )
    
    ################### Count isoforms per gene #####################################################
    
    isoform_count_summary =
      count_nisoforms_per_gene(
        isoforms_per_gene_df
      )
    
    normalized_isoform_count_summary =
      count_nisoforms_per_gene(
        normalized_isoforms_per_gene_df
      )
    
    ################# Output path ###################################################################
    
    output_path = paste0(
      Result_path,
      "Scale_",
      subsample_percent,
      "/Read_thresh_",
      Read_count_thresh,
      "/"
    )
    
    ifelse(
      !dir.exists(output_path),
      dir.create(output_path, recursive = TRUE),
      FALSE
    )
    
    ################# Full distribution ############################################################
    
    make_bar_plot(
      df = isoform_count_summary,
      Colors,
      x_lab = "Number of Isoforms per Gene",
      y_lab = "Number of Genes",
      Title = paste0("Read thresh = ", Read_count_thresh),
      legend_title = "Sample",
      fig_name = paste0(
        output_path,
        "N_isoform_per_gene_barplot.pdf"
      ),
      fig_sz = c(7, 7)
    )
    
    make_bar_plot(
      df = normalized_isoform_count_summary,
      Colors,
      x_lab = "Number of Isoforms per Gene",
      y_lab = "Number of Genes",
      Title = paste0("TPM thresh = ", tpm_thresh),
      legend_title = "Sample",
      fig_name = paste0(
        output_path,
        "N_isoform_per_gene_barplot_norm.pdf"
      ),
      fig_sz = c(7, 7)
    )
    
    ####################################################################
    
    barplot_summary =
      make_bar_plot2(
        df = isoform_count_summary,
        Colors,
        x_lab = "Number of Isoforms per Gene",
        y_lab = "Number of Genes",
        Title = paste0("Read thresh = ", Read_count_thresh),
        legend_title = "Sample",
        fig_name = paste0(
          output_path,
          "more_than_5_isoforms_per_gene_barplot.pdf"
        ),
        fig_sz = c(7, 7)
      )
    
    normalized_barplot_summary =
      make_bar_plot2(
        df = normalized_isoform_count_summary,
        Colors,
        x_lab = "Number of Isoforms per Gene",
        y_lab = "Number of Genes",
        Title = paste0("TPM thresh = ", tpm_thresh),
        legend_title = "Sample",
        fig_name = paste0(
          output_path,
          "more_than_5_isoforms_per_gene_barplot_norm.pdf"
        ),
        fig_sz = c(7, 7)
      )
  }
}