#!/usr/bin/env Rscript
# Mutation Site Distribution Across Transcript Length
# Manuscript-specific figure-generation script.
#
# Purpose:
# Calculate transcript-length normalized mutation-site positions
# and generate density plots illustrating their distribution.
#
# Preserve the original implementation and parameters.

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

library(rtracklayer)  # import()
library(ggplot2)      # ggplot(), geom_density(), theme_*, scale_*, ggsave()
library(cowplot)      # ggdraw(), draw_plot()
library(dplyr)        # group_by(), summarise(), ungroup()
library(grid)         # unit()

########################## Mutation-site distribution plot ##################################################

mutsites_distribution_plot = function(
    dirclip_file,
    tx_lengths,
    Color,
    fig_name,
    fig_sz
){
  
  text_size = 16
  
  Mut_sites <- read.delim(
    dirclip_file,
    header = FALSE
  )
  
  colnames(Mut_sites) =
    c(
      "Transcipt_id",
      "Loc",
      "Mut_rate_diff"
    )
  
  Mut_sites =
    Mut_sites[
      ,
      c(
        "Transcipt_id",
        "Loc"
      )
    ]
  
  ############################################################
  
  Mut_sites$Transcipt_id <-
    sub(
      "^Chr_",
      "",
      Mut_sites$Transcipt_id
    )
  
  ########### Merge mutation sites and transcript lengths ####################
  
  colnames(tx_lengths)[1] =
    "Transcipt_id"
  
  Mut_sites =
    merge(
      Mut_sites,
      tx_lengths,
      by = "Transcipt_id"
    )
  
  ########## Normalize by transcript length ##################################
  
  Mut_sites$norm_loc =
    Mut_sites$Loc /
    Mut_sites$transcript_length
  
  ######### Make density plot ###############################################
  
  Title = ""
  
  density_plot =
    make_density_plot(
      df = Mut_sites,
      Color,
      Title,
      text_size
    )
  
  ######### Save figure #####################################################
  
  ggsave(
    fig_name,
    plot = density_plot,
    width = fig_sz[1] * 2,
    height = fig_sz[2],
    limitsize = FALSE,
    device = cairo_pdf
  )
}


############ Function to add theme ######################################

add_theme = function(
    plot_object,
    text_size
){
  
  plot_object =
    plot_object +
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
  
  return(plot_object)
}

######## Function to set axis limits and ticks ##################################################

set_axis_lim_ticks = function(
    plot_object,
    df
){
  
  max_density =
    max(
      density(df$norm_loc)[["y"]]
    )
  
  y_ticks =
    pretty(
      c(0, max_density),
      n = 5
    )
  
  y_min = min(y_ticks)
  y_max = max(y_ticks)
  
  x_ticks = seq(
    0,
    1,
    by = 0.2
  )
  
  x_min = min(x_ticks)
  x_max = max(x_ticks)
  
  plot_object =
    plot_object +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = x_ticks,
      expand = expansion(mult = c(0, 0.02))
    )
  
  return(plot_object)
}


########## Function to make density plot ###########################################

make_density_plot = function(
    df,
    Color,
    Title,
    text_size
){
  
  y_lab = "Density"
  x_lab = "Mut sites locs normalised by transcript length"
  
  alpha = 0.8
  
  ################## 1. Density plot with filled area #####################
  
  filled_density_plot =
    ggplot(df, aes(x = norm_loc)) +
    geom_density(
      alpha = alpha,
      size = 1.0,
      color = "black",
      fill = Color
    ) +
    labs(
      y = y_lab,
      x = x_lab,
      title = Title
    )
  
  filled_density_plot =
    add_theme(
      filled_density_plot,
      text_size
    )
  
  filled_density_plot =
    set_axis_lim_ticks(
      filled_density_plot,
      df
    )
  
  ################## 2. Density plot with line only #####################
  
  line_density_plot =
    ggplot(df, aes(x = norm_loc)) +
    geom_density(
      alpha = alpha,
      size = 1.0,
      color = Color
    ) +
    labs(
      y = y_lab,
      x = x_lab,
      title = Title
    )
  
  line_density_plot =
    add_theme(
      line_density_plot,
      text_size
    )
  
  line_density_plot =
    set_axis_lim_ticks(
      line_density_plot,
      df
    )
  
  ##########################################################################
  
  combined_plot = ggdraw()
  
  combined_plot =
    combined_plot +
    draw_plot(
      filled_density_plot,
      x = 0,
      y = 0,
      width = 0.5,
      height = 1
    )
  
  combined_plot =
    combined_plot +
    draw_plot(
      line_density_plot,
      x = 0.5,
      y = 0,
      width = 0.5,
      height = 1
    )
  
  return(combined_plot)
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

#############################################################################################################################

dirclip_file =
  "figure_scripts/figure_data/1_Mutation_profiling/2_Isoform_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/2_Inc_del/cov_thresh_50/Pos_mut_thresh_0.03.tsv"

annotation_gtf =
  "genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz"

figure_name =
  "Result/figure/1_Mutation_profiling/SRSF3_Mutsites_distribution_across_length.pdf"

Color = "#66c2a5"

fig_sz = c(7, 7)

#############################################################
gtf_granges =
  import(annotation_gtf)

################# Get transcript lengths ####################

transcript_lengths_df =
  get_transcript_length_from_gtf(
    gtf_granges
  )

#############################################################

########## Make mutation site distribution plot ################################

mutsites_distribution_plot(
  dirclip_file,
  transcript_lengths_df,
  Color,
  fig_name = figure_name,
  fig_sz = fig_sz
)

###########################################################