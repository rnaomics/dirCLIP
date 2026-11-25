#!/usr/bin/env Rscript
# Mutation Site Annotation Analysis
# Manuscript-specific figure-generation script.
#
# Purpose:
# Annotate mutation sites with genomic features (UTRs, exons,
# introns, intergenic regions) and generate stacked-bar summaries.
#
### set work directory ##############
# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")

#### Set library path #########
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

#### Load required libraries ###############
library(rtracklayer)    
library(GenomicRanges)  
library(dplyr)          
library(ggplot2)        
library(reshape2)       
library(grid)    

###### Function to get exons and introns from GTF #############################

get_exons_introns_from_gtf = function(gtf_file, out_file){
  
  # Step 1: Load the GTF file
  
  gtf_granges <- import(gtf_file)
  
  ##### Remove random chromosomes ########################
  
  gtf_granges =
    gtf_granges[
      !grepl("_random", seqnames(gtf_granges))
    ]
  
  ########################################################
  
  transcript_granges =
    gtf_granges[
      gtf_granges$type == "transcript"
    ]
  
  reduced_transcripts <-
    reduce(
      split(
        transcript_granges,
        transcript_granges$transcript_id
      )
    )
  
  exon_granges =
    gtf_granges[
      gtf_granges$type == "exon"
    ]
  
  exons_by_tx =
    split(
      exon_granges,
      exon_granges$transcript_id
    )
  
  reduced_exons =
    reduce(exon_granges)
  
  ########################################################
  # Generate introns
  ########################################################
  
  intron_list = list()
  
  for (transcript_idx in seq_along(reduced_transcripts)) {
    
    print(transcript_idx)
    
    intron_list[transcript_idx] =
      GenomicRanges::setdiff(
        reduced_transcripts[[transcript_idx]],
        reduce(
          exons_by_tx[[transcript_idx]]
        )
      )
  }
  
  intron_granges <-
    do.call(c, intron_list)
  
  ########################################################
  # Remove introns overlapping exons
  ########################################################
  
  intron_granges =
    intron_granges[
      !overlapsAny(
        intron_granges,
        reduced_exons
      )
    ]
  
  reduced_introns =
    reduce(intron_granges)
  
  ########################################################
  # Convert GRanges objects to data frames
  ########################################################
  
  exons_reduced =
    data.frame(reduced_exons)
  
  introns_reduced =
    data.frame(reduced_introns)
  
  ########################################################
  # Save reduced exon and intron annotations
  ########################################################
  
  save(
    exons_reduced,
    introns_reduced,
    file = out_file
  )
}

#################### Function to generate 5'UTR and 3'UTR annotations ####################

get_UTR_from_gtf = function(gtf_file, out_file){
  
  gtf_granges <- import(gtf_file)
  
  gtf_granges =
    gtf_granges[
      !grepl(
        "_random",
        seqnames(gtf_granges)
      )
    ]
  
  ########################################################
  
  utrs =
    gtf_granges[
      gtf_granges$type == "UTR"
    ]
  
  start_codons =
    gtf_granges[
      gtf_granges$type == "start_codon"
    ]
  
  stop_codons =
    gtf_granges[
      gtf_granges$type == "stop_codon"
    ]
  
  ########################################################
  # Assign UTR type
  ########################################################
  
  utrs$UTR_type <- NA
  
  for (transcript_id in unique(utrs$transcript_id)) {
    
    transcript_utrs =
      utrs[
        utrs$transcript_id ==
          transcript_id
      ]
    
    transcript_strand =
      unique(
        as.character(
          strand(transcript_utrs)
        )
      )
    
    transcript_start_codon =
      start_codons[
        start_codons$transcript_id ==
          transcript_id
      ]
    
    transcript_stop_codon =
      stop_codons[
        stop_codons$transcript_id ==
          transcript_id
      ]
    
    if (
      length(transcript_start_codon) == 0 |
      length(transcript_stop_codon) == 0
    ) next
    
    ######################################################
    
    if (transcript_strand == "+") {
      
      start_pos =
        max(
          start(transcript_start_codon)
        )
      
      stop_pos =
        min(
          end(transcript_stop_codon)
        )
      
      transcript_utrs$UTR_type[
        start(transcript_utrs) < start_pos
      ] <- "5UTR"
      
      transcript_utrs$UTR_type[
        end(transcript_utrs) > stop_pos
      ] <- "3UTR"
      
    } else if (transcript_strand == "-") {
      
      start_pos =
        min(
          end(transcript_start_codon)
        )
      
      stop_pos =
        max(
          start(transcript_stop_codon)
        )
      
      transcript_utrs$UTR_type[
        end(transcript_utrs) > start_pos
      ] <- "5UTR"
      
      transcript_utrs$UTR_type[
        start(transcript_utrs) < stop_pos
      ] <- "3UTR"
    }
    
    ######################################################
    
    utrs[
      utrs$transcript_id ==
        transcript_id
    ]$UTR_type =
      transcript_utrs$UTR_type
  }
  
  save(
    utrs,
    file = out_file
  )
}

######### Function to load sites ####################

load_sites = function(site_file){
  
  ############ Load sites #######################
  
  site_df = read.delim(
    site_file,
    header = FALSE
  )
  
  ############ Convert to GRanges ###############
  
  site_granges <- GRanges(
    seqnames = site_df$V1,
    ranges = IRanges(
      start = site_df$V2,
      end = site_df$V2
    )
  )
  
  ############ Remove random chromosomes ########
  
  site_granges =
    site_granges[
      !grepl(
        "_random",
        seqnames(site_granges)
      )
    ]
  
  return(site_granges)
}

######### Function to annotate exons, introns and UTRs ############

annotate_exons_introns_UTRs = function(
    Sites_GR,
    exons_GR,
    introns_GR,
    utrs_GR,
    gtf_GR
){
  
  ######## Check overlaps with genes ############################################
  
  Sites_GR$annotation <- "Intergenic"
  
  gene_overlap_hits <- findOverlaps(
    Sites_GR,
    gtf_GR
  )
  
  Sites_GR$annotation[
    queryHits(gene_overlap_hits)
  ] <- "Genic"
  
  exon_overlap_hits <- findOverlaps(
    Sites_GR,
    exons_GR
  )
  
  Sites_GR$annotation[
    queryHits(exon_overlap_hits)
  ] <- "Exon"
  
  ### Assign remaining genic sites as introns ####
  
  Sites_GR[
    Sites_GR$annotation == "Genic"
  ]$annotation = "Intron"
  
  ### Annotate UTRs ################################
  
  UTR_3_GR =
    utrs_GR[
      utrs_GR$UTR_type == "3UTR"
    ]
  
  UTR_5_GR =
    utrs_GR[
      utrs_GR$UTR_type == "5UTR"
    ]
  
  utr_overlap_hits <- findOverlaps(
    Sites_GR,
    UTR_3_GR
  )
  
  Sites_GR$annotation[
    queryHits(utr_overlap_hits)
  ] <- "3UTR"
  
  utr_overlap_hits <- findOverlaps(
    Sites_GR,
    UTR_5_GR
  )
  
  Sites_GR$annotation[
    queryHits(utr_overlap_hits)
  ] <- "5UTR"
  
  return(Sites_GR)
}

######### Function to annotate gene type ############

annotate_gene_type = function(
    Sites_GR,
    gtf_GR
){
  
  # Find overlaps between sites and genes
  
  gene_overlap_hits <- findOverlaps(
    Sites_GR,
    gtf_GR
  )
  
  # Extract gene types
  
  matched_gene_types <- mcols(
    gtf_GR
  )$gene_type[
    subjectHits(gene_overlap_hits)
  ]
  
  # Assign gene types
  
  mcols(Sites_GR)$gene_type <-
    "Not_annotated"
  
  mcols(Sites_GR)$gene_type[
    queryHits(gene_overlap_hits)
  ] <- matched_gene_types
  
  return(Sites_GR)
}

############################# Make bar plot ##############################################

make_bar_plot = function(
    df,
    Colors,
    Title,
    y_lim,
    y_tick_interval,
    x_lab,
    y_lab,
    legend_title
){
  
  text_size = 16
  
  bar_plot =
    ggplot(
      df,
      aes(
        x = Group,
        y = y,
        fill = sub_group
      )
    ) +
    
    geom_bar(
      stat = "identity",
      width = 0.7,
      color = "black",
      linewidth = 0.6
    ) +
    
    labs(
      x = x_lab,
      y = y_lab,
      fill = legend_title,
      title = Title
    ) +
    
    scale_fill_manual(
      values = Colors
    ) +
    
    theme_minimal(
      base_size = text_size
    ) +
    
    theme(
      panel.grid = element_blank(),
      
      panel.border = element_blank(),
      
      axis.line = element_line(
        color = "black",
        linewidth = 0.8
      ),
      
      legend.position = "right",
      
      axis.title = element_text(
        size = text_size
      ),
      
      axis.text = element_text(
        size = text_size,
        color = "black"
      ),
      
      legend.title = element_text(
        size = text_size
      ),
      
      legend.text = element_text(
        size = text_size
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
  
  ##### Setting y-axis ticks ####################################
  
  if (length(y_lim) > 0){
    
    y_min = y_lim[1]
    y_max = y_lim[2]
    
  } else {
    
    group_summary_df <- df %>%
      group_by(Group) %>%
      summarise(
        total_y = sum(y),
        .groups = "drop"
      )
    
    y_min = 0
    
    y_max =
      ceiling(
        max(group_summary_df$total_y) /
          y_tick_interval
      ) *
      y_tick_interval
  }
  
  y_ticks =
    seq(
      y_min,
      y_max,
      by = y_tick_interval
    )
  
  bar_plot =
    bar_plot +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = expansion(
        mult = c(0, 0.02)
      )
    )
  
  return(bar_plot)
}

###################### Make stacked bar plot #######################################

Mut_sites_annotation_bar_plot = function(
    sites_annot,
    Colors,
    y_lim,
    y_tick_interval,
    fig_name,
    fig_sz
){
  
  annot_sites =
    sites_annot$n_sites_annot_type
  
  annot_sites =
    annot_sites /
    rowSums(annot_sites)
  
  ####################################################
  # Convert to long format
  ####################################################
  
  plot_df <- as.data.frame(annot_sites)
  
  plot_df$Group =
    rownames(plot_df)
  
  plot_df <- melt(
    plot_df,
    id.vars = "Group",
    variable.name = "sub_group",
    value.name = "y"
  )
  
  ####################################################
  # Create bar plot
  ####################################################
  
  Title = ""
  
  bar_plot =
    make_bar_plot(
      plot_df,
      Colors,
      Title,
      y_lim,
      y_tick_interval,
      x_lab = "Sample",
      y_lab = "Fraction of sites",
      legend_title = "Annotation"
    )
  
  ####################################################
  # Save figure
  ####################################################
  
  ggsave(
    fig_name,
    plot = bar_plot,
    width = fig_sz[1],
    height = fig_sz[2],
    limitsize = FALSE,
    device = cairo_pdf
  )
}

##################### Function to annotate the sites #####################################

annotate_sites = function(
    dirclip_files,
    exons_GR,
    introns_GR,
    utrs_GR,
    gtf_GR,
    annot_type,
    gene_type
){
  
  ##############################################################
  
  for (sample_idx in seq_along(dirclip_files)) {
    
    ######## Load mutation sites ##############################################
    
    sites_GR = load_sites(
      site_file = dirclip_files[sample_idx]
    )
    
    ####### Annotate exons, introns, 3'UTR and 5'UTR #########################
    
    sites_GR = annotate_exons_introns_UTRs(
      Sites_GR = sites_GR,
      exons_GR,
      introns_GR,
      utrs_GR,
      gtf_GR
    )
    
    ####### Annotate gene type ###############################################
    
    sites_GR = annotate_gene_type(
      Sites_GR = sites_GR,
      gtf_GR
    )
    
    ###### Count annotation types ############################################
    
    annotation_type_counts <- table(
      sites_GR$annotation
    )
    
    annotation_type_counts =
      ifelse(
        annot_type %in% names(annotation_type_counts),
        annotation_type_counts[annot_type],
        0
      )
    
    ###### Count gene types ##################################################
    
    gene_type_counts <- table(
      sites_GR$gene_type
    )
    
    gene_type_counts =
      ifelse(
        gene_type %in% names(gene_type_counts),
        gene_type_counts[gene_type],
        0
      )
    
    if (sample_idx == 1){
      
      n_sites_annot_type =
        annotation_type_counts
      
      n_sites_gene_type =
        gene_type_counts
      
    } else {
      
      n_sites_annot_type =
        rbind(
          n_sites_annot_type,
          annotation_type_counts
        )
      
      n_sites_gene_type =
        rbind(
          n_sites_gene_type,
          gene_type_counts
        )
    }
  }
  
  ##############################################################
  # Output
  ##############################################################
  
  colnames(n_sites_annot_type) =
    annot_type
  
  rownames(n_sites_annot_type) =
    Samples
  
  colnames(n_sites_gene_type) =
    gene_type
  
  rownames(n_sites_gene_type) =
    Samples
  
  annotation_summary =
    tibble(
      n_sites_annot_type,
      n_sites_gene_type
    )
  
  return(annotation_summary)
}

######################## Create folder ######################################################

create_folder = function(folder_path){
  
  ifelse(
    !dir.exists(folder_path),
    dir.create(
      folder_path,
      recursive = TRUE
    ),
    FALSE
  )
}

################################### Including deletions ###################################

dirclip_files = c(
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/3_Hmga1a_cDNA/ref_P19_WT_TotalRNA_cDNA_no_UV/2_Inc_del/cov_thresh_50/Pos_mut_thresh_0.03.tsv",
  "figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/4_Hmga1b_cDNA/ref_P19_WT_TotalRNA_cDNA_no_UV/2_Inc_del/cov_thresh_50/Pos_mut_thresh_0.03.tsv"
)

Samples = c(
  "Hmga1a",
  "Hmga1b"
)

figure_name =
  "Result/figure/1_Mutation_profiling/Hmga1a_1b_Mut_sites_annotation_stacked_bar.pdf"

annotation_source = "gencode"

gtf_file =
  "genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz"

exons_introns_RData =
  "genome/mm9/Annotation/gencode_exons_and_intron.RData"

UTR_RData =
  "genome/mm9/Annotation/gencode_UTR.RData"

############################################################
# Generate annotation resources
############################################################

get_exons_introns_from_gtf(
  gtf_file = gtf_file,
  out_file = exons_introns_RData
)

get_UTR_from_gtf(
  gtf_file = gtf_file,
  out_file = UTR_RData
)

############################################################
# Load exons and introns
############################################################

load(exons_introns_RData)
exons_GR = GRanges(exons_reduced)
introns_GR = GRanges(introns_reduced)

############################################################
# Load UTRs
############################################################

load(UTR_RData)

utrs_GR =
  utrs[
    !is.na(utrs$UTR_type)
  ]

############################################################
# Load genes
############################################################

gtf_GR <- import(gtf_file)

gtf_GR =
  gtf_GR[
    !grepl(
      "_random",
      seqnames(gtf_GR)
    )
  ]

gtf_GR =
  gtf_GR[
    gtf_GR$type == "gene"
  ]

############################################################

gene_type =
  c(
    unique(gtf_GR$gene_type),
    "Not_annotated"
  )

annot_type =
  c(
    "3UTR",
    "5UTR",
    "Exon",
    "Intergenic",
    "Intron"
  )

Colors = c(
  "3UTR" = "#66c2a5",
  "5UTR" = "#fc8d62",
  "Exon" = "#8da0cb",
  "Intergenic" = "#e78ac3",
  "Intron" = "#a6d854"
)

############################################################
# Annotate mutation sites
############################################################

mutation_site_annotations =
  annotate_sites(
    dirclip_files,
    exons_GR,
    introns_GR,
    utrs_GR,
    gtf_GR,
    annot_type,
    gene_type
  )

############################################################
# Generate stacked bar plot
############################################################

y_tick_interval = 0.2
y_lim = c()

Mut_sites_annotation_bar_plot(
  mutation_site_annotations,
  Colors,
  y_lim,
  y_tick_interval,
  fig_name = figure_name,
  fig_sz = c(6, 7)
)