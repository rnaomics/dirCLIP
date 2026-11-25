#!/usr/bin/env Rscript
# dcDNA and dRNA dirCLIP Overlap Analysis
# Manuscript-specific figure-generation and overlap-analysis script.
#
# Purpose:
# Compare dcDNA and dRNA dirCLIP sites/regions, calculate overlap
# statistics, and generate UpSet visualizations.

### set work directory ##############
# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")

#### Set library path #########
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)   # makeTxDbFromGFF(), exonsBy()
  library(rtracklayer)       # import.chain(), liftOver()
  library(readxl)            # read_excel()
  library(ggplot2)
  library(ggupset)           # scale_x_upset()
  library(dplyr)             # loaded last so that dplyr verbs win over masking
})

################# Function to convert transcript coordinates to genome coordinates ###################

convert_tx_coordinates_to_genome_coordinates = function(txdb, tx_sites){
  
  tx_sites$ref_id <- sub("\\..*", "", tx_sites$ref_id)
  
  ############ Extract transcript structures #############################
  
  exons_by_transcript <- exonsBy(txdb, by = "tx", use.names = TRUE)
  
  names(exons_by_transcript) <- sub(
    "\\..*",
    "",
    names(exons_by_transcript)
  )
  
  ########## Convert transcript coordinates to GRanges ###################
  
  transcript_sites_gr <- GRanges(
    seqnames = tx_sites$ref_id,
    ranges = IRanges(
      start = tx_sites$region_start,
      end = tx_sites$range_end
    )
  )
  
  ###### Map transcript -> genome ########################################
  
  genomic_sites_gr <- mapFromTranscripts(
    transcript_sites_gr,
    exons_by_transcript
  )
  
  unmapped_site_idx <- !(1:length(transcript_sites_gr) %in% genomic_sites_gr$xHits)
  
  unmapped_sites_gr <- transcript_sites_gr[unmapped_site_idx]
  
  return(genomic_sites_gr)
}

################ Function to load dirCLIP sites #####################################################

load_dirCLIP_sites = function(dirCLIP_path, dcdna_type, txdb){
  
  if (dcdna_type == "sites"){
    
    dirclip_df <- read.delim(dirCLIP_path, header = FALSE)
    
    dirclip_df = dirclip_df[, c("V1", "V2")]
    
    dirclip_df$V1 <- sub("^Chr_", "", dirclip_df$V1)
    
    colnames(dirclip_df) = c("chr", "start")
    
    dirclip_df$end = dirclip_df$start
    
  } else if (dcdna_type == "regions"){
    
    dirclip_df <- read.delim(dirCLIP_path, header = TRUE)
    
    dirclip_df$chr <- sub("^Chr_", "", dirclip_df$chr)
  }
  
  dirclip_df$chr <- sub("\\..*", "", dirclip_df$chr)
  
  ########## Convert transcript coordinates to GRanges ###################
  
  transcript_sites_gr <- GRanges(
    seqnames = dirclip_df$chr,
    ranges = IRanges(
      start = dirclip_df$start,
      end = dirclip_df$end
    )
  )
  
  ############ Extract transcript structures #############################
  
  exons_by_transcript <- exonsBy(txdb, by = "tx", use.names = TRUE)
  
  names(exons_by_transcript) <- sub(
    "\\..*",
    "",
    names(exons_by_transcript)
  )
  
  ###### Map transcript → genome #########################################
  
  genomic_sites_gr <- mapFromTranscripts(
    transcript_sites_gr,
    exons_by_transcript
  )
  
  genomic_sites_gr = unique(genomic_sites_gr)
  
  unmapped_site_idx <- !(1:length(transcript_sites_gr) %in% genomic_sites_gr$xHits)
  
  unmapped_sites_gr <- transcript_sites_gr[unmapped_site_idx]
  
  return(genomic_sites_gr)
}

################# Function to get overlap between dcDNA_dirCLIP and dRNA_dirCLIP sites ###############

get_dcDNA_dRNA_overlap_sites = function(dcDNA_dirCLIP, dRNA_dirCLIP){
  
  overlap_hits <- findOverlaps(dcDNA_dirCLIP, dRNA_dirCLIP)
  
  dcDNA_overlap_idx <- unique(queryHits(overlap_hits))
  dRNA_overlap_idx <- unique(subjectHits(overlap_hits))
  
  ###### dcDNA → dRNA overlap ############################################
  
  dcDNA_total_sites <- length(dcDNA_dirCLIP)
  dcDNA_overlapping_sites <- length(dcDNA_overlap_idx)
  
  dcDNA_overlap_pct <- dcDNA_overlapping_sites / dcDNA_total_sites * 100
  
  ###### dRNA → dcDNA overlap ############################################
  
  dRNA_total_sites <- length(dRNA_dirCLIP)
  dRNA_overlapping_sites <- length(dRNA_overlap_idx)
  
  dRNA_overlap_pct <- dRNA_overlapping_sites / dRNA_total_sites * 100
  
  overlap_summary_df <- data.frame(
    dcDNA_total = dcDNA_total_sites,
    dcDNA_overlap_dRNA = dcDNA_overlapping_sites,
    dRNA_total = dRNA_total_sites,
    dRNA_overlap_dcDNA = dRNA_overlapping_sites,
    dcDNA_pct = dcDNA_overlap_pct,
    dRNA_pct = dRNA_overlap_pct
  )
  
  return(overlap_summary_df)
}

########## Function to make ggupset plot #####################################

make_ggupset = function(df, Title, text_size, text_size2, y_lim,
                        y_tick_interval, fig_name, fig_sz){
  
  # Generate list of conditions per row
  group_membership <- apply(df[, -1], 1, function(x) names(x)[x == 1])
  
  upset_plot_df <- tibble(
    ID = df$ID,
    group = group_membership
  )
  
  y_lab = "Intersection size"
  x_lab = "group"
  
  upset_plot = ggplot(upset_plot_df, aes(x = group)) +
    geom_bar() +
    scale_x_upset() +
    geom_text(
      stat = "count",
      aes(label = after_stat(count)),
      vjust = -0.5,
      size = text_size2
    ) +
    labs(
      y = y_lab,
      x = x_lab,
      title = Title
    ) +
    theme_minimal(text_size) +
    theme(
      text = element_text(size = text_size),
      
      panel.grid = element_blank(),
      
      panel.border = element_blank(),
      
      axis.line = element_line(color = "black", linewidth = 0.8),
      
      axis.text.x = element_text(
        size = text_size,
        color = "black",
        angle = 0,
        hjust = 1
      ),
      
      axis.title.x = element_text(size = text_size),
      
      axis.text.y = element_text(
        size = text_size,
        color = "black"
      ),
      
      axis.title.y = element_text(size = text_size),
      
      plot.title = element_text(
        size = text_size,
        hjust = 0.5
      ),
      
      axis.ticks = element_line(
        color = "black",
        size = 0.8
      ),
      
      axis.ticks.length = unit(0.25, "cm"),
      
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black")
    )
  
  if (length(y_lim) == 0){
    
    y_ticks <- pretty(
      c(0, max(table(df[, -1])) * 1.1),
      n = 5
    )
    
    y_min = min(y_ticks)
    y_max = max(y_ticks)
    
  } else {
    
    y_min = floor(y_lim[1] / y_tick_interval) * y_tick_interval
    y_max = ceiling(y_lim[2] / y_tick_interval) * y_tick_interval
    
    y_ticks = seq(
      y_min,
      y_max,
      by = y_tick_interval
    )
  }
  
  upset_plot = upset_plot +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = y_ticks,
      expand = expansion(mult = c(0, 0.02))
    )
  
  ggsave(
    fig_name,
    plot = upset_plot,
    width = fig_sz[1],
    height = fig_sz[2],
    limitsize = FALSE,
    device = cairo_pdf
  )
  
  return(upset_plot)
}

##### Function to get overlaps for upset plot ###########

get_overlaps_for_upset_plot = function(GR_list, types){
  
  combined_gr <- unique(do.call(c, GR_list))
  
  overlap_df <- data.frame(ID = 1:length(combined_gr))
  
  for (group_idx in 1:length(GR_list)){
    
    overlap_hits <- findOverlaps(combined_gr, GR_list[[group_idx]])
    
    overlap_presence <- 1:length(combined_gr) %in% queryHits(overlap_hits)
    
    overlap_df[, types[group_idx]] = overlap_presence
  }
  
  return(overlap_df)
}

############ Load dRNA sites and convert to genomic sites ######################

load_dRNA_sites = function(
    dRNA_dirCLIP_path,
    Nanocompore_sheet,
    Read_termination_sheet,
    txdb
){
  
  nanocompore_df <- read_excel(
    dRNA_dirCLIP_path,
    sheet = Nanocompore_sheet
  )
  
  read_termination_df <- read_excel(
    dRNA_dirCLIP_path,
    sheet = Read_termination_sheet
  )
  
  # Convert columns
  
  nanocompore_df <- nanocompore_df %>%
    mutate(
      LOR = as.numeric(LOR),
      GMM_adj_pval = as.numeric(GMM_adj_pval)
    )
  
  read_termination_df <- read_termination_df %>%
    mutate(
      L2FC = as.numeric(L2FC),
      adj_pval = as.numeric(adj_pval)
    )
  
  nanocompore_sig_sites <- nanocompore_df %>%
    filter(
      abs(LOR) >= sigdev_thresh$lor_thresh,
      GMM_adj_pval < sigdev_thresh$pval_thresh
    ) %>%
    select(ref_id, pos)
  
  read_termination_sig_sites <- read_termination_df %>%
    filter(
      L2FC >= read_end_thresh$lfc_thresh,
      adj_pval < read_end_thresh$pval_thresh
    ) %>%
    select(ref_id, pos)
  
  transcript_sites_df <- bind_rows(
    nanocompore_sig_sites,
    read_termination_sig_sites
  ) %>%
    distinct()
  
  transcript_sites_df$ref_id <- sub(
    "\\..*",
    "",
    transcript_sites_df$ref_id
  )
  
  ############ Extract transcript structures #############################
  
  exons_by_transcript <- exonsBy(txdb, by = "tx", use.names = TRUE)
  
  names(exons_by_transcript) <- sub(
    "\\..*",
    "",
    names(exons_by_transcript)
  )
  
  ########## Convert transcript coordinates to GRanges ###################
  
  transcript_sites_gr <- GRanges(
    seqnames = transcript_sites_df$ref_id,
    ranges = IRanges(
      start = transcript_sites_df$pos,
      end = transcript_sites_df$pos
    )
  )
  
  transcript_sites_gr = unique(transcript_sites_gr)
  
  ###### Map transcript → genome #########################################
  
  genomic_sites_gr <- mapFromTranscripts(
    transcript_sites_gr,
    exons_by_transcript
  )
  
  unmapped_site_idx <- !(1:length(transcript_sites_gr) %in% genomic_sites_gr$xHits)
  
  unmapped_sites_gr <- transcript_sites_gr[unmapped_site_idx]
  
  return(genomic_sites_gr)
}
################################################################################
# 1. SRSF3 dcDNA-dRNA Sites (-UV)
################################################################################

dcDNA_dirCLIP_path =
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03.tsv"

olap_tsv_path =
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03_olap_with_dRNA_sites.tsv"

upset_path =
  "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_TotalRNA_cDNA_no_UV/FDR/dcDNA_dRNA_sites_overlap.pdf"

dcDNA_type = "sites"

dRNA_dirCLIP_path = "Supplementary Table 3_dRNAseq_SRSF3.xlsx"

Nanocompore_sheet     = "-UV_SRSF3_Nanocompore"
Read_termination_sheet = "-UV_SRSF3_Read_Termination"

read_end_thresh = data.frame(
  lfc_thresh  = 2,
  pval_thresh = 0.05
)

sigdev_thresh = data.frame(
  lor_thresh  = 1,
  pval_thresh = 0.005
)

################################################################################
# SRSF3 Regions
################################################################################

# dcDNA_dirCLIP_path =
#   "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03_region.tsv"
#
# olap_tsv_path =
#   "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03_olap_with_dRNA_regions.tsv"
#
# upset_path =
#   "Results/Nanopore/3_SRSF3_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_TotalRNA_cDNA_no_UV/FDR/dcDNA_dRNA_regions_overlap.pdf"
#
# dcDNA_type = "regions"
#
# dRNA_dirCLIP_path  = "Supplementary Table 3_dRNAseq_SRSF3.xlsx"
# dRNA_dirCLIP_sheet = "-UV_SRSF3_Sig_Regions"
#
# cov_thresh  = 50
# diff_thresh = 0.03


################################################################################
# 2. HMGA1A dcDNA-dRNA Sites (-UV)
################################################################################

# dcDNA_dirCLIP_path =
#   "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/Hmga1A_CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03.tsv"
# 
# olap_tsv_path =
#   "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/Hmga1A_CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03_olap_with_dRNA_sites.tsv"
# 
# upset_path =
#   "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/Hmga1A_CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/FDR/dcDNA_dRNA_sites_overlap.pdf"
# 
# dcDNA_type = "sites"
# 
# dRNA_dirCLIP_path = "Supplementary Table 9_dRNAseq_HMGA1.xlsx"
# 
# Nanocompore_sheet      = "noUV_polyA_HMGA1_Nanocompore"
# Read_termination_sheet = "noUVpolyA_HMGA1_Read_Terminatio"
# 
# read_end_thresh = data.frame(
#   lfc_thresh  = 2,
#   pval_thresh = 0.05
# )
# 
# sigdev_thresh = data.frame(
#   lor_thresh  = 1,
#   pval_thresh = 0.05
# )

################################################################################
# HMGA1A Regions
################################################################################

# dcDNA_dirCLIP_path =
#   "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/Hmga1A_CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03_region.tsv"
#
# olap_tsv_path =
#   "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/Hmga1A_CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/cov_thresh_50/FDR/Pos_mut_thresh_0.03_olap_with_dRNA_regions.tsv"
#
# upset_path =
#   "Results/Nanopore/2_Hmga1_DIRCLIP/Mutation_profiling/bambu/mm9_transcriptome/Mauro/Hmga1A_CLIP_cDNA_merge/q_8_Q_8/gencode/subsample_70/2_Inc_del/ref_P19_WT_TotalRNA_cDNA_no_UV/FDR/dcDNA_dRNA_regions_overlap.pdf"
#
# dcDNA_type = "regions"
#
# dRNA_dirCLIP_path  = "Supplementary Table 9_dRNAseq_HMGA1.xlsx"
# dRNA_dirCLIP_sheet = "noUV_polyA_HMGA1_Sig_regions"

cov_thresh  = 50
diff_thresh = 0.03

################################################################################
# Annotation Files
################################################################################

mm9_annotation  =
  "genome/mm9/Annotation/gencode.vM1.annotation.gtf.gz"

mm39_annotation =
  "genome/mm39/Annotation/gencode.vM36.annotation.gtf.gz"

# mm39_to_mm9_chain =
#   "genome/mm39/Annotation/mm39ToMm9.over.chain"

mm9_to_mm39_chain =
  "genome/mm39/Annotation/mm9ToMm39.over.chain"

################################################################################
# Load Annotation and Chain
################################################################################

txdb_mm9  <- makeTxDbFromGFF(mm9_annotation)
txdb_mm39 <- makeTxDbFromGFF(mm39_annotation)

mm9_to_mm39_chain_gr <- import.chain(mm9_to_mm39_chain)

################################################################################
# Load dRNA Sites / Regions
################################################################################
if (dcDNA_type == "sites") {
  
  dRNA_dirCLIP_mm39 = load_dRNA_sites(
    dRNA_dirCLIP_path,
    Nanocompore_sheet,
    Read_termination_sheet,
    txdb_mm39
  )
  
  dRNA_dirCLIP_mm39 <- resize(
    dRNA_dirCLIP_mm39,
    width = 11,
    fix = "center"
  )
  
} else if (dcDNA_type == "regions") {
  
  dRNA_dirCLIP_df <- read_excel(
    dRNA_dirCLIP_path,
    sheet = dRNA_dirCLIP_sheet
  )
  
  dRNA_dirCLIP_mm39 =
    convert_tx_coordinates_to_genome_coordinates(
      txdb_mm39,
      dRNA_dirCLIP_df
    )
}


################################################################################
# Load dcDNA Sites / Regions
################################################################################
dcDNA_dirCLIP_mm9 =
  load_dirCLIP_sites(
    dcDNA_dirCLIP_path,
    dcDNA_type,
    txdb_mm9
  )

dcDNA_dirCLIP_mm9 <- resize(
  dcDNA_dirCLIP_mm9,
  width = 11,
  fix = "center"
)

dcDNA_dirCLIP_mm39 <- liftOver(
  dcDNA_dirCLIP_mm9,
  mm9_to_mm39_chain_gr
)

dcDNA_dirCLIP_mm39 <- unlist(dcDNA_dirCLIP_mm39)

################################################################################
# Overlap Analysis
################################################################################
overlap_summary_df =
  get_dcDNA_dRNA_overlap_sites(
    dcDNA_dirCLIP_mm39,
    dRNA_dirCLIP = dRNA_dirCLIP_mm39
  )

write.table(
  overlap_summary_df,
  file = olap_tsv_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

################################################################################
# UpSet Plot
################################################################################
upset_input_df =
  get_overlaps_for_upset_plot(
    GR_list = list(
      dcDNA_dirCLIP_mm39,
      dRNA_dirCLIP_mm39
    ),
    types = c("dcDNA", "dRNA")
  )

upset_plot =
  make_ggupset(
    upset_input_df,
    Title = paste0(
      "Cov_thresh = ",
      cov_thresh,
      ", Diff_thresh = ",
      diff_thresh
    ),
    text_size = 16,
    text_size2 = 6,
    y_lim = c(),
    y_tick_interval = 10000,
    fig_name = upset_path,
    fig_sz = c(5, 7)
  )

################################################################################################
