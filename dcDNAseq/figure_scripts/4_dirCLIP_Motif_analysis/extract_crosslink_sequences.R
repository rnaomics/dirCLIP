# Extract sequences around dirCLIP crosslinking sites

# setwd("/data/groups/anko/mohamed/oxford_nanopore/Submission/Tools/")
setwd("Tools/")
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

library(writexl)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(scales)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm9)

genome <- BSgenome.Mmusculus.UCSC.mm9

######### Function to load the sites with some offset ####################

load_sites = function(site_file, Left_offset, Right_offset){
  
  ############ Load the CLIP tags #######################
  
  site_df = read.delim(
    site_file,
    header = FALSE
  )
  
  site_df$strand = "*"
  
  ########### Remove random chromosomes #################
  
  random_chr_idx =
    grep("random", site_df$V1)
  
  if (!isEmpty(random_chr_idx)){
    site_df = site_df[-random_chr_idx,]
  }
  
  ######### Creating GRanges object #####################
  
  site_granges <- GRanges(
    seqnames = site_df[, "V1"],
    IRanges(
      site_df[, "V2"] - Left_offset,
      site_df[, "V2"] + Right_offset
    ),
    site_df$strand
  )
  
  ######### Removing negative ranges ####################
  
  valid_range_idx =
    start(site_granges) >= 0
  
  site_granges =
    site_granges[valid_range_idx]
  
  return(site_granges)
}

############## Function to load the regions with some offset ###############

load_regions = function(region_file,
                        Left_offset,
                        Right_offset){
  
  region_df = read.delim(
    region_file,
    header = TRUE
  )
  
  region_df$strand = "*"
  
  region_df =
    region_df[!is.na(region_df$start), ]
  
  random_chr_idx =
    grep("random", region_df$chr)
  
  if (!isEmpty(random_chr_idx)){
    region_df =
      region_df[-random_chr_idx, ]
  }
  
  region_granges <- GRanges(
    seqnames = region_df[, "chr"],
    IRanges(
      region_df[, "start"] - Left_offset,
      region_df[, "end"] + Right_offset
    ),
    region_df$strand
  )
  
  valid_range_idx =
    start(region_granges) > 0
  
  region_granges =
    region_granges[valid_range_idx]
  
  return(region_granges)
}

########## Function to load dirCLIP sites, remove random chromosomes,
########## add strand information and generate site windows ##########

get_site_windows = function(
    dirCLIP_file,
    genes_gr,
    win
){
  
  ###### 1. Load dirCLIP sites ##########################
  
  dirclip_df <- read.table(
    dirCLIP_file,
    sep = "\t"
  )
  
  colnames(dirclip_df) <- c(
    "chr",
    "pos",
    "score"
  )
  
  ###### 2. Remove random chromosomes ##################
  
  filtered_dirclip_df <- dirclip_df[
    !grepl("random", dirclip_df$chr),
  ]
  
  ###### 3. Create GRanges object ######################
  
  site_windows_gr <- GRanges(
    seqnames = filtered_dirclip_df$chr,
    ranges = IRanges(
      start = filtered_dirclip_df$pos,
      end = filtered_dirclip_df$pos
    )
  )
  
  site_windows_gr <- resize(
    site_windows_gr,
    width = win,
    fix = "center"
  )
  
  ######### Remove negative ranges #####################
  
  valid_range_idx =
    start(site_windows_gr) >= 0
  
  site_windows_gr =
    site_windows_gr[valid_range_idx]
  
  ########## 4. Add strand information #################
  
  overlap_hits =
    findOverlaps(
      site_windows_gr,
      genes_gr
    )
  
  site_windows_df =
    data.frame(
      seqnames =
        seqnames(site_windows_gr)[queryHits(overlap_hits)],
      
      start =
        start(site_windows_gr)[queryHits(overlap_hits)],
      
      end =
        end(site_windows_gr)[queryHits(overlap_hits)],
      
      strand =
        strand(genes_gr)[subjectHits(overlap_hits)]
    )
  
  site_windows_df =
    unique(site_windows_df)
  
  ######################################################
  # Assign strand "*" if both strands overlap
  ######################################################
  
  site_windows_df$strand =
    as.character(site_windows_df$strand)
  
  site_windows_df =
    site_windows_df %>%
    group_by(seqnames, start, end) %>%
    summarise(
      strand =
        if (n_distinct(strand) > 1)
          "*"
      else
        unique(strand),
      .groups = "drop"
    )
  
  ######################################################
  # Remove ambiguous strand assignments
  ######################################################
  
  unknown_strand_idx =
    site_windows_df$strand == "*"
  
  site_windows_df =
    site_windows_df[!unknown_strand_idx,]
  
  ######################################################
  # Convert back to GRanges
  ######################################################
  
  site_windows_gr <- GRanges(
    seqnames = site_windows_df$seqnames,
    ranges = IRanges(
      start = site_windows_df$start,
      end = site_windows_df$end
    ),
    strand = site_windows_df$strand
  )
  
  return(site_windows_gr)
}

######################## Create folder ######################################################
create_folder = function(Path){
  ifelse(!dir.exists(Path), dir.create(Path, recursive = TRUE), FALSE)
}

#######################################################################################

paths <- c('figure_scripts/figure_data/1_Mutation_profiling/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/')
subfolders2 <- c('2_Inc_del', '1_No_del')
Res_paths <- c('Result/figure/2_Motif_analysis/1_Genome_Mapping/1_SRSF3_cDNA_Mg/ref_TotalRNA_cDNA_no_UV/')
Coverage_thresholds <- c(50)
Diff_thresholds <- c(0.03)
wins <- c(100)

gtf_file <- 'genome/mm9/Annotation/gencode.vM1.annotation.gtf'

genome_GR <- import(gtf_file)
genome_GR <- genome_GR[!grepl('_random', seqnames(genome_GR))]
genome_GR <- genome_GR[genome_GR$type == 'gene']

for (i in seq_along(paths)) {
  for (subfolder2 in subfolders2) {
    for (cov_thresh in Coverage_thresholds) {

      Mut_path <- paste0(paths[i], subfolder2, '/cov_thresh_', cov_thresh, '/')

      for (mut_thresh in Diff_thresholds) {
        for (win in wins) {
          Res_path <- paste0(Res_paths[i], subfolder2, '/cov_thresh_', cov_thresh, '/sites/win_', win, "/")
  
          dirCLIP_file <- paste0(Mut_path, 'Pos_mut_thresh_', mut_thresh, '.tsv')

          gr <- get_site_windows(dirCLIP_file, genome_GR, win)

          create_folder(Res_path)

          bedfile <- paste0(Res_path,  'Pos_mut_thresh_', mut_thresh, '.bed.gz')
          export(gr, con = bedfile, format = 'bed')

          gr_seq <- getSeq(genome, gr)
          names(gr_seq) <- paste0('ID_', seq_along(gr_seq))

          gr_fa <- paste0(Res_path, 'Pos_mut_thresh_', mut_thresh, '.fasta.gz')

          writeXStringSet(
            gr_seq,
            filepath = gr_fa,
            format = 'fasta',
            compress = TRUE
          )
        }
      }
    }
  }
}
