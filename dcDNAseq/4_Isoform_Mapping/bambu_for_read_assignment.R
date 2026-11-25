############################################################
## Isoform Assignment Pipeline using bambu
#
## Purpose:
##   - Perform isoform-level read assignment using bambu
##   - Map reads to transcripts using provided annotations
##   - Export transcript assignment per read
##
## Input:
##   - BAM files (aligned long-read RNA-seq)
##   - Genome FASTA file
##   - GTF annotation file
##
## Output:
##   - bambu expression output per sample
##   - read-to-transcript assignment table (TSV)

###############################
## 1. Environment Setup
###############################

# Set custom R library path (adjust if needed for your system)
# .libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.6")

# Load required libraries
library(bambu)
library(dplyr)

############################################################
## 2. Helper Functions
############################################################

############################################################
## Convert transcript indices returned by bambu to names
##
## Arguments:
##   Annot : named vector mapping transcript indices → names
##   x     : vector of transcript indices
##
## Returns:
##   Character vector of transcript IDs (or NULL if input is NULL)
############################################################
index_to_transcript_id <- function(Annot, x) {
  if (is.null(x)) {
    return(NULL)
  } else {
    return(Annot[x])
  }
}

############################################################
## 3. Load Parameters
############################################################

# Load key-value parameter file (format: key=value)
params <- read.table(
  "parameters.txt",
  sep = "=",
  stringsAsFactors = FALSE,
  header = FALSE
)

# Convert parameters into a named list
params_list <- setNames(as.list(params$V2), params$V1)

# Extract parameters
Bam_files  <- strsplit(params_list$Bam_files, ",")[[1]]
bambu_dirs <- strsplit(params_list$bambu_dirs, ",")[[1]]

fa.file  <- params_list$genome_file   # Genome FASTA
gtf.file <- params_list$gtf_file      # Annotation GTF

############################################################
## 4. Prepare Annotation for bambu
############################################################

# Prepare annotation object used by bambu
bambuAnnotations <- prepareAnnotations(gtf.file)

############################################################
## 5. Main Processing Loop
############################################################

for (i in seq_along(Bam_files)) {
  
  ###############################
  ## 5.1 Define Input/Output
  ###############################
  
  bam.file   <- Bam_files[i]
  Result_dir <- paste0(bambu_dirs[i], "Out/")
  
  # Create output directory if it does not exist
  ifelse(!dir.exists(Result_dir), dir.create(Result_dir, recursive = TRUE), FALSE)
  
  message("Processing: ", bam.file)
  
  ###############################
  ## 5.2 Run bambu
  ###############################
  
  # Run bambu without novel transcript discovery
  # (set discovery = TRUE if novel isoform detection is required)
  se <- bambu(
    reads = bam.file,
    annotations = bambuAnnotations,
    genome = fa.file,
    trackReads = TRUE,
    discovery = FALSE
  )
  
  
  ###############################
  ## 5.3 Extract Read-to-Transcript Mapping
  ###############################
  
  # Retrieve mapping object
  Transcript_assign <- metadata(se)$readToTranscriptMaps[[1]]
  
  
  ###############################
  ## 5.4 Map Transcript IDs
  ###############################
  
  # Extract transcript name mapping
  Annot <- rowData(se)@listData[["TXNAME"]]
  
  # Remove transcript version suffix (e.g. "-1")
  Annot <- gsub("-[0-9]+", "", Annot)
  
  # Assign names to map indices → transcript IDs
  names(Annot) <- rowData(se)@listData[["txid"]]
  
  
  ###############################
  ## 5.5 Replace Index-Based Matches with IDs
  ###############################
  
  # Convert equalMatches and compatibleMatches indices to transcript IDs
  Transcript_assign$equalMatches <- lapply(
    Transcript_assign$equalMatches,
    function(x) index_to_transcript_id(Annot, x)
  )
  
  Transcript_assign$compatibleMatches <- lapply(
    Transcript_assign$compatibleMatches,
    function(x) index_to_transcript_id(Annot, x)
  )
  
  
  ###############################
  ## 5.6 Flatten Output Table
  ###############################
  
  # Convert list columns into comma-separated strings
  Transcript_assign <- Transcript_assign %>%
    mutate(across(where(is.list), ~ sapply(., toString)))
  
  
  ###############################
  ## 5.7 Write Output Files
  ###############################
  
  # 1. Standard bambu output (counts, transcripts, etc.)
  writeBambuOutput(se, path = Result_dir)
  
  # 2. Read-to-transcript assignment table
  write.table(
    Transcript_assign,
    file = paste0(Result_dir, "readToTranscriptMaps.tsv"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  
  message("Completed: ", bam.file)
}


############################################################
## End of Script
############################################################
