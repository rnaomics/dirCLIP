#### 📦 Set library path ####
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

#### 🧬 Bambu Installation (optional) ####
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("bambu")

# 📖 Reference guide:
# https://www.bioconductor.org/packages/devel/data/experiment/vignettes/NanoporeRNASeq/inst/doc/NanoporeRNASeq.html

#### 📚 Load required libraries ####
library(bambu)
library(dplyr)

#### 🧩 Function to convert index to transcript ID ####
index_to_transcript_id = function(Annot, x) {
  if (is.null(x)) {
    return(NULL)  # Handle NULL cases
  } else {
    return(Annot[x])
  }
}

#### 📁 Load parameters from key-value file ####
params <- read.table("parameters.txt", sep = "=", stringsAsFactors = FALSE, header = FALSE)
params_list <- setNames(as.list(params$V2), params$V1)

# 🔗 Extract individual parameters
Project_path <- params_list$Project_dir
Bam_files <- strsplit(params_list$Bam_files, ",")[[1]]
bambu_dirs <- strsplit(params_list$bambu_dirs, ",")[[1]]
fa.file <- paste0(Project_path, params_list$genome_file)
gtf.file <- paste0(Project_path, params_list$gtf_file)

#### 🏗️ Prepare annotation data ####
bambuAnnotations <- prepareAnnotations(gtf.file)

#### 🔄 Process BAM files iteratively ####
for (i in 1:length(Bam_files)) {
  # ⛓️ Define input/output paths
  bam.file <- paste0(Project_path, Bam_files[i])
  Result_path <- paste0(Project_path, bambu_dirs[i], "Out/")
  
  # Run bambu with discovery disabled
  se <- bambu(reads = bam.file, annotations = bambuAnnotations, genome = fa.file, trackReads = TRUE, discovery = FALSE)
  
  # Extract transcript assignments
  Transcript_assign <- metadata(se)$readToTranscriptMaps[[1]]
  
  #### 🔠 Translate transcript indices to names ####
  Annot <- rowData(se)@listData[["TXNAME"]]
  Annot <- gsub("-[0-9]+", "", Annot)  # Remove suffixes like "-1", "-2", etc.
  names(Annot) <- rowData(se)@listData[["txid"]]
  
  Transcript_assign$equalMatches <- lapply(Transcript_assign$equalMatches, function(x) index_to_transcript_id(Annot, x))
  Transcript_assign$compatibleMatches <- lapply(Transcript_assign$compatibleMatches, function(x) index_to_transcript_id(Annot, x))
  
  # 🧽 Flatten list columns to strings
  Transcript_assign <- Transcript_assign %>% 
    mutate(across(where(is.list), ~ sapply(., toString)))
  
  #### 💾 Write output ####
  if (!dir.exists(Result_path)) dir.create(Result_path, recursive = TRUE)
  writeBambuOutput(se, path = Result_path)
  
  write.table(
    Transcript_assign,
    paste0(Result_path, "readToTranscriptMaps.tsv"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
}
