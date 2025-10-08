################################################################################
# gene_count_pipeline.R
# Author: Mohamed B
# Description:
# Loads project-specific parameters and quantifies gene-level read counts
# using Subread::featureCounts. Supports exon and transcript features.
################################################################################

### 📁 Set Library Path (R 4.3)
.libPaths("/data/groups/anko/mohamed/Install_packages/R_packages/x86_64-pc-linux-gnu-library/4.3")

### 📚 Load Required Libraries
library(Rsubread)
library(writexl)

### 📄 Load External Parameters
Load_from_parameters <- TRUE
if (Load_from_parameters) {
  params <- read.table("parameters.txt", sep = "=", stringsAsFactors = FALSE, header = FALSE)
  params_list <- setNames(as.list(params$V2), params$V1)
}

### 📂 Set Working Directory
Project_path <- if (Load_from_parameters) params_list$Project_dir else "/data/groups/anko/mohamed/oxford_nanopore/"
setwd(Project_path)

### 🧭 Input Paths and Metadata
if (Load_from_parameters) {
  Bam_path        <- params_list$Path
  Res_path        <- params_list$Path
  bam_file_name   <- params_list$bam_file_name
  gene_count_file <- params_list$gene_count_file
  Annot_file      <- params_list$Annot_file
  Gtf_attr_type   <- params_list$Gtf_attr_type
} else {
  Bam_path        <- "Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/CLIP_cDNA_rep1/"
  Res_path        <- "Results/Nanopore/2_Hmga1_DIRCLIP/Alignment/minimap/mm9/Mauro/CLIP_cDNA_rep1/"
  bam_file_name   <- "sorted.bam"
  gene_count_file <- "Genewise_readcounts.tsv"
  Annot_file      <- "genome/mm9/Annotation/gencode.vM1.annotation.gtf"
  Gtf_attr_type   <- "gene_name"  # or "gene_id" depending on GTF
}

### 🧬 Features to Quantify
features <- c("exon", "transcript")

### 🧮 Initialize Result Matrix
Count_matrix <- NULL

### 🔁 Loop Over Features and Run featureCounts
for (feature in features) {
  bam_file <- file.path(Bam_path, bam_file_name)
  
  Gene_wise_counts <- featureCounts(
    files               = bam_file,
    annot.ext           = Annot_file,
    isPairedEnd         = FALSE,
    isGTFAnnotationFile = TRUE,
    GTF.featureType     = feature,
    GTF.attrType        = Gtf_attr_type,
    allowMultiOverlap   = TRUE,
    fraction            = TRUE,
    isLongRead          = TRUE,
    primaryOnly         = TRUE
  )
  
  # Round counts and rename column
  Counts <- round(Gene_wise_counts[["counts"]])
  colnames(Counts) <- paste0(feature, "_reads")
  
  # Combine into matrix
  Count_matrix <- cbind(Count_matrix, Counts)
  
  # Save annotation for final output
  Annotation <- Gene_wise_counts$annotation
}

### 📊 Assemble Final Data Frame
df <- data.frame(
  Gene_id = rownames(Count_matrix),
  Count_matrix,
  Start   = Annotation$Start,
  End     = Annotation$End,
  Strand  = Annotation$Strand,
  Length  = Annotation$Length,
  Chr     = Annotation$Chr
)

### 💾 Save Output Table
write.table(
  df,
  file      = file.path(Res_path, gene_count_file),
  sep       = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE
)

###################################################################


