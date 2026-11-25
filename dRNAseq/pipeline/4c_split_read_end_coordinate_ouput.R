# 4c_split_read_end_coordinate_ouput.R: subset the bedtools read-coordinate file by replicate.
# Brushed: paths supplied as arguments (or edit the defaults below). Logic unchanged.
suppressMessages(library(tidyverse))

# ---- arguments: Rscript 4c_*.R <read_ids_file> <bedtools_output> <out_file> ----
args <- commandArgs(trailingOnly = TRUE)
read_ids_file <- if (length(args) >= 1) args[1] else "path/to/split_replicate_files/r00"
bedtools_file <- if (length(args) >= 2) args[2] else "bam_to_bed_output/sample_read_coordinates_file.txt"
out_file      <- if (length(args) >= 3) args[3] else "sample_read_coordinates_r00.txt"

read_ids <- read.table(read_ids_file, header = FALSE, sep = "\t")
bedtools_output_full <- read.table(bedtools_file, header = FALSE, sep = "\t")
colnames(bedtools_output_full) <- c("Transcript", "Start", "End", "Read_ID", "Strand")

bedtools_rep <- bedtools_output_full %>% filter(Read_ID %in% read_ids$V1)
write.table(bedtools_rep, out_file, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE)  # no header row: script 8 reads this file with its own column names (corrected 10 Oct 2026)
cat("Written", nrow(bedtools_rep), "rows to", out_file, "\n")
