#!/usr/bin/env bash
# 3_read_coordinate_extraction.sh: extract per-read start/end coordinates from the BAM
# BUGFIX: original used Snakemake-style {input}/{output} braces inside a plain bash command,
#         which bash does not expand. Replaced with proper "$input"/"$output" variables.
set -euo pipefail

# ---------------- CONFIG ----------------
input="${input:-sample_bam_file.sorted.bam}"           # sorted BAM from step 2
output="${output:-sample_read_coordinates_file.txt}"   # read-coordinate table
# ----------------------------------------

# columns kept: 1=transcript, 2=read start, 3=read end, 4=read name, 6=strand
bedtools bamtobed -i "${input}" | cut -f 1,2,3,4,6 > "${output}"
