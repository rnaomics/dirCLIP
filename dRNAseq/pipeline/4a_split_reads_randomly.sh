#!/usr/bin/env bash
# 4a_split_reads_randomly.sh: shuffle read IDs and split into 3 pseudo-replicates.
# Can be applied to the nanopolish output or the extracted read-coordinates file.
# BUGFIX: original shuffled "${file}_reads.txt" which was never created (the dedup step
#         wrote "${file}_unique_reads.txt"). The chain now uses one consistent name.
set -euo pipefail

# ---------------- CONFIG ----------------
file="${file:-sample_name_nano}"        # basename WITHOUT extension; expects ${file}.tsv to exist
split_dir="${split_dir:-${file}_split}"  # directory for split intermediates/outputs
# ----------------------------------------
mkdir -p "${split_dir}"

# read-id column (4) from the nanopolish TSV
cut -f 4 -d $'\t' "${file}.tsv"                       > "${split_dir}/${file}_read_ids.tsv"
# drop header
tail -n +2 "${split_dir}/${file}_read_ids.tsv"       > "${split_dir}/${file}_no_header.tsv"
# sort + dedup read ids
sort "${split_dir}/${file}_no_header.tsv"            > "${split_dir}/${file}_no_header_sorted.tsv"
uniq "${split_dir}/${file}_no_header_sorted.tsv"    > "${split_dir}/${file}_unique_reads.txt"
# shuffle (now reads the dedup output) then split into 3 equal groups -> r00 r01 r02
shuf "${split_dir}/${file}_unique_reads.txt" -o "${split_dir}/${file}_reads_shuff.txt"
split -n l/3 -d "${split_dir}/${file}_reads_shuff.txt" "${split_dir}/r"
