#!/usr/bin/env bash
# 4b_split_nanopolish_ouput.sh: subset the nanopolish output into 3 pseudo-replicates
# using the shuffled read-id groups (r00/r01/r02) from 4a. Logic unchanged; paths parameterised.
set -euo pipefail

# ---------------- CONFIG ----------------
input_dir="${input_dir:-path/to/split_replicate_files}"             # dir holding r00/r01/r02
nanopolish_file="${nanopolish_file:-nanopolish_file_path/sample_name_nano.tsv}"
sample_name="${sample_name:-sample_name}"                            # output prefix
# ----------------------------------------
cd "${input_dir}"

for rfile in r00 r01 r02; do
  # keep header (line 1) then grep records whose read id is in this replicate group
  (head -n 1 "${nanopolish_file}" && grep -F -f "${rfile}" -w "${nanopolish_file}") \
    > "./${sample_name}_split${rfile}.tsv"
done
