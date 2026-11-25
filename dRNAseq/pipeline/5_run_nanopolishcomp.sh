#!/usr/bin/env bash
# 5_run_nanopolishcomp.sh: NanopolishComp Eventalign_collapse (event-level -> kmer-level)
# Brushed: paths/threads parameterised. Logic unchanged.
# Docs: https://a-slide.github.io/NanopolishComp/demo/Eventalign_collapse_usage/
set -euo pipefail

# ---------------- CONFIG ----------------
NANOPOLISHCOMP="${NANOPOLISHCOMP:-NanopolishComp}"   # binary (or rely on $PATH)
INPUT="${INPUT:-sample_name_split_r00.tsv}"          # nanopolish file from step 2/4b
OUTDIR="${OUTDIR:-path_to_output/sample_name_split_r00}"  # output dir (name preserves replicate)
THREADS="${THREADS:-10}"
# ----------------------------------------

"${NANOPOLISHCOMP}" Eventalign_collapse -i "${INPUT}" -o "${OUTDIR}" -t "${THREADS}"
