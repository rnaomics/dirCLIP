#!/usr/bin/env bash
# 6_run_nanocompore.sh: Nanocompore sampcomp (GMM signal-deviation calling)
# Brushed: file lists/paths/labels parameterised. Statistical settings unchanged.
set -euo pipefail

# ---------------- CONFIG ----------------
NANOCOMPORE="${NANOCOMPORE:-nanocompore}"            # binary (or rely on $PATH)
# Control (condition 1) NanopolishComp Eventalign_collapse outputs:
noUV1="${noUV1:-noUV_tsp_filtered_split_r00/out_eventalign_collapse.tsv}"
noUV2="${noUV2:-noUV_tsp_filtered_split_r01/out_eventalign_collapse.tsv}"
noUV3="${noUV3:-noUV_tsp_filtered_split_r02/out_eventalign_collapse.tsv}"
# Test (condition 2) outputs:
SRSF3_1="${SRSF3_1:-SRSF3_rep1/out_eventalign_collapse.tsv}"
SRSF3_2="${SRSF3_2:-SRSF3_rep2/out_eventalign_collapse.tsv}"
SRSF3_3="${SRSF3_3:-SRSF3_rep3/out_eventalign_collapse.tsv}"
LABEL1="${LABEL1:-noUV_tsp}"
LABEL2="${LABEL2:-SRSF3}"
FASTA="${FASTA:-Mus_musculus.GRCm39.cdna.all.fa}"    # see README §reference re: coordinate frame
OUTPATH="${OUTPATH:-nanocompore_res/noUV_SRSF3}"
OUTPREFIX="${OUTPREFIX:-noUV_vs_SRSF3}"
MIN_COVERAGE="${MIN_COVERAGE:-5}"
MAX_INVALID_KMERS_FREQ="${MAX_INVALID_KMERS_FREQ:-0.6}"
# ----------------------------------------

# Non-default parameter rationale:
#  --min_coverage 5            : skip sites with <5 supporting reads
#  --comparison_methods GMM    : GMM integrates signal intensity + dwell time
#  --max_invalid_kmers_freq 0.6: permissive (adduct-driven extreme signal -> more invalid kmers)
#  --logit                     : logistic regression of cluster membership on condition (Nanocompore 1.0.x;
#                                comment corrected 10 Oct 2026, see ../nanocompore_pvalue/)
"${NANOCOMPORE}" sampcomp \
  --file_list1 "${noUV1},${noUV2},${noUV3}" \
  --file_list2 "${SRSF3_1},${SRSF3_2},${SRSF3_3}" \
  --label1 "${LABEL1}" \
  --label2 "${LABEL2}" \
  --fasta "${FASTA}" \
  --outpath "${OUTPATH}" \
  --outprefix "${OUTPREFIX}" \
  --min_coverage "${MIN_COVERAGE}" --comparison_methods GMM \
  --max_invalid_kmers_freq "${MAX_INVALID_KMERS_FREQ}" --log_level debug --logit --overwrite
