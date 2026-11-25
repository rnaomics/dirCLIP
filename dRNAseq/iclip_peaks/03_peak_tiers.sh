#!/bin/bash
# 03_peak_tiers.sh: the iCLIP peak sets by crosslink score, as made on 6 Jul 2026 inside
# run_R3.2_highconf.slurm on the cluster (the awk filter below is that script's, unchanged).
# Tiers: all peaks; score >= 15 (about the top 10 percent); score >= 79 (about the top 1 percent).
# Input: iclip.full.sorted.bed and iclip.centre.sorted.bed from 02_liftover_and_intersect.sh.
set -euo pipefail
OUT="${OUT:-.}"
for TIER in 15 79; do
  for V in full centre; do
    awk -F'\t' -v c=$TIER 'BEGIN{OFS="\t"} $5>=c' "$OUT/iclip.$V.sorted.bed" > "$OUT/iclip.$V.hc$TIER.sorted.bed"
    echo "  tier>=$TIER $V : peaks=$(wc -l < "$OUT/iclip.$V.hc$TIER.sorted.bed")"
  done
done
