#!/bin/bash
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
# Aggregate per-sample primary-aligned read counts into one table.
set -euo pipefail
BASE="${BASE:?set BASE to the pipeline working directory}"
OUT=$BASE/analysis/mapping
echo -e "sample\ttotal_records\tprimary_aligned\treference" > "$OUT/aligned_read_counts.tsv"
cat "$OUT"/bam/*.count.tsv >> "$OUT/aligned_read_counts.tsv"
column -t "$OUT/aligned_read_counts.tsv"
