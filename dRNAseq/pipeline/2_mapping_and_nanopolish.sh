#!/usr/bin/env bash
# 2_mapping_and_nanopolish.sh: minimap2 mapping + nanopolish eventalign (direct-RNA)
# Brushed for reproducibility: paths/tools/threads are variables; logic unchanged.
set -euo pipefail

# ---------------- CONFIG (edit or override via env) ----------------
BAM_NAME="${BAM_NAME:-sample_bam_file}"                 # output BAM basename
FASTQ="${FASTQ:-merged_all_reads.fastq.gz}"             # reads from step 1
REF="${REF:-Mus_musculus.GRCm39.combined.fa}"          # reference transcriptome (see README §reference)
SEQ_SUM="${SEQ_SUM:-sequencing_summary.txt}"            # MinKNOW/guppy sequencing_summary
FAST5_DIR="${FAST5_DIR:-path_to_fast5/fast5}"
NANO_OUT="${NANO_OUT:-nanopolish_file_path/sample_name_nano.tsv}"
NANOPOLISH="${NANOPOLISH:-nanopolish}"                  # nanopolish binary (or rely on $PATH)
MAP_THREADS="${MAP_THREADS:-48}"
EVENTALIGN_THREADS="${EVENTALIGN_THREADS:-30}"
# -------------------------------------------------------------------

# map to reference transcriptome (noisy Nanopore direct RNA-seq)
# -ax map-ont : optimised for dRNA; mapping to transcriptome so no splice mode needed
# -k14        : minimum k-mer length 14, better for dRNA
minimap2 -ax map-ont -k14 -uf "${REF}" "${FASTQ}" > "${BAM_NAME}.sam"

# remove unmapped (-F 4) / keep primary alignments (-F 2324) / sort to BAM
samtools view -b -@ "${MAP_THREADS}" -F 4 -u "${BAM_NAME}.sam" \
  | samtools view -@ "${MAP_THREADS}" -b -F 2324 \
  | samtools sort -@ "${MAP_THREADS}" > "${BAM_NAME}.sorted.bam"
samtools index "${BAM_NAME}.sorted.bam"

# index reads against signal (required for signal-level analysis)
"${NANOPOLISH}" index -s "${SEQ_SUM}" -d "${FAST5_DIR}" "${FASTQ}"

# eventalign: align raw signal to reference -> per-read signal-level alignment
# --scale-events normalise signal; --print-read-names keep IDs (needed for replicate splitting);
# --samples write raw samples into the TSV
"${NANOPOLISH}" eventalign -t "${EVENTALIGN_THREADS}" --reads "${FASTQ}" \
  --bam "${BAM_NAME}.sorted.bam" --genome "${REF}" \
  --scale-events --print-read-names --samples > "${NANO_OUT}"
