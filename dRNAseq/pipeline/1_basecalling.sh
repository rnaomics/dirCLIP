#!/usr/bin/env bash
# 1_basecalling.sh: Guppy basecalling of direct-RNA FAST5 (FLO-MIN106 / SQK-RNA002)
# Brushed for reproducibility: all paths are variables, overridable from the environment.
# Scientific behaviour unchanged from the original.
set -euo pipefail

# ---------------- CONFIG (edit or override via env) ----------------
GUPPY="${GUPPY:-guppy_basecaller}"          # path to guppy_basecaller binary (or rely on $PATH)
FAST5_DIR="${FAST5_DIR:-path_to_fast5s/fast5}"   # directory containing FAST5 files
FASTQ_DIR="${FASTQ_DIR:-path_to_fastq/fastq}"    # output directory for FASTQ
FLOWCELL="${FLOWCELL:-FLO-MIN106}"
KIT="${KIT:-SQK-RNA002}"
DEVICE="${DEVICE:-cuda:all:100%}"
NUM_CALLERS="${NUM_CALLERS:-4}"
RUNNERS_PER_DEVICE="${RUNNERS_PER_DEVICE:-2}"
CHUNKS_PER_RUNNER="${CHUNKS_PER_RUNNER:-512}"
CHUNK_SIZE="${CHUNK_SIZE:-1000}"
# -------------------------------------------------------------------

# basecall with guppy
"${GUPPY}" --input_path "${FAST5_DIR}" \
  --save_path "${FASTQ_DIR}" --flowcell "${FLOWCELL}" --kit "${KIT}" \
  --device "${DEVICE}" --compress_fastq --recursive \
  --num_callers "${NUM_CALLERS}" \
  --gpu_runners_per_device "${RUNNERS_PER_DEVICE}" \
  --chunks_per_runner "${CHUNKS_PER_RUNNER}" \
  --chunk_size "${CHUNK_SIZE}"

# merge pass + fail FASTQ into a single file containing all reads
cat "${FASTQ_DIR}"/pass/*.fastq.gz "${FASTQ_DIR}"/fail/*.fastq.gz \
  > "${FASTQ_DIR}/merged_all_reads.fastq.gz"
