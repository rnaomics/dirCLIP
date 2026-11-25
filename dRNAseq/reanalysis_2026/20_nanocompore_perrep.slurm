#!/bin/bash
#SBATCH --job-name=dirclip_ncpr
# no partition is named here: sbatch -p <partition>, or export SBATCH_PARTITION
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/ncpr_%A_%a.out
#SBATCH --error=logs/ncpr_%A_%a.err
#SBATCH --array=0-4
#
# R3.11 site-level reproducibility: call signal-deviation sites with EACH dirCLIP test
# replicate SEPARATELY (against the SAME pooled control group), so per-replicate site
# sets can be overlapped (UpSet + Jaccard). Statistical settings are IDENTICAL to the
# group-level run (10_nanocompore_array.slurm): --min_coverage 5 --comparison_methods GMM
# --max_invalid_kmers_freq 0.6 --logit. Only the test group membership changes (1 rep).
#
# Tasks:
#   0 SRSF3  test=CLIP_rep1   1 SRSF3  test=CLIP_rep2
#   2 HMGA1a test=CLIP_rep1   3 HMGA1a test=CLIP_rep2   4 HMGA1a test=CLIP_rep3

set +u
# CLUSTER PATH. CONDA_SH is the etc/profile.d/conda.sh of your conda installation; CONDA_ENV is the
# environment that holds the tools of the top README, section 2.
source "${CONDA_SH:?set CONDA_SH to the conda.sh of your conda installation}"
conda activate "${CONDA_ENV:?set CONDA_ENV to the conda environment with the tools of the top README}"
set -euo pipefail

# CLUSTER PATH. Set BASE to the pipeline working directory on your own system.
BASE="${BASE:?set BASE to the pipeline working directory}"
EA=$BASE/analysis/eventalign
REFDIR=$BASE/ref
OUT=$BASE/analysis/nanocompore_perrep
mkdir -p "$OUT" logs
CDNA=$REFDIR/Mus_musculus.GRCm39.cdna.all.fa
COMB=$REFDIR/Mus_musculus.GRCm39.combined.fa
C() { echo "$EA/$1/collapsed/out_eventalign_collapse.tsv"; }

case "$SLURM_ARRAY_TASK_ID" in
  0) RBP=SRSF3;  FASTA=$CDNA; L1=SRSF3_noUV; L2=SRSF3_CLIP_rep1
     G1=(SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep1 SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2)
     G2=(SRSF3_DIRCLIP_dRNA_rep1) ;;
  1) RBP=SRSF3;  FASTA=$CDNA; L1=SRSF3_noUV; L2=SRSF3_CLIP_rep2
     G1=(SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep1 SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2)
     G2=(SRSF3_DIRCLIP_dRNA_rep2) ;;
  2) RBP=HMGA1a; FASTA=$COMB; L1=HMGA1a_WTivpolyA; L2=HMGA1a_CLIP_rep1
     G1=(WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1); G2=(HMGA1a_DIRCLIP_iv-polyA_dRNA_rep1) ;;
  3) RBP=HMGA1a; FASTA=$COMB; L1=HMGA1a_WTivpolyA; L2=HMGA1a_CLIP_rep2
     G1=(WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1); G2=(HMGA1a_DIRCLIP_iv-polyA_dRNA_rep2) ;;
  4) RBP=HMGA1a; FASTA=$COMB; L1=HMGA1a_WTivpolyA; L2=HMGA1a_CLIP_rep3
     G1=(WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1); G2=(HMGA1a_DIRCLIP_iv-polyA_dRNA_rep3) ;;
esac

PREFIX=${RBP}_${L2}_vs_${L1}
OUTPATH=$OUT/$PREFIX
mkdir -p "$OUTPATH"
echo "[$(date)] per-rep nanocompore $PREFIX  node=$(hostname)"

build_list() { local out="" s f; for s in "$@"; do f=$(C "$s")
  if [[ -s "$f" ]]; then out="${out:+$out,}$f"; echo "  include: $s" >&2
  else echo "  SKIP (no collapsed): $s" >&2; fi; done; printf '%s' "$out"; }
echo "[ctrl=$L1]" >&2; LST1=$(build_list "${G1[@]}")
echo "[test=$L2]" >&2; LST2=$(build_list "${G2[@]}")
[[ -z "$LST1" || -z "$LST2" ]] && { echo "FATAL no reps $PREFIX" >&2; exit 1; }

nanocompore sampcomp \
  --file_list1 "$LST1" --file_list2 "$LST2" \
  --label1 "$L1" --label2 "$L2" --fasta "$FASTA" \
  --outpath "$OUTPATH" --outprefix "$PREFIX" \
  --nthreads "${SLURM_CPUS_PER_TASK:-16}" \
  --min_coverage 5 --comparison_methods GMM \
  --max_invalid_kmers_freq 0.6 --logit --overwrite --log_level info
echo "[$(date)] DONE $PREFIX"
