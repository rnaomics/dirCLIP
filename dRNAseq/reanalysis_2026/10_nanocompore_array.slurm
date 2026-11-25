#!/bin/bash
#SBATCH --job-name=dirclip_nc
# no partition is named here: sbatch -p <partition>, or export SBATCH_PARTITION
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/nc_%A_%a.out
#SBATCH --error=logs/nc_%A_%a.err
#SBATCH --array=0-1
#
# dirCLIP signal deviation site calling (pipeline script 6: nanocompore sampcomp, GMM).
# Statistical settings unchanged from the original pipeline:
#   --min_coverage 5 --comparison_methods GMM --max_invalid_kmers_freq 0.6 --logit
#
# Comparison design (CONTROL = file_list1, TEST = file_list2):
#   [0] SRSF3 : control = SRSF3 -UV total-RNA transcriptome (rep1,rep2)
#               test    = SRSF3 dirCLIP (rep1,rep2)            ref = cdna.all
#   [1] HMGA1a: control = WT non-UV iv-polyA transcriptome (rep1)
#               test    = HMGA1a dirCLIP iv-polyA (rep1,rep2,rep3) ref = combined
# ** Control choice is the one scientific decision here; see the README. **

# environment FIRST without strict mode (conda activate trips set -u)
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
OUT=$BASE/analysis/nanocompore
mkdir -p "$OUT" logs

CDNA=$REFDIR/Mus_musculus.GRCm39.cdna.all.fa
COMB=$REFDIR/Mus_musculus.GRCm39.combined.fa
C() { echo "$EA/$1/collapsed/out_eventalign_collapse.tsv"; }

case "$SLURM_ARRAY_TASK_ID" in
  0)
    LABEL1=SRSF3_noUV; LABEL2=SRSF3_CLIP; FASTA=$CDNA
    G1=(SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep1 SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2)
    G2=(SRSF3_DIRCLIP_dRNA_rep1 SRSF3_DIRCLIP_dRNA_rep2)
    OUTPATH=$OUT/SRSF3_noUV_vs_CLIP; PREFIX=SRSF3_noUV_vs_CLIP ;;
  1)
    LABEL1=HMGA1a_WTivpolyA; LABEL2=HMGA1a_CLIP; FASTA=$COMB
    G1=(WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1)
    G2=(HMGA1a_DIRCLIP_iv-polyA_dRNA_rep1 HMGA1a_DIRCLIP_iv-polyA_dRNA_rep2 HMGA1a_DIRCLIP_iv-polyA_dRNA_rep3)
    OUTPATH=$OUT/HMGA1a_WTivpolyA_vs_CLIP; PREFIX=HMGA1a_WTivpolyA_vs_CLIP ;;
esac

echo "[$(date)] nanocompore $LABEL1 vs $LABEL2  node=$(hostname)"
mkdir -p "$OUTPATH"

# Build comma lists from only the replicates whose collapsed eventalign exists.
# This tolerates a missing/incomplete per-sample upload (e.g. a truncated FAST5
# zip) by dropping that replicate. Control/test DESIGN and all GMM statistical
# settings are unchanged; only data-absent replicates are omitted, and every
# inclusion/exclusion is logged so the realised design is explicit in the log.
build_list() {
  local out="" s f
  for s in "$@"; do
    f=$(C "$s")
    if [[ -s "$f" ]]; then out="${out:+$out,}$f"; echo "  include: $s" >&2
    else echo "  SKIP (no collapsed eventalign): $s" >&2; fi
  done
  printf '%s' "$out"
}
echo "[group1 control = $LABEL1]" >&2; L1=$(build_list "${G1[@]}")
echo "[group2 test    = $LABEL2]" >&2; L2=$(build_list "${G2[@]}")
if [[ -z "$L1" || -z "$L2" ]]; then
  echo "FATAL: control or test group has no usable replicates; aborting $PREFIX" >&2
  exit 1
fi

nanocompore sampcomp \
  --file_list1 "$L1" \
  --file_list2 "$L2" \
  --label1 "$LABEL1" --label2 "$LABEL2" \
  --fasta "$FASTA" \
  --outpath "$OUTPATH" --outprefix "$PREFIX" \
  --nthreads "${SLURM_CPUS_PER_TASK:-16}" \
  --min_coverage 5 --comparison_methods GMM \
  --max_invalid_kmers_freq 0.6 --logit --overwrite --log_level info

echo "[$(date)] DONE nanocompore $PREFIX"
