#!/bin/bash
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
#SBATCH --job-name=dirclip_map
# no partition is named here: sbatch -p <partition>, or export SBATCH_PARTITION
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/map_%A_%a.out
#SBATCH --error=logs/map_%A_%a.err
#SBATCH --array=0-8

source /etc/profile 2>/dev/null || true
export COMPILER=gcc; export COMPILER_VER=14.3.0
module load samtools/1.21
set -euo pipefail

BASE="${BASE:?set BASE to the pipeline working directory}"
REFDIR=$BASE/ref
FQDIR=$BASE/fastq
OUT=$BASE/analysis/mapping/bam
MINIMAP2="${MINIMAP2:-$BASE/tools/minimap2-2.29_x64-linux/minimap2}"
mkdir -p "$OUT"

CDNA=$REFDIR/Mus_musculus.GRCm39.cdna.all.fa
COMB=$REFDIR/Mus_musculus.GRCm39.combined.fa

# per-RBP: SRSF3->cdna ; HMGA1a/WT (iv-polyA)->combined cdna+ncrna
SAMPLES=(
"HMGA1a_DIRCLIP_iv-polyA_dRNA_rep1 $COMB"
"HMGA1a_DIRCLIP_iv-polyA_dRNA_rep2 $COMB"
"HMGA1a_DIRCLIP_iv-polyA_dRNA_rep3 $COMB"
"SRSF3_DIRCLIP_dRNA_rep1 $CDNA"
"SRSF3_DIRCLIP_dRNA_rep2 $CDNA"
"SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep1 $CDNA"
"SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2 $CDNA"
"SRSF3-EGFP_BAC_Total_RNA_UV_dRNA_rep1 $CDNA"
"WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1 $COMB"
)

ENTRY=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(echo "$ENTRY" | cut -d" " -f1)
REF=$(echo "$ENTRY" | cut -d" " -f2)
FQ=$FQDIR/${SAMPLE}.fastq.gz

echo "[$(date)] sample=$SAMPLE ref=$(basename "$REF")"
"$MINIMAP2" -t ${SLURM_CPUS_PER_TASK} -ax map-ont -k14 -uf "$REF" "$FQ" \
  | samtools sort -@ 4 -o "$OUT/${SAMPLE}.sorted.bam" -
samtools index "$OUT/${SAMPLE}.sorted.bam"

PRIMARY=$(samtools view -c -F 0x904 "$OUT/${SAMPLE}.sorted.bam")
TOTAL=$(samtools view -c "$OUT/${SAMPLE}.sorted.bam")
echo -e "${SAMPLE}\t${TOTAL}\t${PRIMARY}\t$(basename "$REF")" > "$OUT/${SAMPLE}.count.tsv"
echo "[$(date)] done $SAMPLE total=$TOTAL primary_aligned=$PRIMARY"
