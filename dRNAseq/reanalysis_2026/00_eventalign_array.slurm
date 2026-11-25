#!/bin/bash
#SBATCH --job-name=dirclip_ea
# no partition is named here: sbatch -p <partition>, or export SBATCH_PARTITION
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/ea_%A_%a.out
#SBATCH --error=logs/ea_%A_%a.err
#SBATCH --array=0-8
#
# dirCLIP FAST5 signal pipeline, per-sample stage (deterministic; needed for every sample).
# Faithful to pipeline scripts 2 and 5: nanopolish index -> nanopolish eventalign
#   (--scale-events --print-read-names --samples) -> NanopolishComp Eventalign_collapse.
# Engine = nanopolish 0.14.0 (NOT f5c) to match the original published pipeline.
# Extraction uses python zipfile because the FAST5 zips are >4 GB streamed ZIP64
# archives that system UnZip 6.00 cannot read.
#
# Each step is idempotent: a sample that already has out_eventalign_collapse.tsv is skipped.

# --- environment FIRST, WITHOUT strict mode (profile.d/conda reference unset vars; set -u
#     would abort the job at conda activate) ---
set +u
source /etc/profile 2>/dev/null || true
export COMPILER=gcc; export COMPILER_VER=14.3.0
module load samtools/1.21 2>/dev/null || true
# CLUSTER PATH. CONDA_SH is the etc/profile.d/conda.sh of your conda installation; CONDA_ENV is the
# environment that holds the tools of the top README, section 2.
source "${CONDA_SH:?set CONDA_SH to the conda.sh of your conda installation}"
conda activate "${CONDA_ENV:?set CONDA_ENV to the conda environment with the tools of the top README}"
# now safe to be strict
set -euo pipefail

# CLUSTER PATH. Set BASE to the pipeline working directory on your own system.
BASE="${BASE:?set BASE to the pipeline working directory}"
REFDIR=$BASE/ref
FQDIR=$BASE/fastq
BAMDIR=$BASE/analysis/mapping/bam
F5ZIP=$BASE/fast5
RRZIP=$BASE/Run_reports
WORKROOT=$BASE/analysis/eventalign
mkdir -p "$WORKROOT" logs

CDNA=$REFDIR/Mus_musculus.GRCm39.cdna.all.fa
COMB=$REFDIR/Mus_musculus.GRCm39.combined.fa

# index: "sample reference"  (per-RBP ref, identical to map_dRNA.slurm)
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
BAM=$BAMDIR/${SAMPLE}.sorted.bam
WD=$WORKROOT/$SAMPLE
mkdir -p "$WD"
NPROC=${SLURM_CPUS_PER_TASK:-16}

echo "[$(date)] START $SAMPLE  ref=$(basename "$REF")  cpus=$NPROC  node=$(hostname)"

COLLAPSED=$WD/collapsed/out_eventalign_collapse.tsv
if [[ -s "$COLLAPSED" ]]; then
  echo "[$(date)] $SAMPLE already collapsed -> skip"; exit 0
fi

# ---- 1. sequencing_summary (from the clean Run_reports zip; python for uniformity) ----
SEQSUM=$WD/sequencing_summary.txt
if [[ ! -s "$SEQSUM" ]]; then
  echo "[$(date)] extract sequencing_summary"
  python3 - "$RRZIP/${SAMPLE}.zip" "$SEQSUM" <<'PY'
import zipfile, sys, shutil
zp, out = sys.argv[1], sys.argv[2]
zf = zipfile.ZipFile(zp)
cand = [n for n in zf.namelist()
        if "sequencing_summary" in n and n.endswith(".txt") and "_Error" not in n]
if not cand:
    sys.exit("NO sequencing_summary in %s" % zp)
with zf.open(cand[0]) as fi, open(out, "wb") as fo:
    shutil.copyfileobj(fi, fo)
print("  ->", cand[0])
PY
fi

# ---- 2. extract FAST5 (pass+fail) via python zipfile (system unzip can't read these) ----
F5DIR=$WD/fast5
if [[ ! -f "$WD/.fast5_done" ]]; then
  echo "[$(date)] extract FAST5"
  mkdir -p "$F5DIR"
  python3 - "$F5ZIP/${SAMPLE}.zip" "$F5DIR" <<'PY'
import zipfile, sys, os, shutil
zp, outdir = sys.argv[1], sys.argv[2]
zf = zipfile.ZipFile(zp)
n = 0
for m in zf.namelist():
    if m.endswith(".fast5"):
        dest = os.path.join(outdir, os.path.basename(m))
        with zf.open(m) as fi, open(dest, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1024*1024)
        n += 1
print("  extracted %d fast5" % n)
PY
  touch "$WD/.fast5_done"
fi

# ---- 3. eventalign-ready BAM: drop unmapped/reverse/secondary/supplementary (-F 2324) ----
EABAM=$WD/${SAMPLE}.ea.bam
if [[ ! -s "$EABAM.bai" ]]; then
  echo "[$(date)] filter BAM -F 2324"
  samtools view -b -F 2324 -@ "$NPROC" "$BAM" > "$EABAM"
  samtools index "$EABAM"
fi

# ---- 4. nanopolish index (read_id <-> raw signal) ----
if [[ ! -s "${FQ}.index.readdb" ]]; then
  echo "[$(date)] nanopolish index"
  nanopolish index -s "$SEQSUM" -d "$F5DIR" "$FQ"
fi

# ---- 5. eventalign | NanopolishComp Eventalign_collapse (piped: no giant TSV) ----
echo "[$(date)] eventalign -> collapse"
mkdir -p "$WD/collapsed.tmp"
nanopolish eventalign -t "$NPROC" --reads "$FQ" --bam "$EABAM" --genome "$REF" \
    --scale-events --print-read-names --samples \
  | NanopolishComp Eventalign_collapse -t "$NPROC" -o "$WD/collapsed.tmp"
mv "$WD/collapsed.tmp" "$WD/collapsed"

echo "[$(date)] DONE $SAMPLE -> $COLLAPSED"
