#!/bin/bash
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
# R3.2 orthogonal test, step 2 of 3 (run on the cluster).
#
# 1. R2Dtool liftover: transcript-coordinate Nanocompore sites -> genomic BED6+payload (GRCm39).
#    (skipped if sites_genomic.sorted.bed already present; the DRS side is assembly-correct.)
# 2. Cross-assembly liftOver of iCount iCLIP peaks mm9 (NCBI37) -> mm39 (GRCm39)
#    [UCSC liftOver + mm9ToMm39 chain], THEN harmonise UCSC 'chr1' -> Ensembl '1' (chrM->MT),
#    for BOTH peak-handling variants: CENTRE (peak centre, 1 nt) and FULL (whole peak).
# 3. Stranded bedtools intersect to label each DRS site as iCLIP-overlapping vs not.
#
# ASSEMBLY FIX (2026-07-06, coauthor correction): the iCLIP peaks are mm9-based;
# earlier runs intersected them against mm39 DRS coordinates with only chromosome-name
# harmonisation and NO cross-assembly liftOver. This version lifts the peaks mm9->mm39 first.
#
# R2Dtool output column layout (verified):
#   1 chromosome  2 start  3 end  4 name  5 score  6 strand   <- genomic BED6 (prepended)
#   7 transcript  8 start  9 end  10 strand  11 ref_kmer  12 GMM_logit_pvalue
#   13 Logit_LOR  14 c1_mean_intensity  15 c2_mean_intensity  16 abs_d_intensity
#   17 c1_mean_dwell  18 c2_mean_dwell   19 abs_log2_dwell_fc  20 test_depth
set -euo pipefail

BASE="${BASE:?set BASE to the pipeline working directory}"
R2D="${R2D:-$BASE/tools/R2Dtool/target/release/r2d}"
GTF="${GTF:-$BASE/ref/Mus_musculus.GRCm39.114.gtf}"
ICLIP_BED="${ICLIP_BED:-$BASE/analysis/R3.2_iCLIP_signal/SRSF3_merged.3nt_3nt.peaks.bed.gz}"
LIFTOVER="${LIFTOVER:-$BASE/tools/bin/liftOver}"
CHAIN="${CHAIN:-$BASE/ref/chains/mm9ToMm39.over.chain.gz}"
INPUT="${INPUT:-r2dtool_input.tsv}"
OUT="${OUT:-.}"
mkdir -p "$OUT"

# ---- 1. R2Dtool liftover (DRS side, GRCm39): reuse if already present ---------
if [[ -s "$OUT/sites_genomic.sorted.bed" ]]; then
  echo "[$(date)] reuse existing sites_genomic.sorted.bed ($(wc -l < "$OUT/sites_genomic.sorted.bed") sites)"
else
  echo "[$(date)] R2Dtool liftover"
  "$R2D" liftover -H -g "$GTF" -i "$INPUT" -o "$OUT/sites_genomic.raw.tsv"
  echo "[$(date)] format lifted sites"
  tail -n +2 "$OUT/sites_genomic.raw.tsv" | awk -F'\t' 'BEGIN{OFS="\t"}
  {
    chr=$1; if(chr=="chrM"||chr=="chrMT"||chr=="M"){chr="MT"} else {sub(/^chr/,"",chr)}
    name=$7"|"$11"|"$12"|"$13"|"$14"|"$15"|"$16"|"$17"|"$18"|"$19"|"$20
    print chr, $2, $3, name, ".", $6
  }' | sort -k1,1 -k2,2n > "$OUT/sites_genomic.sorted.bed"
  echo "  lifted sites: $(wc -l < "$OUT/sites_genomic.sorted.bed")"
fi

# ---- 2a. iCLIP peaks -> plain BED6, KEEP UCSC 'chr' naming for liftOver ---------
echo "[$(date)] iCLIP peaks -> BED6 (mm9, UCSC naming)"
zcat "$ICLIP_BED" | awk -F'\t' 'BEGIN{OFS="\t"}
/^track/||/^#/{next}
{ print $1,$2,$3,($4==""?".":$4),($5==""?".":$5),($6==""?"+":$6) }' > "$OUT/iclip.mm9.bed"
echo "  peaks (mm9): $(wc -l < "$OUT/iclip.mm9.bed")"

# ---- 2b. cross-assembly liftOver mm9 -> mm39 -----------------------------------
echo "[$(date)] liftOver mm9 -> mm39"
"$LIFTOVER" "$OUT/iclip.mm9.bed" "$CHAIN" "$OUT/iclip.mm39.bed" "$OUT/iclip.mm9_unmapped.bed"
mapped=$(wc -l < "$OUT/iclip.mm39.bed")
unmapped=$(grep -vc '^#' "$OUT/iclip.mm9_unmapped.bed" 2>/dev/null || echo 0)
echo "  mapped: $mapped   unmapped: $unmapped"

# ---- 2c. mm39 UCSC -> Ensembl naming; CENTRE and FULL variants -----------------
echo "[$(date)] harmonise lifted peaks -> Ensembl naming"
awk -F'\t' -v full="$OUT/iclip.full.bed" -v centre="$OUT/iclip.centre.bed" 'BEGIN{OFS="\t"}
{
  chr=$1; if(chr=="chrM"||chr=="chrMT"||chr=="M"){chr="MT"} else {sub(/^chr/,"",chr)}
  s=$2; e=$3; name=$4; score=$5; str=$6
  print chr,s,e,name,score,str       > full
  mid=int((s+e)/2)
  print chr,mid,mid+1,name,score,str > centre
}' "$OUT/iclip.mm39.bed"
sort -k1,1 -k2,2n "$OUT/iclip.full.bed"   > "$OUT/iclip.full.sorted.bed"
sort -k1,1 -k2,2n "$OUT/iclip.centre.bed" > "$OUT/iclip.centre.sorted.bed"

# ---- 3. stranded intersect: count iCLIP overlaps per DRS site ------------------
echo "[$(date)] bedtools intersect"
for V in centre full; do
  bedtools intersect -s -c \
    -a "$OUT/sites_genomic.sorted.bed" \
    -b "$OUT/iclip.${V}.sorted.bed" \
    > "$OUT/sites_labeled.${V}.bed"   # last col = # iCLIP overlaps (0 = background)
  fg=$(awk -F'\t' '$NF>0' "$OUT/sites_labeled.${V}.bed" | wc -l)
  echo "  [$V] foreground sites (iCLIP-overlapping): $fg"
done
echo "[$(date)] DONE -> $OUT/sites_labeled.{centre,full}.bed"
