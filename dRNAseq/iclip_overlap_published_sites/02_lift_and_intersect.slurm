#!/bin/bash
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
#SBATCH --job-name=r8_iclip
# no partition is named here: sbatch -p <partition>, or export SBATCH_PARTITION
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=logs/r8_%j.out
#SBATCH --error=logs/r8_%j.err
#
# Step 2 (cluster). Lift every tested position of the published Nanocompore sheets to GRCm39,
# widen to the 11 nt window of the coauthor's R script (position +/- 5 nt, genomic), and intersect on
# the same strand with the lifted SRSF3 iCLIP peaks of 6 Jul 2026 (all, score >= 15, score >= 79).
# Reuses those peak files and the GTF unchanged. No liftOver is repeated.
set +u
# CLUSTER PATH. CONDA_SH is the etc/profile.d/conda.sh of your conda installation; CONDA_ENV is the
# environment that holds the tools of the top README, section 2.
source "${CONDA_SH:?set CONDA_SH to the conda.sh of your conda installation}"
conda activate "${CONDA_ENV:?set CONDA_ENV to the conda environment with the tools of the top README}"
set -euo pipefail

BASE="${BASE:?set BASE to the pipeline working directory}"
W=$BASE/analysis/R8_published_iclip
P=$BASE/analysis/R3.2_iCLIP_signal
R2D="${R2D:-$BASE/tools/R2Dtool/target/release/r2d}"
GTF="${GTF:-$BASE/ref/Mus_musculus.GRCm39.114.gtf}"
IN=$W/kat_positions_r2d_input.tsv
cd "$W"

echo "[$(date)] versions"; "$R2D" --version; bedtools --version
md5sum "$IN" "$P/iclip.full.sorted.bed" "$P/iclip.full.hc15.sorted.bed" "$P/iclip.full.hc79.sorted.bed"
echo "input rows (without header): $(( $(wc -l < "$IN") - 1 ))"

# ---- 1. liftover (no -t: match on version stripped identifiers, as in August) ----
echo "[$(date)] R2Dtool liftover"
"$R2D" liftover -H -g "$GTF" -i "$IN" -o kat_positions_genomic.raw.tsv
echo "lifted rows (without header): $(( $(wc -l < kat_positions_genomic.raw.tsv) - 1 ))"
head -2 kat_positions_genomic.raw.tsv

# ---- 2. 11 nt window, BED6, name = arm|transcript|pos|LOR|p|gpos ----
# raw columns: 1 chromosome 2 start 3 end 4 name 5 score 6 strand 7 transcript 8 start 9 end
#              10 strand 11 arm 12 pos 13 LOR 14 GMM_adj_pval
tail -n +2 kat_positions_genomic.raw.tsv | awk -F'\t' 'BEGIN{OFS="\t"}
$1=="" || $1=="NA" || $2=="" || $2=="NA" {bad++; next}
{
  chr=$1; if(chr=="chrM"||chr=="chrMT"||chr=="M"){chr="MT"} else {sub(/^chr/,"",chr)}
  s=$2-5; if(s<0){s=0}
  print chr, s, $2+6, $11"|"$7"|"$12"|"$13"|"$14"|"$2, ".", $6
}
END{print "rows without a genomic position: " bad+0 > "/dev/stderr"}' | sort -k1,1 -k2,2n > kat_windows.sorted.bed
echo "window rows: $(wc -l < kat_windows.sorted.bed)"

# ---- 3. counts per tier ----
for T in "sorted:all" "hc15.sorted:hc15" "hc79.sorted:hc79"; do
  SUF=${T%%:*}; TAG=${T##*:}
  echo "[$(date)] intersect -s -c ${TAG}"
  bedtools intersect -s -c -a kat_windows.sorted.bed -b "$P/iclip.full.${SUF}.bed" | awk -F'\t' '{print $NF}' > "count_${TAG}.txt"
done
{
  printf 'arm\ttranscript\tpos\tLOR\tGMM_adj_pval\tchr\tgpos0\tstrand\tn_all\tn_hc15\tn_hc79\n'
  awk -F'\t' 'BEGIN{OFS="\t"}{split($4,a,"|"); print a[1],a[2],a[3],a[4],a[5],$1,a[6],$6}' kat_windows.sorted.bed \
    | paste - count_all.txt count_hc15.txt count_hc79.txt
} | gzip > r8_positions.tsv.gz

# ---- 4. position to peak pairs (all peaks, with score), for peak level recovery ----
echo "[$(date)] intersect -s -wa -wb all"
bedtools intersect -s -wa -wb -a kat_windows.sorted.bed -b "$P/iclip.full.sorted.bed" \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "arm","transcript","pos","peak","score"}
    {split($4,a,"|"); print a[1],a[2],a[3],$7":"$8"-"$9":"$12,$11}' | gzip > r8_pairs.tsv.gz

# ---- 5. transcript versions: table against GTF ----
echo "[$(date)] transcript versions"
tail -n +2 "$IN" | cut -f1 | sort -u > kat_transcripts.txt
awk -F'\t' '$3=="transcript"' "$GTF" | grep -o 'transcript_id "[^"]*"; transcript_version "[^"]*"' \
  | sed 's/transcript_id "//; s/"; transcript_version "/\t/; s/"$//' | sort -u > gtf_tx_versions.tsv
awk -F'\t' 'NR==FNR{v[$1]=$2; next}
  {split($1,a,"."); g=(a[1] in v)?v[a[1]]:"absent"; st=(g=="absent")?"absent":((g==a[2])?"same":"differs"); print $1"\t"g"\t"st}' \
  gtf_tx_versions.tsv kat_transcripts.txt > r8_transcript_versions.tsv
cut -f3 r8_transcript_versions.tsv | sort | uniq -c

echo "[$(date)] sanity"
zcat r8_positions.tsv.gz | awk -F'\t' 'NR>1{n[$1]++; if($9>0) a[$1]++} END{for(k in n) print k, n[k], a[k], a[k]/n[k]*100}'
zcat r8_pairs.tsv.gz | wc -l
ls -la
echo "[$(date)] DONE"
