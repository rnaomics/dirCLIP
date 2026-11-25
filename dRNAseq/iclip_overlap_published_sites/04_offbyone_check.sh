#!/bin/bash
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
# Diagnostic (cluster). The coauthor's R function load_iclip_peaks() builds IRanges(start = BED start, end = BED end)
# without adding 1 to the 0 based BED start, so every peak is read one nucleotide wider on its low
# coordinate side. Test: does widening each peak by 1 nt on that side reproduce the percentages of that script?
# Equivalent here: extend each 11 nt window by 1 nt on its high coordinate side.
set +u
# CLUSTER PATH. CONDA_SH is the etc/profile.d/conda.sh of your conda installation; CONDA_ENV is the
# environment that holds the tools of the top README, section 2.
source "${CONDA_SH:?set CONDA_SH to the conda.sh of your conda installation}"
conda activate "${CONDA_ENV:?set CONDA_ENV to the conda environment with the tools of the top README}"
set -euo pipefail
BASE="${BASE:?set BASE to the pipeline working directory}"
W=$BASE/analysis/R8_published_iclip
P=$BASE/analysis/R3.2_iCLIP_signal
cd "$W"
tally() {  # stdin: BED6 window + count; prints the 2x2 per arm and threshold
  awk -F'\t' '{split($4,a,"|"); arm=a[1]; lor=a[4]+0; if(lor<0) lor=-lor; p=a[5]+0; hit=($NF>0);
       for(t=1;t<=2;t++){thr=(t==1)?0.005:0.05; s=(p<thr && lor>=1)?"sig":"other"; n[arm,thr,s]++; if(hit) h[arm,thr,s]++}}
     END{split("+UV -UV",A," "); for(i=1;i<=2;i++) for(t=1;t<=2;t++){thr=(t==1)?0.005:0.05; arm=A[i];
       printf "%s p<%s  sig %d/%d = %.2f%%  other %d/%d = %.2f%%\n", arm, thr, h[arm,thr,"sig"], n[arm,thr,"sig"], h[arm,thr,"sig"]/n[arm,thr,"sig"]*100, h[arm,thr,"other"], n[arm,thr,"other"], h[arm,thr,"other"]/n[arm,thr,"other"]*100}}'
}
echo "== BED semantics (peak = [start, end), as computed in step 2)"
bedtools intersect -s -c -a kat_windows.sorted.bed -b "$P/iclip.full.sorted.bed" | tally
echo "== peaks read one nucleotide wider on the low coordinate side (BED start taken as 1 based)"
awk -F'\t' 'BEGIN{OFS="\t"}{$3=$3+1; print}' kat_windows.sorted.bed \
  | bedtools intersect -s -c -a - -b "$P/iclip.full.sorted.bed" | tally
echo "== covered peaks (all tier), BED semantics, per arm"
for ARM in "+UV" "-UV"; do
  n=$(awk -F'\t' -v A="$ARM" '{split($4,a,"|"); if(a[1]==A) print}' kat_windows.sorted.bed | bedtools intersect -s -u -a "$P/iclip.full.sorted.bed" -b - | wc -l)
  m=$(awk -F'\t' -v A="$ARM" 'BEGIN{OFS="\t"}{split($4,a,"|"); if(a[1]==A){$3=$3+1; print}}' kat_windows.sorted.bed | bedtools intersect -s -u -a "$P/iclip.full.sorted.bed" -b - | wc -l)
  echo "$ARM covered peaks: $n ; with the one nucleotide wider reading: $m"
done
