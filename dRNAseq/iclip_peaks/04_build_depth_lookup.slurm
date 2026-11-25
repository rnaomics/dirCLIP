#!/bin/bash
#SBATCH --job-name=a4a_lookup
# no partition is named here: sbatch -p <partition>, or export SBATCH_PARTITION
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/a4a_%j.out
#SBATCH --error=logs/a4a_%j.err
#
# A4a step 2: build a (transcript, transcript position) -> iCLIP overlap lookup.
#
# Why this exists. The corrected R3.2 run already lifted all 344,464 DRS positions
# from transcript space into GRCm39 and intersected them against the mm9 to mm39
# lifted SRSF3 iCLIP peaks. The awk formatter in 02_liftover_and_intersect.sh packs
# the transcript ID but drops the transcript POSITION, so the existing
# sites_labeled.*.bed files cannot be joined to any other run's position list.
#
# This job repeats only the formatting and the intersect, reusing the same
# sites_genomic.raw.tsv and the same lifted iCLIP BEDs, and keeps the transcript
# position. No liftover is repeated, so the mm9 to mm39 correction is carried
# through verbatim and no new coordinate error can be introduced.
#
# Output: a4a_lookup.tsv.gz, one row per DRS position, carrying our own GMM p and
# LOR alongside the six iCLIP overlap counts (full and centre, at the all, hc15
# and hc79 peak confidence tiers).
set +u
# CLUSTER PATH. CONDA_SH is the etc/profile.d/conda.sh of your conda installation; CONDA_ENV is the
# environment that holds the tools of the top README, section 2.
source "${CONDA_SH:?set CONDA_SH to the conda.sh of your conda installation}"
conda activate "${CONDA_ENV:?set CONDA_ENV to the conda environment with the tools of the top README}"
set -euo pipefail

# CLUSTER PATH. Set B to the R3.2 iCLIP signal working directory on your system.
B="${B:?set B to the R3.2 iCLIP signal working directory}"
W=$B/a4a_lookup
mkdir -p "$W" logs
cd "$W"

RAW=$B/sites_genomic.raw.tsv

# ---- 1. reformat, keeping the transcript position ($8) in the packed name ------
# name = transcript|tx_pos|ref_kmer|GMM_logit_pvalue|Logit_LOR|abs_d_intensity|abs_log2_dwell_fc|test_depth
echo "[$(date)] format"
tail -n +2 "$RAW" | awk -F'\t' 'BEGIN{OFS="\t"}
{
  chr=$1; if(chr=="chrM"||chr=="chrMT"||chr=="M"){chr="MT"} else {sub(/^chr/,"",chr)}
  name=$7"|"$8"|"$11"|"$12"|"$13"|"$16"|"$19"|"$20
  print chr, $2, $3, name, ".", $6
}' | sort -k1,1 -k2,2n > sites_pos.sorted.bed
echo "  rows: $(wc -l < sites_pos.sorted.bed)"

# ---- 2. intersect against each peak tier, keeping only the count column --------
for V in full centre; do
  for T in "sorted:all" "hc15.sorted:hc15" "hc79.sorted:hc79"; do
    SUF=${T%%:*}; TAG=${T##*:}
    echo "[$(date)] intersect ${V} ${TAG}"
    bedtools intersect -s -c -a sites_pos.sorted.bed -b "$B/iclip.${V}.${SUF}.bed" \
      | awk -F'\t' '{print $NF}' > "count_${V}_${TAG}.txt"
  done
done

# ---- 3. assemble one compact table --------------------------------------------
echo "[$(date)] assemble"
awk -F'\t' 'BEGIN{OFS="\t"}{n=split($4,a,"|"); print a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]}' \
  sites_pos.sorted.bed > payload.tsv

{
  printf 'transcript\ttx_pos\tref_kmer\tour_GMM_p\tour_LOR\tabs_d_intensity\tabs_log2_dwell_fc\ttest_depth'
  printf '\tov_full_all\tov_full_hc15\tov_full_hc79\tov_centre_all\tov_centre_hc15\tov_centre_hc79\n'
  paste payload.tsv count_full_all.txt count_full_hc15.txt count_full_hc79.txt \
                    count_centre_all.txt count_centre_hc15.txt count_centre_hc79.txt
} | gzip > a4a_lookup.tsv.gz

echo "[$(date)] DONE -> $W/a4a_lookup.tsv.gz"
zcat a4a_lookup.tsv.gz | wc -l
