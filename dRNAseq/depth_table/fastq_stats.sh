#!/usr/bin/env bash
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
# R1.3 / R3.3 sequencing-depth table.
# Streams each gzipped FASTQ member straight out of fastq.zip (no extraction)
# and computes: reads, total bases, mean/median length, N50, max length.
set -euo pipefail

ZIP="${ZIP:-fastq.zip}"      # zip archive holding one .fastq.gz per library
OUT="${OUT:-.}"
mkdir -p "$OUT/per_sample"

# awk: read every 2nd line of each 4-line FASTQ record (the sequence), tally lengths
read -r -d '' AWKPROG <<'EOF' || true
NR%4==2 {
  L=length($0); n++; tot+=L; len[n]=L; if(L>mx)mx=L
}
END{
  if(n==0){print "0\t0\t0\t0\t0\t0"; exit}
  asort(len)
  # median
  if(n%2){med=len[(n+1)/2]} else {med=(len[n/2]+len[n/2+1])/2}
  # N50
  half=tot/2; acc=0; n50=0
  for(i=n;i>=1;i--){acc+=len[i]; if(acc>=half){n50=len[i]; break}}
  printf "%d\t%d\t%.1f\t%d\t%d\t%d\n", n, tot, tot/n, med, n50, mx
}
EOF

stat_one () {
  local member="$1"
  local name
  name=$(basename "$member" .fastq.gz)
  unzip -p "$ZIP" "$member" | gunzip -c | gawk "$AWKPROG" > "$OUT/per_sample/${name}.tsv"
  echo "done: $name"
}
export -f stat_one
export ZIP OUT AWKPROG

MEMBERS=$(unzip -Z1 "$ZIP" | grep '\.fastq\.gz$')
echo "$MEMBERS" | xargs -I{} -P 9 bash -c 'stat_one "$@"' _ {}

echo "=== combining ==="
printf "sample\treads\ttotal_bases\tmean_len\tmedian_len\tN50\tmax_len\n" > "$OUT/depth_table.tsv"
for f in "$OUT"/per_sample/*.tsv; do
  s=$(basename "$f" .tsv)
  printf "%s\t%s\n" "$s" "$(cat "$f")" >> "$OUT/depth_table.tsv"
done
column -t -s$'\t' "$OUT/depth_table.tsv"
