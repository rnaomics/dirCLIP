#!/usr/bin/env python3
"""
R3.2 orthogonal test, step 1 of 3.

Build the R2Dtool liftover input from the SRSF3 Nanocompore output, carrying the
signal-deviation metrics (the Fig 1D quantities) as payload columns so they survive
the transcript->genome liftover.

Input  (transcript coordinates, ENSMUST):
  <NANO>/..._nanocompore_results.tsv       (GMM stats)
  <NANO>/..._nanocompore_shift_stats.tsv   (per-position intensity/dwell)

Output (transcriptomic BED-like, with header) for R2Dtool `liftover`:
  cols 1-4: transcript  start  end  strand   (transcript space; strand always '+')
  payload : ref_kmer GMM_logit_pvalue Logit_LOR
            c1_mean_intensity c2_mean_intensity abs_d_intensity
            c1_mean_dwell c2_mean_dwell abs_log2_dwell_fc
            test_depth   (dirCLIP/test coverage, from cluster_counts; for depth matching)

Coordinate convention:
  Nanocompore `pos` is the 0-based transcript position of the first base of the 5-mer.
  We emit the full k-mer footprint start=pos, end=pos+5 (0-based half-open). The k-mer
  centre (pos+2) is derived downstream for the centre-only variant.

NOTE: verify the exact column-name expectations of the installed R2Dtool with
`r2d liftover --help` on the cluster. R2Dtool preserves trailing columns; the first
three (transcript,start,end) are the binding ones.
"""
import sys, math, argparse, csv

KMER = 5

def parse_test_depth(cluster_counts):
    """Sum dirCLIP (test = *_CLIP_*) read counts from the cluster_counts string.
    Format: 'sample:a/b__sample:a/b__...' where a/b are the two-cluster counts."""
    if cluster_counts in ("NA", "", None):
        return 0
    tot = 0
    for tok in cluster_counts.split("__"):
        if "_CLIP_" not in tok:
            continue
        try:
            ab = tok.split(":", 1)[1]
            a, b = ab.split("/")
            tot += int(a) + int(b)
        except (IndexError, ValueError):
            pass
    return tot

def fnum(x):
    return None if x in ("NA", "", None) else float(x)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--shift", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # index shift_stats by (ref_id, pos)
    shift = {}
    with open(a.shift) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            shift[(row["ref_id"], row["pos"])] = row

    n_out = n_miss = 0
    with open(a.results) as fh, open(a.out, "w") as out:
        r = csv.DictReader(fh, delimiter="\t")
        header = ["transcript", "start", "end", "strand",
                  "ref_kmer", "GMM_logit_pvalue", "Logit_LOR",
                  "c1_mean_intensity", "c2_mean_intensity", "abs_d_intensity",
                  "c1_mean_dwell", "c2_mean_dwell", "abs_log2_dwell_fc",
                  "test_depth"]
        out.write("\t".join(header) + "\n")
        for row in r:
            tid, pos = row["ref_id"], row["pos"]
            s = shift.get((tid, pos))
            if s is None:
                n_miss += 1
                continue
            start = int(pos)
            end = start + KMER
            c1i, c2i = fnum(s["c1_mean_intensity"]), fnum(s["c2_mean_intensity"])
            c1d, c2d = fnum(s["c1_mean_dwell"]), fnum(s["c2_mean_dwell"])
            d_int = abs(c2i - c1i) if (c1i is not None and c2i is not None) else "NA"
            if c1d and c2d and c1d > 0 and c2d > 0:
                dwell_fc = abs(math.log2(c2d / c1d))
            else:
                dwell_fc = "NA"
            fields = [tid, str(start), str(end), "+",
                      row["ref_kmer"], row["GMM_logit_pvalue"], row["Logit_LOR"],
                      s["c1_mean_intensity"], s["c2_mean_intensity"],
                      f"{d_int}" if d_int == "NA" else f"{d_int:.6f}",
                      s["c1_mean_dwell"], s["c2_mean_dwell"],
                      f"{dwell_fc}" if dwell_fc == "NA" else f"{dwell_fc:.6f}",
                      str(parse_test_depth(row["cluster_counts"]))]
            out.write("\t".join(fields) + "\n")
            n_out += 1

    sys.stderr.write(f"[prep] wrote {n_out} sites; {n_miss} results rows had no shift_stats match\n")

if __name__ == "__main__":
    main()
