#!/usr/bin/env python3
"""
Extract significant Nanocompore signal deviation sites and collapse them into regions,
using the rule that actually produced Supplementary Table 1.

WHY THIS FILE EXISTS
--------------------
An earlier extractor (not included here) was written to be faithful to the archived R script 7.
It was faithful, and the archived R is wrong. Established 16 Aug 2026 on the 2026 rerun.

The archived R collapses within 5 nt, in script 7 and again in script 10.
Applied to the published position sheet (Supplementary Table 1) that yields 1,712 regions with signal_dev
recall 0.848 and precision 0.678. It does not reproduce the published table.

A single pass collapse within 10 nt with a floored centre reproduces it exactly:
  -UV   1,395 regions produced against 1,395 reported carrying signal evidence,
        signal_dev recall 1,369 of 1,369 = 1.0000
  +UV   363 produced against 363 reported, signal_dev recall 359 of 359 = 1.0000,
        held out, and an independent refit of the gap on +UV alone also lands on 10
  negative control, read_ends only regions reproduced from the Nanocompore sheet:
        0 of 498 on -UV, 0 of 379 on +UV

Effect of the correction on our own runs:
  control rep1 and rep2   704 regions at gap 5  ->  584 at gap 10 with floor
  control rep2 only       799 regions at gap 5  ->  591 at gap 10 with floor

THE RULE
--------
1. recentre position to the 5mer centre, pos + 3
2. keep GMM_logit_pvalue < P_THRESH (default 0.005; the column is Benjamini Hochberg adjusted
   by Nanocompore 1.0.x) AND abs(Logit_LOR) >= 1
3. merge surviving positions on a transcript when they are within 10 nt
4. region centre = floor((min(pos) + max(pos)) / 2)
5. report the region as centre minus 5 to centre plus 5, a fixed 10 nt window

Steps 4 and 5 are what Supplementary Table 1 reports. Every region in both halves of that
workbook is exactly 10 nt wide with centre equal to start plus 5, 1,893 of 1,893 and
742 of 742. The merged span itself is not preserved in the published table, so it is carried here in
separate columns rather than being lost.

Usage:  extract_sig_sites_v2.py <nanocompore_results.tsv> <outdir> [prefix] [p_threshold]
Outputs in <outdir>:
    <prefix>_v2_significant_sites_raw.tsv       per position, pos = kmer centre
    <prefix>_v2_sig_sites_collapsed_regions.tsv regions, gap 10, floored centre
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

P_THRESH = float(sys.argv[4]) if len(sys.argv) > 4 else 0.005   # optional 4th argument (10 Oct 2026); 0.05 for HMGA1a
LOR_THRESH = 1.0
MERGE_GAP = 10
HALF_WINDOW = 5


def main():
    res_path = Path(sys.argv[1])
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = sys.argv[3] if len(sys.argv) > 3 else res_path.stem

    df = pd.read_csv(res_path, sep="\t",
                     usecols=["pos", "ref_id", "ref_kmer",
                              "GMM_logit_pvalue", "Logit_LOR"])
    df["GMM_logit_pvalue"] = pd.to_numeric(df["GMM_logit_pvalue"], errors="coerce")
    df["Logit_LOR"] = pd.to_numeric(df["Logit_LOR"], errors="coerce")
    df["pos"] = df["pos"] + 3

    sig = df[(df.GMM_logit_pvalue < P_THRESH) &
             (df.Logit_LOR.abs() >= LOR_THRESH)].copy()
    sig.to_csv(outdir / f"{prefix}_v2_significant_sites_raw.tsv", sep="\t", index=False)

    sig = sig.sort_values(["ref_id", "pos"])
    regions = []
    for ref, g in sig.groupby("ref_id"):
        pos = g["pos"].to_numpy()
        lor = g["Logit_LOR"].abs().to_numpy()
        pval = g["GMM_logit_pvalue"].to_numpy()
        grp = np.concatenate([[0], (np.diff(pos) > MERGE_GAP).cumsum()])
        for gi in np.unique(grp):
            m = grp == gi
            lo, hi = int(pos[m].min()), int(pos[m].max())
            centre = int(np.floor((lo + hi) / 2))
            regions.append((ref, centre - HALF_WINDOW, centre + HALF_WINDOW, centre,
                            lo, hi, int(m.sum()),
                            float(lor[m].max()), float(pval[m].min())))
    reg = pd.DataFrame(regions, columns=["ref_id", "region_start", "region_end",
                                         "region_center", "merged_span_start",
                                         "merged_span_end", "site_count",
                                         "max_abs_lor", "min_raw_pval"])
    reg.to_csv(outdir / f"{prefix}_v2_sig_sites_collapsed_regions.tsv", sep="\t", index=False)
    print(f"[{prefix}] gap={MERGE_GAP} floored centre  "
          f"significant positions={len(sig)}  collapsed regions={len(reg)} "
          f"on {reg.ref_id.nunique()} transcripts")


if __name__ == "__main__":
    main()
