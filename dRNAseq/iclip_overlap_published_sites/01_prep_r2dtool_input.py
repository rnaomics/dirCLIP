#!/usr/bin/env python3
"""R8 step 1 (local). Build the R2Dtool liftover input for every tested position of the two
published Nanocompore sheets of Supplementary Table 1.

Usage: python3 01_prep_r2dtool_input.py <SupT1_tsv_dir> <out.tsv>
<SupT1_tsv_dir> is the output of ../supplementary_table1_checks/00_extract_sup_table1.py.

Coordinate convention. `pos` in Supplementary Table 1 is the 1 based transcript coordinate of the
kmer centre (Nanocompore pos + 3; pipeline script 7). R2Dtool takes a 0 based start and lifts that
single nucleotide, so start = pos - 1 and end = pos.
Universe: rows with a finite LOR and a finite GMM_adj_pval (the coauthor's R script drops LOR == "NC").
"""
import math, sys

src, out = sys.argv[1], sys.argv[2]


def num(s):
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


n_total = n_kept = 0
with open(out, 'w') as fo:
    fo.write('transcript\tstart\tend\tstrand\tarm\tpos\tLOR\tGMM_adj_pval\n')
    for arm in ('+UV', '-UV'):
        n_arm = k_arm = lor_only = p_only = 0
        with open(f'{src}/{arm}_SRSF3_Nanocompore.tsv') as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            assert hdr == ['pos', 'ref_id', 'LOR', 'GMM_adj_pval'], hdr
            for line in fh:
                f = line.rstrip('\n').split('\t')
                n_arm += 1
                lor, p = num(f[2]), num(f[3])
                if lor is None or p is None:
                    lor_only += lor is not None
                    p_only += p is not None
                    continue
                pos = int(float(f[0]))
                assert pos >= 1
                fo.write(f'{f[1]}\t{pos - 1}\t{pos}\t+\t{arm}\t{pos}\t{f[2]}\t{f[3]}\n')
                k_arm += 1
        print(f'{arm}: rows {n_arm:,}; finite LOR and p {k_arm:,}; finite LOR only {lor_only}; finite p only {p_only}')
        n_total += n_arm
        n_kept += k_arm
print(f'written {n_kept:,} of {n_total:,} rows to {out}')
