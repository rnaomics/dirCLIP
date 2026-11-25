#!/usr/bin/env python3
"""R8 diagnostic (local). Positions on the two transcripts whose version differs between
Supplementary Table 1 and the Ensembl 114 GTF were lifted with the newer exon structure.
How many are they, and does dropping them move the crude 2x2? Also lists what was not lifted.
Usage: python3 06_version_sensitivity.py <r8_positions.tsv.gz> <r8_transcript_versions.tsv> <SupT1_tsv_dir>"""
import sys, math
import pandas as pd
from scipy.stats import fisher_exact
d = pd.read_csv(sys.argv[1], sep='\t', dtype={'transcript': str, 'chr': str})
v = pd.read_csv(sys.argv[2], sep='\t', header=None, names=['transcript', 'gtf_version', 'status'], dtype=str)
print(v.status.value_counts().to_dict(), '| transcripts examined', len(v))
assert len(v) > 900
differs = set(v.loc[v.status == 'differs', 'transcript']); absent = set(v.loc[v.status == 'absent', 'transcript'])
d['sig'] = (d.GMM_adj_pval < 0.005) & (d.LOR.abs() >= 1); d['hit'] = d.n_all > 0
for arm in ('+UV', '-UV'):
    da = d[d.arm == arm]
    k = pd.read_csv(f'{sys.argv[3]}/{arm}_SRSF3_Nanocompore.tsv', sep='\t', dtype=str)
    k['p'] = pd.to_numeric(k.GMM_adj_pval, errors='coerce'); k['lor'] = pd.to_numeric(k.LOR, errors='coerce')
    k = k[k.p.notna() & k.lor.notna()]
    ka = k[k.ref_id.isin(absent)]
    print(f'{arm}: not lifted (transcript absent from the GTF): {len(ka):,} positions on {ka.ref_id.nunique()} transcripts, '
          f'of which significant {int(((ka.p < 0.005) & (ka.lor.abs() >= 1)).sum())}')
    dm = da[da.transcript.isin(differs)]
    print(f'   on version differing transcripts: {len(dm):,} positions, significant {int(dm.sig.sum())}, in peak {int(dm.hit.sum())}')
    for label, x in (('all lifted', da), ('without version differing transcripts', da[~da.transcript.isin(differs)])):
        a = int((x.sig & x.hit).sum()); b = int((x.sig & ~x.hit).sum()); c = int((~x.sig & x.hit).sum()); dd = int((~x.sig & ~x.hit).sum())
        orr, p = fisher_exact([[a, b], [c, dd]])
        print(f'   {label}: sig {a}/{a+b} = {a/(a+b)*100:.2f}%  other {c}/{c+dd} = {c/(c+dd)*100:.2f}%  OR {orr:.3f}  p {p:.3g}')
