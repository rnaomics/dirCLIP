# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
"""Independent check, written from scratch, sharing no code with A2 script 05.
Question: is Nanocompore's GMM_logit_pvalue column raw or BH adjusted?
Also: is the 'BH fixed point' test used on 10 Sep 2026 able to tell them apart?"""
import gzip, math
import numpy as np

import sys
F = sys.argv[1]      # a Nanocompore 1.0.x results table, .tsv or .tsv.gz

def my_bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p, kind='mergesort')
    ranked = p[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(ranked, 1.0)
    return out

def norm_sf2(z):  # two sided normal tail
    return math.erfc(abs(z) / math.sqrt(2))

rep, raw, lor_rep, lor_mine = [], [], [], []
bad = 0
with (gzip.open(F, 'rt') if F.endswith('.gz') else open(F)) as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    ip, ic, il = hdr.index('GMM_logit_pvalue'), hdr.index('cluster_counts'), hdr.index('Logit_LOR')
    for line in fh:
        f = line.rstrip('\n').split('\t')
        try:
            q = float(f[ip]); L = float(f[il])
        except ValueError:
            continue
        if not (math.isfinite(q) and math.isfinite(L)):
            continue
        a = b = c = d = 1.0
        ok = True
        for part in f[ic].split('__'):
            name, cnt = part.rsplit(':', 1)
            k0, k1 = (int(x) for x in cnt.split('/'))
            if 'noUV' in name: a += k0; b += k1
            elif 'CLIP' in name: c += k0; d += k1
            else: ok = False
        if not ok:
            bad += 1; continue
        l = math.log(a * d / (b * c))
        se = math.sqrt(1/a + 1/b + 1/c + 1/d)
        rep.append(q); raw.append(norm_sf2(l / se)); lor_rep.append(L); lor_mine.append(l)

rep = np.array(rep); raw = np.array(raw)
print(f'rows used {len(rep):,}   rows skipped for unexpected sample names {bad}')
print(f'max |my LOR - reported LOR| = {np.max(np.abs(np.array(lor_mine) - np.array(lor_rep))):.2e}')
eq_raw = np.mean(np.isclose(rep, raw, rtol=1e-6, atol=0))
print(f'reported == closed form raw p (rtol 1e-6): {eq_raw*100:.2f} %')
bh = my_bh(raw)
eq_bh = np.mean(np.isclose(rep, bh, rtol=1e-6, atol=0))
print(f'reported == BH(closed form raw p) (rtol 1e-6): {eq_bh*100:.2f} %')
print(f'fraction reported > raw: {np.mean(rep > raw*(1+1e-9))*100:.2f} %')
# The 10 Sep 'fixed point' test, applied to a vector KNOWN to be BH adjusted
q = bh
fp = np.mean(np.isclose(my_bh(q), q, rtol=1e-9, atol=0))
print(f'fixed point test on a vector KNOWN to be BH adjusted: BH(q)==q in {fp*100:.2f} % of rows')
print(f'  fraction of that vector exactly at its maximum: {np.mean(q == q.max())*100:.2f} %')
fr = np.mean(np.isclose(my_bh(raw), raw, rtol=1e-9, atol=0))
print(f'fixed point test on the raw vector: {fr*100:.2f} %')
print('site rule counts on reported column, |LOR|>=1:')
L = np.abs(np.array(lor_rep))
for thr in (0.005, 0.05):
    print(f'  reported < {thr}: {int(np.sum((rep < thr) & (L >= 1))):,}    raw < {thr}: {int(np.sum((raw < thr) & (L >= 1))):,}')
