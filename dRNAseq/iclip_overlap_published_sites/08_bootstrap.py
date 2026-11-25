#!/usr/bin/env python3
"""Step 8. Intervals for the crude and the Mantel Haenszel odds ratios by resampling transcripts and
by resampling loci.

A locus is a set of transcripts connected by overlap, on the same chromosome and strand, of the
genomic spans of their tested positions. Run once.

Usage: python3 08_bootstrap.py <r8_positions.tsv.gz> <a4a_lookup.tsv.gz> <r8_r32_table.csv of step 3> <out.csv>
"""
import sys
import numpy as np, pandas as pd

POS, LOOK, REF, OUT = sys.argv[1:5]
B, SEED, THR = 10000, 0, 0.005
TIERS = (('all', 'n_all'), ('hc15', 'n_hc15'), ('hc79', 'n_hc79'))

d = pd.read_csv(POS, sep='\t', dtype={'transcript': str, 'chr': str})
lk = pd.read_csv(LOOK, sep='\t', usecols=['transcript', 'tx_pos', 'our_GMM_p', 'our_LOR', 'test_depth'], low_memory=False)
lk['our_GMM_p'] = pd.to_numeric(lk.our_GMM_p, errors='coerce'); lk['our_LOR'] = pd.to_numeric(lk.our_LOR, errors='coerce')
ours = lk[lk.our_GMM_p.notna() & lk.our_LOR.notna()]
assert len(lk) == 344464 and len(ours) == 226360
edges = np.unique(np.quantile(ours.test_depth.values, np.linspace(0, 1, 11))).astype(float); edges[-1] += 1   # as in step 3
depth = dict(zip(zip(lk.transcript, lk.tx_pos), lk.test_depth))
d['depth'] = [depth.get((t, p - 3), np.nan) for t, p in zip(d.transcript, d.pos)]
d.loc[d['depth'] <= 0, 'depth'] = np.nan      # run 2 (addendum 2): 0 is the lookup's placeholder for no depth value
d['sig'] = (d.GMM_adj_pval < THR) & (d.LOR.abs() >= 1)
ref = pd.read_csv(REF)
rng = np.random.default_rng(SEED)
rows = []


def loci_of(da):
    """transcript -> locus id: connected components of overlapping genomic spans on one chromosome and strand"""
    sp = da.groupby('transcript').agg(chr=('chr', 'first'), strand=('strand', 'first'), lo=('gpos0', 'min'), hi=('gpos0', 'max'),
                                      nchr=('chr', 'nunique'), nstr=('strand', 'nunique')).reset_index()
    assert (sp.nchr == 1).all() and (sp.nstr == 1).all()
    sp = sp.sort_values(['chr', 'strand', 'lo', 'hi']).reset_index(drop=True)
    loc, cur, last = [], -1, None
    for r in sp.itertuples():
        if last is None or (r.chr, r.strand) != last[:2] or r.lo > last[2]:
            cur += 1; last = (r.chr, r.strand, r.hi)
        else:
            last = (r.chr, r.strand, max(last[2], r.hi))
        loc.append(cur)
    return dict(zip(sp.transcript, loc))


def tables(df, flag, strata, unit):
    """counts[u, s, :] = (a, b, c, d) per resampling unit and stratum"""
    names, ucode = np.unique(unit, return_inverse=True)
    cell = np.where(df.sig.values, np.where(df[flag].values, 0, 1), np.where(df[flag].values, 2, 3))
    out = np.zeros((len(names), strata.max() + 1, 4))
    np.add.at(out, (ucode, strata, cell), 1)
    return out


def mh_from(tab):
    a, b, c, dd = tab[..., 0], tab[..., 1], tab[..., 2], tab[..., 3]
    n = a + b + c + dd
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.where(n > 0, a * dd / n, 0.0).sum(axis=-1); q = np.where(n > 0, b * c / n, 0.0).sum(axis=-1)
    return r / q


def boot(tab):
    U = tab.shape[0]
    out = np.empty(B)
    for i in range(0, B, 1000):                                             # in blocks, to bound memory
        w = rng.multinomial(U, np.full(U, 1.0 / U), size=1000).astype(float)
        out[i:i + 1000] = mh_from(np.tensordot(w, tab, axes=(1, 0)))
    return out


for arm in ('+UV', '-UV'):
    da = d[d.arm == arm].copy()
    locus = loci_of(da)
    da['locus'] = da.transcript.map(locus)
    n_tx, n_loc = da.transcript.nunique(), da.locus.nunique()
    multi = da.groupby('locus').transcript.nunique()
    print(f'===== {arm}: transcripts {n_tx}, loci {n_loc}; loci with more than one transcript {int((multi > 1).sum())}, largest {int(multi.max())} transcripts')
    assert n_loc < n_tx and n_loc > 300
    for tier, col in TIERS:
        da['flag'] = da[col] > 0
        sub = da[da.depth.notna()]
        s_idx = np.clip(np.searchsorted(edges, sub.depth.values, side='right') - 1, 0, len(edges) - 2)
        zero = np.zeros(len(da), dtype=int)
        r = ref[(ref.comparison == arm) & (ref.peak_tier == tier) & (ref.adj_p_threshold == THR)].iloc[0]
        res = {}
        for unit in ('transcript', 'locus'):
            t_mh = tables(sub, 'flag', s_idx, sub[unit].values); t_cr = tables(da, 'flag', zero, da[unit].values)
            point_mh, point_cr = float(mh_from(t_mh.sum(axis=0))), float(mh_from(t_cr.sum(axis=0)))
            assert abs(point_mh - r.MH_OR) < 1e-5 and abs(point_cr - r.crude_OR) < 1e-5, (arm, tier, unit, point_mh, r.MH_OR, point_cr, r.crude_OR)
            bm, bc = boot(t_mh), boot(t_cr)
            res[unit] = (np.percentile(bc, [2.5, 97.5]), np.percentile(bm, [2.5, 97.5]), float((bm > 1).mean()), t_cr.shape[0], t_mh.shape[0])
        (ct, mt, _, _, _), (cl, ml, share, ul, uld) = res['transcript'], res['locus']
        rows.append(dict(comparison=arm, peak_tier=tier, loci_all=ul, loci_depth=uld,
                         crude_OR=r.crude_OR, crude_indep_lo=r.crude_lo, crude_indep_hi=r.crude_hi,
                         crude_tx_lo=round(float(ct[0]), 6), crude_tx_hi=round(float(ct[1]), 6),
                         crude_locus_lo=round(float(cl[0]), 6), crude_locus_hi=round(float(cl[1]), 6),
                         MH_OR=r.MH_OR, MH_indep_lo=r.MH_lo, MH_indep_hi=r.MH_hi,
                         MH_tx_lo=round(float(mt[0]), 6), MH_tx_hi=round(float(mt[1]), 6),
                         MH_locus_lo=round(float(ml[0]), 6), MH_locus_hi=round(float(ml[1]), 6),
                         MH_locus_share_above_1=round(share, 4)))
        print(f'{arm} {tier:5} crude {r.crude_OR:.3f}: independent {r.crude_lo:.3f} to {r.crude_hi:.3f}; transcripts {ct[0]:.3f} to {ct[1]:.3f}; loci {cl[0]:.3f} to {cl[1]:.3f}')
        print(f'          MH    {r.MH_OR:.3f}: independent {r.MH_lo:.3f} to {r.MH_hi:.3f}; transcripts {mt[0]:.3f} to {mt[1]:.3f}; loci {ml[0]:.3f} to {ml[1]:.3f}; share of locus resamples above 1: {share:.3f}')
pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'rows written {len(rows)} (expected 6); resamples {B}; seed {SEED}')
assert len(rows) == 6
