#!/usr/bin/env python3
"""Step 3. Crude and depth adjusted iCLIP overlap of the published SRSF3 site sets,
on the 11 nt genomic window of the coauthor's R script. Design: the preregistered design of 5 Oct 2026 (project record), sections 1 and 2.

Usage: python3 03_r32_stats.py <r8_positions.tsv.gz> <a4a_lookup.tsv.gz> <out_prefix>
"""
import math, sys
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, norm

POS, LOOK, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
TIERS = (('all', 'n_all'), ('hc15', 'n_hc15'), ('hc79', 'n_hc79'))
EXPECT_FINITE = {'+UV': 273052, '-UV': 248662}          # rows written by 01_prep_r2dtool_input.py
MOHAMED = {  # sigsites_overlap_with_iclip.xlsx, 29 Sep 2026: (sig %, nonsig %, Fisher p)
    ('-UV', 0.05): (34.95, 30.33, 1.744274e-13), ('-UV', 0.005): (36.3, 30.38, 8.830071e-10),
    ('+UV', 0.05): (29.34, 29.78, 0.7174285), ('+UV', 0.005): (27.64, 29.78, 0.3237127)}

d = pd.read_csv(POS, sep='\t', dtype={'transcript': str, 'chr': str})
d['absLOR'] = d.LOR.abs()
print(f'positions with a genomic window: {len(d):,} ' + str(d.arm.value_counts().to_dict()))
for arm, n in EXPECT_FINITE.items():
    k = int((d.arm == arm).sum())
    print(f'  {arm}: {k:,} of {n:,} tested positions lifted ({k / n * 100:.2f}%), {n - k} not lifted')
    assert 0.99 * n < k <= n

lk = pd.read_csv(LOOK, sep='\t', usecols=['transcript', 'tx_pos', 'our_GMM_p', 'our_LOR', 'test_depth'], low_memory=False)
lk['our_GMM_p'] = pd.to_numeric(lk.our_GMM_p, errors='coerce'); lk['our_LOR'] = pd.to_numeric(lk.our_LOR, errors='coerce')
ours = lk[lk.our_GMM_p.notna() & lk.our_LOR.notna()]
assert len(lk) == 344464 and len(ours) == 226360, (len(lk), len(ours))
depth = dict(zip(zip(lk.transcript, lk.tx_pos), lk.test_depth))
d['depth'] = [depth.get((t, p - 3), np.nan) for t, p in zip(d.transcript, d.pos)]
d.loc[d['depth'] <= 0, 'depth'] = np.nan      # run 2 (addendum 2): 0 is the lookup's placeholder for no depth value


def edges_for(n_strata):
    """Construction of 27 Sep 2026: quantile edges of test_depth in the 2026 rerun universe."""
    e = np.unique(np.quantile(ours.test_depth.values, np.linspace(0, 1, n_strata + 1)))
    e = e.astype(float); e[-1] += 1
    return e


def mh(df, flag, e):
    s = np.clip(np.searchsorted(e, df.depth.values, side='right') - 1, 0, len(e) - 2)
    R = S = A = E = V = 0.0
    P_R = P_S = PQ = Q_S = 0.0   # Robins Breslow Greenland terms
    used = 0
    for st in np.unique(s):
        g = df[s == st]
        a = int((g.sig & g[flag]).sum()); b = int((g.sig & ~g[flag]).sum())
        c = int((~g.sig & g[flag]).sum()); dd = int((~g.sig & ~g[flag]).sum())
        n = a + b + c + dd
        if n < 2 or (a + b) == 0 or (c + dd) == 0:
            continue
        used += 1
        r, q = a * dd / n, b * c / n
        R += r; S += q
        P, Q = (a + dd) / n, (b + c) / n
        P_R += P * r; P_S += P * q; PQ += Q * r; Q_S += Q * q
        m1, n1 = a + c, a + b
        A += a; E += n1 * m1 / n; V += n1 * (n - n1) * m1 * (n - m1) / (n * n * (n - 1))
    orr = R / S
    var = P_R / (2 * R * R) + (P_S + PQ) / (2 * R * S) + Q_S / (2 * S * S)
    lo, hi = math.exp(math.log(orr) - 1.96 * math.sqrt(var)), math.exp(math.log(orr) + 1.96 * math.sqrt(var))
    z = (A - E) / math.sqrt(V)
    return orr, lo, hi, float(norm.sf(z)), A, E, used


def crude(df, flag):
    a = int((df.sig & df[flag]).sum()); b = int((df.sig & ~df[flag]).sum())
    c = int((~df.sig & df[flag]).sum()); dd = int((~df.sig & ~df[flag]).sum())
    orr, p = fisher_exact([[a, b], [c, dd]])
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / dd)
    return a, b, c, dd, orr, math.exp(math.log(orr) - 1.96 * se), math.exp(math.log(orr) + 1.96 * se), p


rows = []
e10, e5, e20 = edges_for(10), edges_for(5), edges_for(20)
print(f'depth strata edges (10): {e10.tolist()}')
for arm in ('+UV', '-UV'):
    da = d[d.arm == arm].copy()
    for _, col in TIERS:
        da[col + '_f'] = da[col] > 0
    has = da.depth.notna()
    print(f'\n===== {arm}: {len(da):,} positions; with a depth value {int(has.sum()):,} ({has.mean() * 100:.1f}%)')
    for thr in (0.005, 0.05):
        da['sig'] = (da.GMM_adj_pval < thr) & (da.absLOR >= 1)
        ns, nsd = int(da.sig.sum()), int((da.sig & has).sum())
        print(f'\n--- {arm}, adjusted p < {thr}, |LOR| >= 1: significant {ns:,} (with depth {nsd:,}); '
              f'median depth sig {da.loc[da.sig & has, "depth"].median():.0f} vs other {da.loc[~da.sig & has, "depth"].median():.0f}')
        for tier, col in TIERS:
            f = col + '_f'
            a, b, c, dd, orr, lo, hi, p = crude(da, f)
            sub = da[has]
            a2, b2, c2, d2, or2, lo2, hi2, p2 = crude(sub, f)
            m10 = mh(sub, f, e10); m5 = mh(sub, f, e5); m20 = mh(sub, f, e20)
            exp_pct = m10[5] / (a2 + b2) * 100
            print(f'  {tier:5} ALL POSITIONS  sig {a}/{a + b} = {a / (a + b) * 100:.2f}%  other {c}/{c + dd} = {c / (c + dd) * 100:.2f}%  OR {orr:.3f} ({lo:.3f} to {hi:.3f})  Fisher two sided p {p:.3g}')
            print(f'        WITH DEPTH     sig {a2}/{a2 + b2} = {a2 / (a2 + b2) * 100:.2f}%  other {c2}/{c2 + d2} = {c2 / (c2 + d2) * 100:.2f}%  OR {or2:.3f} ({lo2:.3f} to {hi2:.3f})  Fisher p {p2:.3g}')
            print(f'        DEPTH ADJUSTED MH OR {m10[0]:.3f} ({m10[1]:.3f} to {m10[2]:.3f}), one sided CMH p {m10[3]:.3g}, observed {m10[4]:.0f} expected {m10[5]:.1f} ({exp_pct:.2f}%), strata used {m10[6]};  5 strata {m5[0]:.3f}, 20 strata {m20[0]:.3f}')
            rows.append(dict(comparison=arm, adj_p_threshold=thr, peak_tier=tier,
                             n_sig=a + b, sig_in_peak=a, sig_pct=round(a / (a + b) * 100, 2),
                             n_other=c + dd, other_in_peak=c, other_pct=round(c / (c + dd) * 100, 2),
                             crude_OR=round(orr, 6), crude_lo=round(lo, 6), crude_hi=round(hi, 6), fisher_two_sided_p=p,
                             n_sig_depth=a2 + b2, sig_in_peak_depth=a2, sig_pct_depth=round(a2 / (a2 + b2) * 100, 2),
                             other_pct_depth=round(c2 / (c2 + d2) * 100, 2), crude_OR_depth=round(or2, 6),
                             depth_matched_expected_pct=round(exp_pct, 2),
                             MH_OR=round(m10[0], 6), MH_lo=round(m10[1], 6), MH_hi=round(m10[2], 6), CMH_one_sided_p=m10[3],
                             MH_OR_5=round(m5[0], 6), MH_OR_20=round(m20[0], 6)))
            if tier == 'all':
                ms, mn, mp = MOHAMED[(arm, thr)]
                print(f'        coauthor table, 29 Sep: sig {ms}%  other {mn}%  Fisher p {mp:.3g}   | difference here minus his: sig {a / (a + b) * 100 - ms:+.2f}, other {c / (c + dd) * 100 - mn:+.2f} points')
pd.DataFrame(rows).to_csv(OUT + '_r32_table.csv', index=False)
print(f'\nrows written: {len(rows)} (expected 12) -> {OUT}_r32_table.csv')
assert len(rows) == 12
