#!/usr/bin/env python3
"""Step 9. Two checks that share no code with steps 3 and 5.
(1) Mantel Haenszel odds ratio, its Robins Breslow Greenland interval and the CMH statistic from
    statsmodels.stats.contingency_tables.StratifiedTable, against r8_r32_table.csv.
(2) Covered and recovered peaks by a pandas join, against r8_r34_table.csv.
Usage: python3 09_independent_checks.py <outputs_dir> <a4a_lookup.tsv.gz>"""
import sys
import numpy as np, pandas as pd
from scipy.stats import norm
from statsmodels.stats.contingency_tables import StratifiedTable

O, LOOK = sys.argv[1], sys.argv[2]
d = pd.read_csv(f'{O}/r8_positions.tsv.gz', sep='\t', dtype={'transcript': str, 'chr': str})
pr = pd.read_csv(f'{O}/r8_pairs.tsv.gz', sep='\t', dtype={'transcript': str, 'peak': str})
ref = pd.read_csv(f'{O}/r8_r32_table.csv'); r34 = pd.read_csv(f'{O}/r8_r34_table.csv')
lk = pd.read_csv(LOOK, sep='\t', usecols=['transcript', 'tx_pos', 'our_GMM_p', 'our_LOR', 'test_depth'], low_memory=False)
fin = pd.to_numeric(lk.our_GMM_p, errors='coerce').notna() & pd.to_numeric(lk.our_LOR, errors='coerce').notna()
edges = np.unique(np.quantile(lk.test_depth[fin].values, np.linspace(0, 1, 11))).astype(float); edges[-1] += 1
m = d.merge(lk[['transcript', 'tx_pos', 'test_depth']].assign(pos=lk.tx_pos + 3), how='left', on=['transcript', 'pos'])
assert len(m) == len(d)
m.loc[m.test_depth <= 0, 'test_depth'] = np.nan       # run 2: placeholder zeros are not depths
m['sig'] = (m.GMM_adj_pval < 0.005) & (m.LOR.abs() >= 1)
n_checked = 0
for arm in ('+UV', '-UV'):
    for tier, col in (('all', 'n_all'), ('hc15', 'n_hc15'), ('hc79', 'n_hc79')):
        x = m[(m.arm == arm) & m.test_depth.notna()]
        st = np.clip(np.searchsorted(edges, x.test_depth.values, side='right') - 1, 0, len(edges) - 2)
        tabs = []
        for s in np.unique(st):
            g = x[st == s]; hit = g[col] > 0
            tabs.append([[int((g.sig & hit).sum()), int((g.sig & ~hit).sum())], [int((~g.sig & hit).sum()), int((~g.sig & ~hit).sum())]])
        T = StratifiedTable(tabs)
        lo, hi = T.oddsratio_pooled_confint()
        z = np.sqrt(T.test_null_odds(correction=False).statistic) * np.sign(T.oddsratio_pooled - 1)
        r = ref[(ref.comparison == arm) & (ref.peak_tier == tier) & (ref.adj_p_threshold == 0.005)].iloc[0]
        ok = abs(T.oddsratio_pooled - r.MH_OR) < 1e-5 and abs(lo - r.MH_lo) < 1e-4 and abs(hi - r.MH_hi) < 1e-4 and abs(norm.sf(z) - r.CMH_one_sided_p) < 1e-6 + 1e-3 * r.CMH_one_sided_p
        print(f'{arm} {tier:5} statsmodels MH {T.oddsratio_pooled:.6f} ({lo:.6f} to {hi:.6f}), one sided p {norm.sf(z):.4g} | step 3: {r.MH_OR:.6f} ({r.MH_lo:.6f} to {r.MH_hi:.6f}), {r.CMH_one_sided_p:.4g} | {"agree" if ok else "DIFFER"}')
        assert ok; n_checked += 1
        # peaks by join
        cut = {'all': 0, 'hc15': 15, 'hc79': 79}[tier]
        pa = pr[(pr.arm == arm) & (pr.score >= cut)].merge(m[m.arm == arm][['transcript', 'pos', 'sig', 'test_depth']], on=['transcript', 'pos'], how='left')
        assert pa.sig.notna().all()
        cov, rec = pa.peak.nunique(), pa[pa.sig].peak.nunique()
        cov_d, rec_d = pa[pa.test_depth.notna()].peak.nunique(), pa[pa.sig & pa.test_depth.notna()].peak.nunique()
        q = r34[(r34.comparison == arm) & (r34.peak_tier == tier)].iloc[0]
        ok2 = (cov, rec, cov_d, rec_d) == (q.covered, q.recovered, q.covered_depth, q.recovered_depth)
        print(f'          peaks by join: covered {cov}, recovered {rec}; with depth {cov_d}, {rec_d} | step 5: {q.covered}, {q.recovered}; {q.covered_depth}, {q.recovered_depth} | {"agree" if ok2 else "DIFFER"}')
        assert ok2; n_checked += 1
print(f'checks passed: {n_checked} (expected 12)')
assert n_checked == 12
