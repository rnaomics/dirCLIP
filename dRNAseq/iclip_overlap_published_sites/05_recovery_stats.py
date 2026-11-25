#!/usr/bin/env python3
"""Step 5. Peak level recovery of SRSF3 iCLIP peaks by the published site sets,
on the 11 nt genomic windows. Design: the preregistered design of 5 Oct 2026 (project record), section 3.

Usage: python3 05_recovery_stats.py <r8_positions.tsv.gz> <r8_pairs.tsv.gz> <a4a_lookup.tsv.gz> <out_prefix>
"""
import sys
import numpy as np, pandas as pd

POS, PAIRS, LOOK, OUT = sys.argv[1:5]
N_PEAKS = {'all': 370675, 'hc15': 38724, 'hc79': 3754}   # lifted peak sets, line counts on Kaya
CUT = {'all': 0, 'hc15': 15, 'hc79': 79}
KAYA_COVERED_ALL = {'+UV': 4312, '-UV': 4307}             # bedtools intersect -u, 04_offbyone_check_OUTPUT.txt
REPS, SEED, N_STRATA = 1000, 0, 10

d = pd.read_csv(POS, sep='\t', dtype={'transcript': str, 'chr': str})
pr = pd.read_csv(PAIRS, sep='\t', dtype={'transcript': str, 'peak': str})
lk = pd.read_csv(LOOK, sep='\t', usecols=['transcript', 'tx_pos', 'test_depth'], low_memory=False)
depth = dict(zip(zip(lk.transcript, lk.tx_pos), lk.test_depth))
d['depth'] = [depth.get((t, p - 3), np.nan) for t, p in zip(d.transcript, d.pos)]
d.loc[d['depth'] <= 0, 'depth'] = np.nan      # run 2 (addendum 2): 0 is the lookup's placeholder for no depth value
d['sig'] = (d.GMM_adj_pval < 0.005) & (d.LOR.abs() >= 1)
print(f'positions {len(d):,}; position to peak pairs {len(pr):,}; distinct peaks in pairs {pr.peak.nunique():,}')
assert len(pr) > 100000 and pr.score.min() >= 1

rows = []
rng = np.random.default_rng(SEED)
for arm in ('+UV', '-UV'):
    da = d[d.arm == arm].reset_index(drop=True)
    key = {k: i for i, k in enumerate(zip(da.transcript, da.pos))}
    pa = pr[pr.arm == arm]
    pidx = np.array([key[k] for k in zip(pa.transcript, pa.pos)])
    peak_codes, peak_names = pd.factorize(pa.peak)
    score = pa.score.values
    # gate: positions flagged in the count table are exactly the positions present in the pairs table
    assert set(np.unique(pidx)) == set(np.flatnonzero(da.n_all.values > 0)), 'pairs and counts disagree'
    order = np.argsort(pidx, kind='stable')
    pidx_s, codes_s, score_s = pidx[order], peak_codes[order], score[order]
    starts = np.searchsorted(pidx_s, np.arange(len(da)), side='left')
    ends = np.searchsorted(pidx_s, np.arange(len(da)), side='right')

    def peaks_of(pos_idx, cut):
        """distinct peak codes (score >= cut) hit by the windows of the given positions, and number of positions hitting any"""
        hit, npos = [], 0
        for i in pos_idx:
            s, e = starts[i], ends[i]
            if e > s:
                c = codes_s[s:e][score_s[s:e] >= cut]
                if len(c):
                    hit.append(c); npos += 1
        return (np.unique(np.concatenate(hit)) if hit else np.array([], dtype=int)), npos

    has = da.depth.notna().values
    sig = da.sig.values
    pool = np.flatnonzero(has)                       # D2: the whole eligible pool, significant positions included
    edges = np.unique(np.quantile(da.depth.values[pool], np.linspace(0, 1, N_STRATA + 1)))
    stratum = np.digitize(da.depth.values, edges[1:-1], right=True)
    sig_d = np.flatnonzero(sig & has)
    need = np.bincount(stratum[sig_d], minlength=len(edges) - 1)
    pools = [pool[stratum[pool] == s] for s in range(len(edges) - 1)]
    assert all(len(pools[s]) >= need[s] for s in range(len(need))), 'a depth stratum is exhausted'
    print(f'\n===== {arm}: tested {len(da):,}; significant {int(sig.sum()):,}; with depth: tested {len(pool):,}, significant {len(sig_d):,}; '
          f'strata {len(edges) - 1}, smallest pool {min(len(p) for p in pools):,}, largest requirement {need.max()}')
    for tier in ('all', 'hc15', 'hc79'):
        cut = CUT[tier]
        cov_all, _ = peaks_of(np.arange(len(da)), cut)
        rec_all, npos_all = peaks_of(np.flatnonzero(sig), cut)
        cov_d, _ = peaks_of(pool, cut)
        rec_d, npos_d = peaks_of(sig_d, cut)
        if tier == 'all':
            assert len(cov_all) == KAYA_COVERED_ALL[arm], (len(cov_all), KAYA_COVERED_ALL[arm])
        null_rec, null_pos = np.empty(REPS), np.empty(REPS)
        for r in range(REPS):
            draw = np.concatenate([rng.choice(pools[s], size=need[s], replace=False) for s in range(len(need)) if need[s]])
            pk, npos = peaks_of(draw, cut)
            null_rec[r], null_pos[r] = len(pk), npos
        emp = (1 + int((null_rec >= len(rec_d)).sum())) / (1 + REPS)
        lo, hi = np.percentile(null_rec, [2.5, 97.5])
        print(f'  {tier:5} peaks {N_PEAKS[tier]:,}; covered {len(cov_all):,} ({len(cov_all) / N_PEAKS[tier] * 100:.2f}%); '
              f'recovered {len(rec_all):,} = {len(rec_all) / len(cov_all) * 100:.2f}% of covered')
        print(f'        with depth: covered {len(cov_d):,}; recovered {len(rec_d):,} = {len(rec_d) / len(cov_d) * 100:.2f}%; '
              f'depth matched null {null_rec.mean():.1f} = {null_rec.mean() / len(cov_d) * 100:.2f}% (95% {lo:.0f} to {hi:.0f}); one sided empirical p {emp:.3g}')
        print(f'        positions in a peak per recovered peak: observed {npos_d}/{len(rec_d)} = {npos_d / max(len(rec_d), 1):.3f}; null {null_pos.mean():.1f}/{null_rec.mean():.1f} = {null_pos.mean() / null_rec.mean():.3f}')
        rows.append(dict(comparison=arm, peak_tier=tier, peaks_total=N_PEAKS[tier], covered=len(cov_all),
                         covered_pct=round(len(cov_all) / N_PEAKS[tier] * 100, 2), recovered=len(rec_all),
                         recovered_pct_of_covered=round(len(rec_all) / len(cov_all) * 100, 2),
                         covered_depth=len(cov_d), recovered_depth=len(rec_d),
                         recovered_pct_depth=round(len(rec_d) / len(cov_d) * 100, 2),
                         null_mean=round(float(null_rec.mean()), 1), null_pct=round(float(null_rec.mean()) / len(cov_d) * 100, 2),
                         null_lo=float(lo), null_hi=float(hi), empirical_p_one_sided=emp,
                         sig_positions_in_peak=npos_d, positions_per_recovered_peak=round(npos_d / max(len(rec_d), 1), 3),
                         null_positions_per_peak=round(float(null_pos.mean() / null_rec.mean()), 3)))
pd.DataFrame(rows).to_csv(OUT + '_r34_table.csv', index=False)
print(f'\nrows written: {len(rows)} (expected 6) -> {OUT}_r34_table.csv')
assert len(rows) == 6
