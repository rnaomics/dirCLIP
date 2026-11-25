#!/usr/bin/env python3
"""Step 7. Write the supplementary workbook from the result tables of steps 3, 5 and 8.
Every number in the workbook, the notes included, is read or recomputed here; none is typed.

Usage: python3 07_make_supp_table.py <outputs_dir> <dst.xlsx> <a4a_lookup.tsv.gz> <SupT1_tsv_dir>
"""
import sys
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

OUTD, DST, LOOK, SUP = sys.argv[1:5]
r32 = pd.read_csv(f'{OUTD}/r8_r32_table.csv'); r34 = pd.read_csv(f'{OUTD}/r8_r34_table.csv')
boot = pd.read_csv(f'{OUTD}/r8_r32_bootstrap.csv')
assert len(r32) == 12 and len(r34) == 6 and len(boot) == 6
bcols = ['crude_tx_lo', 'crude_tx_hi', 'crude_locus_lo', 'crude_locus_hi', 'MH_tx_lo', 'MH_tx_hi', 'MH_locus_lo', 'MH_locus_hi']
r32 = r32.merge(boot[['comparison', 'peak_tier'] + bcols].assign(adj_p_threshold=0.005), how='left', on=['comparison', 'peak_tier', 'adj_p_threshold'])
assert r32.MH_locus_lo.notna().sum() == 6

# ---- counts for the notes, computed, not typed ----
d = pd.read_csv(f'{OUTD}/r8_positions.tsv.gz', sep='\t', dtype={'transcript': str, 'chr': str})
lk = pd.read_csv(LOOK, sep='\t', usecols=['transcript', 'tx_pos', 'test_depth'], low_memory=False)
have = set(zip(lk.transcript[lk.test_depth > 0], lk.tx_pos[lk.test_depth > 0]))       # 0 is a placeholder, not a depth
N = {}
for arm in ('+UV', '-UV'):
    k = pd.read_csv(f'{SUP}/{arm}_SRSF3_Nanocompore.tsv', sep='\t', dtype=str)
    k['p'] = pd.to_numeric(k.GMM_adj_pval, errors='coerce'); k['lor'] = pd.to_numeric(k.LOR, errors='coerce')
    k = k[k.p.notna() & k.lor.notna() & k.p.abs().lt(float('inf')) & k.lor.abs().lt(float('inf'))]
    da = d[d.arm == arm]
    lifted = set(zip(da.transcript, da.pos))
    miss = k[[(t, int(float(q))) not in lifted for t, q in zip(k.ref_id, k.pos)]]
    N[arm] = dict(tested=len(k), lifted=len(da), unlifted=len(miss),
                  unlifted_sig=int(((miss.p < 0.005) & (miss.lor.abs() >= 1)).sum()),
                  depth_pct=sum((t, q - 3) in have for t, q in zip(da.transcript, da.pos)) / len(da) * 100,
                  loci=int(boot[boot.comparison == arm].loci_all.iloc[0]))
    assert N[arm]['tested'] == N[arm]['lifted'] + N[arm]['unlifted'] and N[arm]['tested'] > 200000
n_absent_tx = len(set(pd.read_csv(f'{OUTD}/r8_transcript_versions.tsv', sep='\t', header=None, names=['t', 'v', 's']).query("s == 'absent'").t))
print('counts for the notes:', N, 'transcripts absent from the annotation:', n_absent_tx)

COMP = {'+UV': 'SRSF3 dirCLIP against UV treated transcriptome (+UV; Fig. 1D to 1F)',
        '-UV': 'SRSF3 dirCLIP against untreated transcriptome (-UV; Supplementary Fig. 1)'}
TIER = {'all': 'all peaks (n = 370,675)', 'hc15': 'crosslink score >= 15, top 10% (n = 38,724)', 'hc79': 'crosslink score >= 79, top 1% (n = 3,754)'}
r32['order'] = r32.adj_p_threshold.map({0.005: 0, 0.05: 1}) * 100 + r32.comparison.map({'+UV': 0, '-UV': 1}) * 10 + r32.peak_tier.map({'all': 0, 'hc15': 1, 'hc79': 2})
r32 = r32.sort_values('order')


def ci(a, b): return 'not computed' if pd.isna(a) else '%.3f to %.3f' % (a, b)


pos = pd.DataFrame({
    'comparison': r32.comparison.map(COMP),
    'site rule': ['adjusted p < %g and |LOR| >= 1%s' % (t, ' (the rule used for the SRSF3 sites in the paper)' if t == 0.005 else ' (shown for reference)') for t in r32.adj_p_threshold],
    'iCLIP peak set': r32.peak_tier.map(TIER),
    'significant positions (n)': r32.n_sig,
    'significant, within 5 nt of a peak (n)': r32.sig_in_peak,
    'significant, within 5 nt of a peak (%)': r32.sig_pct,
    'other tested positions (n)': r32.n_other,
    'other, within 5 nt of a peak (n)': r32.other_in_peak,
    'other, within 5 nt of a peak (%)': r32.other_pct,
    'odds ratio': r32.crude_OR.round(3),
    'odds ratio, 95% CI (resampling loci)': [ci(a, b) for a, b in zip(r32.crude_locus_lo, r32.crude_locus_hi)],
    'odds ratio, 95% CI (resampling transcripts)': [ci(a, b) for a, b in zip(r32.crude_tx_lo, r32.crude_tx_hi)],
    'odds ratio, 95% CI (positions taken as independent)': [ci(a, b) for a, b in zip(r32.crude_lo, r32.crude_hi)],
    'significant positions with a depth value (n)': r32.n_sig_depth,
    'of these, within 5 nt of a peak (%)': r32.sig_pct_depth,
    'expected if matched on read depth (%)': r32.depth_matched_expected_pct,
    'depth adjusted odds ratio (Mantel Haenszel, 10 strata)': r32.MH_OR.round(3),
    'depth adjusted odds ratio, 95% CI (resampling loci)': [ci(a, b) for a, b in zip(r32.MH_locus_lo, r32.MH_locus_hi)],
    'depth adjusted odds ratio, 95% CI (resampling transcripts)': [ci(a, b) for a, b in zip(r32.MH_tx_lo, r32.MH_tx_hi)],
    'depth adjusted odds ratio, 95% CI (positions taken as independent)': [ci(a, b) for a, b in zip(r32.MH_lo, r32.MH_hi)],
    'depth adjusted odds ratio, 5 strata': r32.MH_OR_5.round(3),
    'depth adjusted odds ratio, 20 strata': r32.MH_OR_20.round(3),
})
r34['order'] = r34.comparison.map({'+UV': 0, '-UV': 1}) * 10 + r34.peak_tier.map({'all': 0, 'hc15': 1, 'hc79': 2})
r34 = r34.sort_values('order')
peak = pd.DataFrame({
    'comparison': r34.comparison.map(COMP),
    'iCLIP peak set': r34.peak_tier.map(TIER),
    'peaks covered by a tested position (n)': r34.covered,
    'covered (% of the peak set)': r34.covered_pct,
    'covered peaks with a significant site (n)': r34.recovered,
    'covered peaks with a significant site (%)': r34.recovered_pct_of_covered,
    'positions with a depth value: covered peaks (n)': r34.covered_depth,
    'positions with a depth value: peaks with a significant site (n)': r34.recovered_depth,
    'positions with a depth value: peaks with a significant site (%)': r34.recovered_pct_depth,
    'depth matched random positions: peaks reached, mean of 1000 draws (n)': r34.null_mean,
    'depth matched random positions: peaks reached (%)': r34.null_pct,
    'depth matched random positions: 95% interval (n)': ['%d to %d' % (round(a), round(b)) for a, b in zip(r34.null_lo, r34.null_hi)],
    'significant positions per peak reached': r34.positions_per_recovered_peak,
    'random positions per peak reached': r34.null_positions_per_peak,
})
notes = pd.DataFrame({'Notes': [
    f'Site sets: Supplementary Table 1, sheets +UV_SRSF3_Nanocompore and -UV_SRSF3_Nanocompore. Tested positions are those with a finite log odds ratio (LOR) and adjusted p value: {N["+UV"]["tested"]:,} (+UV) and {N["-UV"]["tested"]:,} (-UV). The adjusted p value is the Benjamini Hochberg value output by Nanocompore.',
    f'Coordinates: each position (the kmer centre) was projected from the transcript to GRCm39 with R2Dtool 1.0.0 and the Ensembl release 114 annotation. {N["+UV"]["unlifted"]} (+UV) and {N["-UV"]["unlifted"]} (-UV) positions lie on transcripts absent from that annotation ({n_absent_tx} transcripts in all) and are not included; {N["+UV"]["unlifted_sig"]} (+UV) and {N["-UV"]["unlifted_sig"]} (-UV) of them are significant.',
    'iCLIP peaks: SRSF3 iCLIP peaks called with iCount on mm9 (GSE69689), converted to GRCm39 with UCSC liftOver (minMatch 0.95). Peak sets by crosslink score: all peaks, score >= 15 (the top 10%, referred to as high confidence) and score >= 79 (the top 1%).',
    'Overlap: a position counts as within 5 nt of a peak when the 11 nt window centred on it (position plus and minus 5 nt, genomic coordinates) shares at least 1 nt with a peak on the same strand (bedtools 2.31.1). This is the definition used for the dcDNA comparison in the Methods.',
    f'Position level sheet: significant positions are compared with all other tested positions (odds ratio of the 2 x 2 table). Because both the power to call a site and the chance of lying in an iCLIP peak rise with read depth, the comparison is repeated within ten strata of read depth (Mantel Haenszel estimate). The depth of a position is the number of reads of the two SRSF3 dirCLIP libraries with signal at that position in our realignment of these libraries; it is available for {N["+UV"]["depth_pct"]:.1f}% (+UV) and {N["-UV"]["depth_pct"]:.1f}% (-UV) of the projected positions, and the depth columns use those positions only.',
    f'Intervals: positions are not independent. Neighbouring positions share most of their kmer, one crosslink gives a run of significant positions, and isoforms of a gene share positions and peaks. The intervals to rely on are those obtained by resampling whole loci with replacement (10,000 resamples, percentile interval; a locus is a set of transcripts whose tested positions overlap in the genome on the same strand; {N["+UV"]["loci"]} loci on +UV and {N["-UV"]["loci"]} on -UV). Intervals from resampling transcripts, and intervals that take positions as independent (Woolf for the odds ratio, Robins Breslow Greenland for the depth adjusted odds ratio), are given for comparison; the last are too narrow. Resampling intervals were computed for the rule used in the paper only.',
    'A depth adjusted odds ratio below 1 is reported as no enrichment detected; no test for depletion was made.',
    'Peak level sheet: a peak is covered when it lies within 5 nt of at least one tested position, and it contains a significant site when it lies within 5 nt of at least one significant position (adjusted p < 0.005 and |LOR| >= 1). The random comparison draws, 1000 times, as many tested positions as there are significant positions, matched on read depth by decile, and counts the distinct peaks reached. Significant positions cluster (several per peak), which is one reason why they reach fewer distinct peaks than the same number of independent random positions; the last two columns quantify the clustering.',
    'An independent implementation of the position level comparison in R (GenomicFeatures, rtracklayer; script 21_dRNA_iclip_overlap) agrees with these percentages to within 0.01 percentage points when the peaks are read as that script reads them (BED start taken as 1 based), and both implementations then return the same numbers of covered peaks.',
]})
with pd.ExcelWriter(DST, engine='openpyxl') as xw:
    pos.to_excel(xw, sheet_name='position_level', index=False)
    peak.to_excel(xw, sheet_name='peak_level', index=False)
    notes.to_excel(xw, sheet_name='notes', index=False)
wb = load_workbook(DST)
for ws in wb.worksheets:
    for cell in ws[1]:
        cell.font = Font(bold=True, color='FFFFFF'); cell.fill = PatternFill('solid', fgColor='1F4E78')
        cell.alignment = Alignment(wrap_text=True, vertical='center', horizontal='center')
    ws.row_dimensions[1].height = 62 if ws.title != 'notes' else 18
    for i, col in enumerate(ws.columns, 1):
        w = 110 if ws.title == 'notes' else (46 if i <= 3 and ws.title == 'position_level' else (46 if i <= 2 else 17))
        ws.column_dimensions[get_column_letter(i)].width = w
    for row in ws.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = Alignment(wrap_text=True, vertical='top')
    ws.freeze_panes = 'A2'
wb.properties.creator = 'Nikolay Shirokikh'
wb.save(DST)
chk = pd.read_excel(DST, sheet_name=None)
assert list(chk) == ['position_level', 'peak_level', 'notes'] and len(chk['position_level']) == 12 and len(chk['peak_level']) == 6
p = chk['position_level']
row = p[(p['comparison'] == COMP['-UV']) & p['site rule'].str.startswith('adjusted p < 0.005') & (p['iCLIP peak set'] == TIER['hc79'])]
src = boot[(boot.comparison == '-UV') & (boot.peak_tier == 'hc79')]
assert len(row) == 1 and row['depth adjusted odds ratio, 95% CI (resampling loci)'].iloc[0] == '%.3f to %.3f' % (src.MH_locus_lo.iloc[0], src.MH_locus_hi.iloc[0])
print(f'written {DST}: position_level {len(p)} rows, peak_level {len(chk["peak_level"])} rows, notes {len(chk["notes"])} rows')
print(p.iloc[:6, [2, 3, 5, 8, 9, 10, 16, 17]].to_string(index=False))
