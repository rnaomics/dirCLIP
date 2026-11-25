#!/usr/bin/env python3
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
"""Supplementary figure, pore usage in the direct RNA runs (referee 3, comment 1).

Final version, 27 Sep 2026. Supersedes Revision_3/analysis/R3.1_pore_usage (June 2026),
which used three runs. This version uses every direct RNA run whose MinKNOW duty time log
survives (8 of 9; SRSF3 -UV transcriptome rep1, FAS46766, has no duty time file).

Inputs are read only: MinKNOW duty_time_*.csv, per minute channel state time summed over
all channels. Productive sequencing = strand + adapter states.
"""
import glob, io, os, re, zipfile
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

RR = os.environ.get('RUN_REPORTS', 'Run_reports')             # directory of MinKNOW run report zips, one per run
KAT = os.environ.get('KAT_RUN_FILES', 'Nanopore_run_files')   # unzipped run folders (UV treated transcriptome run)
OUT = os.environ.get('OUT', os.path.dirname(os.path.abspath(__file__)))

RUNS = [  # label, group, source
    ('SRSF3 -UV rep2', 'Transcriptome', ('zip', 'SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2.zip')),
    ('SRSF3 +UV rep1', 'Transcriptome', ('dir', 'UV_tsp')),
    ('WT -UV (in vitro polyA)', 'Transcriptome', ('zip', 'WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1.zip')),
    ('SRSF3 IP rep1', 'dirCLIP IP', ('zip', 'SRSF3_DIRCLIP_dRNA_rep1.zip')),
    ('SRSF3 IP rep2', 'dirCLIP IP', ('zip', 'SRSF3_DIRCLIP_dRNA_rep2.zip')),
    ('HMGA1a IP rep1', 'dirCLIP IP', ('zip', 'HMGA1a_DIRCLIP_iv-polyA_dRNA_rep1.zip')),
    ('HMGA1a IP rep2', 'dirCLIP IP', ('zip', 'HMGA1a_DIRCLIP_iv-polyA_dRNA_rep2.zip')),
    ('HMGA1a IP rep3', 'dirCLIP IP', ('zip', 'HMGA1a_DIRCLIP_iv-polyA_dRNA_rep3.zip')),
]
GROUP_COL = {'Transcriptome': '#2a78d6', 'dirCLIP IP': '#eb6834'}
STATE_GROUPS = [  # stack order bottom to top; colours validated (dataviz validator, light)
    ('Sequencing, strand', ['strand'], '#4a3aa7'),
    ('Sequencing, adapter', ['adapter'], '#1baf7a'),
    ('Open pore, idle', ['pore'], '#eda100'),
    ('Blocked or recovering', ['unavailable', 'unblocking', 'pending_mux_change', 'pending_manual_reset'], '#e87ba4'),
]
WINDOW_MIN = 360


def load(src):
    kind, name = src
    if kind == 'zip':
        z = zipfile.ZipFile(f'{RR}/{name}')
        member = [n for n in z.namelist() if re.search(r'(^|/)duty_time_.*\.csv$', n)][0]
        raw = z.read(member)
        fs = z.read([n for n in z.namelist() if 'final_summary' in n][0]).decode()
    else:
        member = glob.glob(f'{KAT}/{name}/duty_time_*.csv')[0]
        raw = open(member, 'rb').read()
        fs = open(glob.glob(f'{KAT}/{name}/final_summary_*.txt')[0]).read()
    df = pd.read_csv(io.BytesIO(raw))
    st, mi, sa = df.columns[:3]
    piv = df.pivot_table(index=mi, columns=st, values=sa, aggfunc='sum').fillna(0)
    return piv, re.search(r'flow_cell_id=(\S+)', fs).group(1), os.path.basename(member)


rows, curves, comp = [], {}, {}
for label, grp, src in RUNS:
    piv, fc, member = load(src)
    tot = piv.sum(axis=1)
    prod = (piv.get('strand', 0) + piv.get('adapter', 0)) / tot * 100
    curves[label] = prod.rolling(30, min_periods=1, center=True).mean()
    early = piv[piv.index <= WINDOW_MIN]
    et = early.values.sum()
    shares = {}
    used = set()
    for name, states, _ in STATE_GROUPS:
        shares[name] = sum(early[s].sum() for s in states if s in early.columns) / et * 100
        used.update(states)
    shares['No functional pore or other'] = 100 - sum(shares.values())
    comp[label] = shares
    rows.append(dict(run=label, group=grp, flow_cell=fc, duty_time_file=member,
                     run_hours=round(piv.index.max() / 60, 1),
                     productive_first6h_pct=round(float(prod[prod.index <= WINDOW_MIN].mean()), 1),
                     productive_whole_run_pct=round(float((piv.get('strand', 0).sum() + piv.get('adapter', 0).sum()) / tot.sum() * 100), 1),
                     **{f'first6h_{k}': round(v, 2) for k, v in shares.items()}))
src = pd.DataFrame(rows)
assert len(src) == 8, 'expected eight runs'
src.to_csv(os.path.join(OUT, 'pore_usage_source_data.csv'), index=False)
summ = src.groupby('group')['productive_first6h_pct'].agg(['mean', 'min', 'max', 'count']).round(1)
summ.to_csv(os.path.join(OUT, 'pore_usage_group_summary.csv'))
print(src[['run', 'group', 'flow_cell', 'run_hours', 'productive_first6h_pct', 'productive_whole_run_pct']].to_string(index=False))
print(summ.to_string())

# ---------------- figure ----------------
plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'axes.edgecolor': '#52514e',
                     'xtick.color': '#52514e', 'ytick.color': '#52514e', 'pdf.fonttype': 42})
fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[1.0, 1.4], hspace=0.42, wspace=0.42)
axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, 0]); axC = fig.add_subplot(gs[1, 1])
for ax in (axA, axB, axC):
    ax.spines[['top', 'right']].set_visible(False)

# A. time course
axA.axvspan(0, 6, color='#f0efec', zorder=0)
axA.text(3, 97, 'first 6 h', ha='center', va='top', color='#52514e', fontsize=7)
for label, grp, _ in RUNS:
    c = curves[label]
    axA.plot(c.index / 60, c.values, color=GROUP_COL[grp], lw=1.0, alpha=0.9)
axA.set_xlim(0, 48); axA.set_ylim(0, 100)
axA.set_xlabel('Run time (h)'); axA.set_ylabel('Channels sequencing (%)')
n_t = int((src.group == 'Transcriptome').sum()); n_i = int((src.group == 'dirCLIP IP').sum())
axA.text(35.9, 84, f'Transcriptome, n = {n_t}', color='#0b0b0b', fontsize=8, va='center')
axA.plot([33.6, 35.4], [84, 84], color=GROUP_COL['Transcriptome'], lw=2)
axA.text(35.9, 74, f'dirCLIP IP, n = {n_i}', color='#0b0b0b', fontsize=8, va='center')
axA.plot([33.6, 35.4], [74, 74], color=GROUP_COL['dirCLIP IP'], lw=2)
axA.grid(axis='y', color='#e6e5e0', lw=0.5); axA.set_axisbelow(True)
axA.text(-0.085, 1.06, 'A', transform=axA.transAxes, fontsize=11, fontweight='bold')

# B. first 6 h productive fraction per run, labels repelled vertically with leader lines
xpos = {'Transcriptome': 0, 'dirCLIP IP': 1}
MIN_GAP = 5.2
for grp, g in src.groupby('group'):
    x0 = xpos[grp]
    g = g.sort_values('productive_first6h_pct')
    vals = g.productive_first6h_pct.values
    ys = list(vals)
    for i in range(1, len(ys)):
        if ys[i] - ys[i - 1] < MIN_GAP:
            ys[i] = ys[i - 1] + MIN_GAP
    shift = (np.array(ys) - vals).mean()
    ys = np.array(ys) - shift / 2
    for (i, r), yl in zip(g.iterrows(), ys):
        axB.scatter(x0, r.productive_first6h_pct, s=26, color=GROUP_COL[grp], edgecolor='white', linewidth=0.8, zorder=3)
        axB.plot([x0 + 0.05, x0 + 0.16], [r.productive_first6h_pct, yl], color='#b8b7b1', lw=0.5, zorder=1)
        axB.text(x0 + 0.18, yl, r.run.replace(' (in vitro polyA)', ''), va='center', fontsize=6.3, color='#52514e')
    m = vals.mean()
    axB.plot([x0 - 0.2, x0 - 0.04], [m, m], color='#0b0b0b', lw=1.2, zorder=2)
    axB.text(x0 - 0.22, m, f'{m:.0f}', ha='right', va='center', fontsize=7, color='#0b0b0b')
axB.set_xticks([0, 1]); axB.set_xticklabels([f'Transcriptome\n(n = {n_t})', f'dirCLIP IP\n(n = {n_i})'])
axB.set_xlim(-0.5, 1.95); axB.set_ylim(0, 100)
axB.set_ylabel('Channels sequencing, first 6 h (%)')
axB.grid(axis='y', color='#e6e5e0', lw=0.5); axB.set_axisbelow(True)
axB.text(-0.24, 1.06, 'B', transform=axB.transAxes, fontsize=11, fontweight='bold')

# C. channel state composition, first 6 h
labels = [r[0] for r in RUNS]
x = np.arange(len(labels)) + np.array([0, 0, 0, 0.6, 0.6, 0.6, 0.6, 0.6])
bottom = np.zeros(len(labels))
names = [s[0] for s in STATE_GROUPS] + ['No functional pore or other']
for name in names:
    vals = np.array([comp[l][name] for l in labels])
    if name == 'No functional pore or other':
        axC.bar(x, vals, bottom=bottom, width=0.72, color='white', edgecolor='#b8b7b1', hatch='////', linewidth=0.4)
    else:
        colr = dict((s[0], s[2]) for s in STATE_GROUPS)[name]
        axC.bar(x, vals, bottom=bottom, width=0.72, color=colr, edgecolor='white', linewidth=0.8)
    for xi, b, v in zip(x, bottom, vals):
        if v >= 9 and name != 'No functional pore or other':
            txtcol = 'white' if name == 'Sequencing, strand' else '#0b0b0b'
            axC.text(xi, b + v / 2, f'{v:.0f}', ha='center', va='center', fontsize=6.3, color=txtcol)
    bottom += vals
axC.set_xticks(x); axC.set_xticklabels([l.replace(' (in vitro polyA)', '') for l in labels], rotation=45, ha='right', fontsize=6.8)
axC.set_ylim(0, 100); axC.set_ylabel('Channel time, first 6 h (%)')
axC.text(1, 103, 'Transcriptome', ha='center', fontsize=7, color='#52514e')
axC.text(x[3:].mean(), 103, 'dirCLIP IP', ha='center', fontsize=7, color='#52514e')
handles = [Patch(facecolor=s[2], edgecolor='white', label=s[0]) for s in STATE_GROUPS] + \
          [Patch(facecolor='white', edgecolor='#b8b7b1', hatch='////', label='No functional pore or other')]
axC.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.45, -0.36), ncol=2, frameon=False, fontsize=6.8,
           handlelength=1.2, columnspacing=1.2)
axC.text(-0.16, 1.06, 'C', transform=axC.transAxes, fontsize=11, fontweight='bold')

fig.savefig(os.path.join(OUT, 'Supplementary_Figure_pore_usage.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'Supplementary_Figure_pore_usage.png'), dpi=300, bbox_inches='tight')
print('written', OUT)
