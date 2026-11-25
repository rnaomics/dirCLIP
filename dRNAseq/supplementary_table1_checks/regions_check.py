#!/usr/bin/env python3
"""Supplementary Table 1 region sheets: the duplicated -UV region, neighbouring regions of the two
channels that were not merged, and which published regions the stated rules (signal deviation 10 nt,
read termination 5 nt with log2FC > 2, cross channel 5 nt) do not reproduce.

Usage: python3 regions_check.py <SupT1_tsv_dir>
The merge functions are those of cross_channel_merge_test.py, copied unchanged so that this file
runs on its own.
"""
import math, sys
from collections import defaultdict, Counter
D = sys.argv[1]


def num(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except ValueError:
        return None


def read(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            yield dict(zip(hdr, line.rstrip('\n').split('\t')))


def collapse(rows, gap):
    by = defaultdict(list)
    for ref, pos in rows:
        by[ref].append(pos)
    out = []
    for ref, ps in by.items():
        ps.sort()
        start = prev = ps[0]
        for p in ps[1:]:
            if p - prev > gap:
                out.append((ref, start, prev)); start = p
            prev = p
        out.append((ref, start, prev))
    return [(ref, (s + e) / 2) for ref, s, e in out]


def cross(sig, rt, d):
    by = defaultdict(list)
    for ref, mid in sig: by[ref].append((mid, 'signal_dev'))
    for ref, mid in rt: by[ref].append((mid, 'read_ends'))
    out = []
    for ref, items in by.items():
        items.sort(key=lambda t: t[0])
        grp = [items[0]]

        def flush(g):
            lo = min(m for m, _ in g) - 5; hi = max(m for m, _ in g) + 5
            srcs = []
            for _, s in g:
                if s not in srcs: srcs.append(s)
            out.append((ref, math.floor((lo + hi) / 2), 'read_ends+signal_dev' if len(srcs) == 2 else srcs[0]))
        for it in items[1:]:
            if it[0] - grp[-1][0] > d:
                flush(grp); grp = [it]
            else:
                grp.append(it)
        flush(grp)
    return out


for arm in ('+UV', '-UV'):
    kat = [(r['ref_id'], int(float(r['region_center'])), r['source']) for r in read(f'{D}/{arm}_SRSF3_Sig_Regions.tsv')]
    cnt = Counter(s for _, _, s in kat)
    keys = Counter((r, c) for r, c, _ in kat)
    dups = [k for k, n in keys.items() if n > 1]
    print(f'\n===== {arm}: published rows {len(kat)}; signal_dev {cnt["signal_dev"]}, read_ends {cnt["read_ends"]}, both {cnt["read_ends+signal_dev"]}; '
          f'unique (transcript, centre) {len(keys)}; duplicated {dups}')
    assert len(kat) in (742, 1893)
    # neighbouring published regions of different single channels, centres within 5 nt
    by = defaultdict(list)
    for r, c, s in kat: by[r].append((c, s))
    near = []
    for r, items in by.items():
        items.sort()
        for (c1, s1), (c2, s2) in zip(items, items[1:]):
            if c2 - c1 <= 5:
                near.append((r, c1, s1, c2, s2))
    print(f'  published neighbours with centres <= 5 nt apart: {len(near)}')
    for n in near: print('    ', n)
    # what the stated rules give
    nano = [(r['ref_id'], int(float(r['pos']))) for r in read(f'{D}/{arm}_SRSF3_Nanocompore.tsv')
            if num(r['GMM_adj_pval']) is not None and num(r['LOR']) is not None and num(r['GMM_adj_pval']) < 0.005 and abs(num(r['LOR'])) >= 1]
    rt = [(r['ref_id'], int(float(r['pos']))) for r in read(f'{D}/{arm}_SRSF3_Read_Termination.tsv')
          if num(r['L2FC']) is not None and num(r['adj_pval']) is not None and num(r['adj_pval']) < 0.05 and num(r['L2FC']) > 2]
    pred = cross(collapse(nano, 10), collapse(rt, 5), 5)
    pc = Counter(s for _, _, s in pred)
    K, P = set(kat), set(pred)
    print(f'  stated rules: {len(pred)} regions; signal_dev {pc["signal_dev"]}, read_ends {pc["read_ends"]}, both {pc["read_ends+signal_dev"]}; '
          f'identical to a published region (transcript, centre, source): {len(K & P)} of {len(K)} distinct published regions')
    only_k = sorted(K - P); only_p = sorted(P - K)
    print(f'  published, not reproduced: {len(only_k)}'); [print('    ', x) for x in only_k[:40]]
    print(f'  produced by the rules, not in the published sheet: {len(only_p)}'); [print('    ', x) for x in only_p[:40]]
    if arm == '-UV':
        # the same sheet with the duplicated region relabelled as found by both channels
        fixed = [k for k in kat if not ((k[0], k[1]) in dups)] + [(r, c, 'read_ends+signal_dev') for r, c in dups]
        fc = Counter(s for _, _, s in fixed)
        print(f'  with the duplicate merged and labelled read_ends+signal_dev: {len(fixed)} regions; signal_dev {fc["signal_dev"]}, read_ends {fc["read_ends"]}, both {fc["read_ends+signal_dev"]}; '
              f'identical to the rules: {len(set(fixed) & P)}')
