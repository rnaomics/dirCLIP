"""Cross channel merge distance: 4 nt against 5 nt, both arms, exact and overlap scoring.
Independent Python rewrite of the logic in a coauthor's R test (test_kat_param.R, 23 Sep 2026)."""
import math, sys
from collections import defaultdict
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
            f = line.rstrip('\n').split('\t')
            yield dict(zip(hdr, f))

def collapse(rows, gap):
    """rows: list of (ref, pos). R: arrange(ref,pos); group=cumsum(c(0,diff(pos))>gap)."""
    by = defaultdict(list)
    for ref, pos in rows:
        by[ref].append(pos)
    out = []
    for ref, ps in by.items():
        ps.sort()
        start = prev = ps[0]; n = 1
        for p in ps[1:]:
            if p - prev > gap:
                out.append((ref, start, prev, n)); start = p; n = 0
            prev = p; n += 1
        out.append((ref, start, prev, n))
    # region_middle unfloored, as in script 7 and in test_kat_param.R
    return [(ref, (s + e) / 2) for ref, s, e, n in out]

def cross(sig, rt, d):
    by = defaultdict(list)
    for ref, mid in sig: by[ref].append((mid, 'signal_dev'))
    for ref, mid in rt: by[ref].append((mid, 'read_ends'))
    out = []
    for ref, items in by.items():
        items.sort(key=lambda t: t[0])  # stable, like dplyr arrange
        grp = [items[0]]
        def flush(g):
            lo = min(m for m, _ in g) - 5; hi = max(m for m, _ in g) + 5
            srcs = []
            for _, s in g:
                if s not in srcs: srcs.append(s)
            label = 'read_ends+signal_dev' if len(srcs) == 2 else srcs[0]
            out.append((ref, math.floor((lo + hi) / 2), label))
        for it in items[1:]:
            if it[0] - grp[-1][0] > d:
                flush(grp); grp = [it]
            else:
                grp.append(it)
        flush(grp)
    return out

def score(pred, kat):
    """exact: same ref, same centre, same source. overlap: the coauthor's metric, per class,
    fraction of published regions [c-5, c+5] hit by any predicted window of the same class."""
    res = {}
    for cls in ('signal_dev', 'read_ends', 'read_ends+signal_dev'):
        P = [(r, c) for r, c, s in pred if s == cls]
        K = [(r, c) for r, c, s in kat if s == cls]
        Pset = set(P); Kset = set(K)
        exact = len(Pset & Kset)
        byref = defaultdict(list)
        for r, c in P: byref[r].append(c)
        ov = sum(1 for r, c in K if any(abs(pc - c) <= 10 for pc in byref.get(r, ())))
        res[cls] = (len(P), len(K), exact, ov)
    allP = set((r, c) for r, c, s in pred); allK = set((r, c) for r, c, s in kat)
    allPs = set(pred); allKs = set(kat)
    return res, len(allPs & allKs), len(allKs), len(allP & allK)

for arm in ('-UV', '+UV'):
    nano = [(r['ref_id'], int(float(r['pos']))) for r in read(f'{D}/{arm}_SRSF3_Nanocompore.tsv')
            if num(r['GMM_adj_pval']) is not None and num(r['LOR']) is not None
            and num(r['GMM_adj_pval']) < 0.005 and abs(num(r['LOR'])) >= 1]
    rtrows = [r for r in read(f'{D}/{arm}_SRSF3_Read_Termination.tsv')
              if num(r['L2FC']) is not None and num(r['adj_pval']) is not None]
    kat_rows = list(read(f'{D}/{arm}_SRSF3_Sig_Regions.tsv'))
    kat = [(r['ref_id'], int(float(r['region_center'])), r['source']) for r in kat_rows]
    kc = defaultdict(int)
    for _, _, s in kat: kc[s] += 1
    dup = len(kat) - len(set((r, c) for r, c, s in kat))
    print(f'\n=== {arm}   published rows {len(kat)}: ' + ', '.join(f'{k} {v}' for k, v in sorted(kc.items())) + f'   duplicate (ref,centre) pairs {dup}')
    sig = collapse(nano, 10)
    for rule in ('>2', '>=2'):
        rt_pos = [(r['ref_id'], int(float(r['pos']))) for r in rtrows
                  if num(r['adj_pval']) < 0.05 and (num(r['L2FC']) > 2 if rule == '>2' else num(r['L2FC']) >= 2)]
        rt = collapse(rt_pos, 5)
        for d in (3, 4, 5, 6):
            pred = cross(sig, rt, d)
            res, exact_all, nK, exact_nosrc = score(pred, kat)
            parts = []
            for cls, short in (('signal_dev', 'sig'), ('read_ends', 'rt'), ('read_ends+signal_dev', 'both')):
                nP, nK_, ex, ov = res[cls]
                parts.append(f'{short} {nP}/{nK_} exact {ex} ovl {ov} ({ov/nK_*100:.1f}%)')
            print(f'  L2FC{rule:3} d={d}: total {len(pred)} | ' + ' | '.join(parts) +
                  f' || exact (centre+source) {exact_all}/{nK} = {exact_all/nK:.4f}; centre only {exact_nosrc}')
