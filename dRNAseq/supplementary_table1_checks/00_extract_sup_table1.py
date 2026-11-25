# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
"""Stream every sheet of Supplementary Table 1 (xlsx) to TSV. Read only on the source.
Usage: python3 00_extract_sup_table1.py <Sup_table1.xlsx> <out_dir>"""
import zipfile, re, sys, time
from xml.etree.ElementTree import iterparse
SRC = sys.argv[1]
OUT = sys.argv[2]
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
z = zipfile.ZipFile(SRC)
wb = z.read('xl/workbook.xml').decode()
rels = z.read('xl/_rels/workbook.xml.rels').decode()
rid2t = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels)) or {}
for m in re.finditer(r'Target="([^"]+)"[^>]*Id="(rId\d+)"', rels):
    rid2t[m.group(2)] = m.group(1)
sheets = re.findall(r'<sheet [^>]*name="([^"]+)"[^>]*r:id="(rId\d+)"', wb)
# shared strings
ss = []
with z.open('xl/sharedStrings.xml') as fh:
    for ev, el in iterparse(fh, events=('end',)):
        if el.tag == NS + 'si':
            ss.append(''.join(t.text or '' for t in el.iter(NS + 't')))
            el.clear()
print('shared strings', len(ss), flush=True)
def colidx(ref):
    c = re.match(r'([A-Z]+)', ref).group(1); n = 0
    for ch in c: n = n * 26 + ord(ch) - 64
    return n - 1
for name, rid in sheets:
    t0 = time.time()
    target = rid2t[rid].lstrip('/')
    if not target.startswith('xl/'): target = 'xl/' + target
    safe = re.sub(r'[^A-Za-z0-9_+-]', '_', name)
    n = 0
    with z.open(target) as fh, open(f'{OUT}/{safe}.tsv', 'w') as out:
        for ev, el in iterparse(fh, events=('end',)):
            if el.tag == NS + 'row':
                vals = {}
                for c in el.findall(NS + 'c'):
                    v = c.find(NS + 'v'); typ = c.get('t')
                    if typ == 's' and v is not None: val = ss[int(v.text)]
                    elif typ == 'inlineStr':
                        val = ''.join(t.text or '' for t in c.iter(NS + 't'))
                    else: val = v.text if v is not None else ''
                    vals[colidx(c.get('r'))] = val
                if vals:
                    mx = max(vals); out.write('\t'.join(vals.get(i, '') for i in range(mx + 1)) + '\n'); n += 1
                el.clear()
    print(f'{name}: {n} rows, {time.time()-t0:.0f}s', flush=True)
