#!/usr/bin/env python3
# Brushed copy for the repository, 10 Oct 2026: paths are variables or arguments; the logic is as run.
# Revised 18 Nov 2026: notes 2 and 8; HMGA1a site yield from Supplementary Table 3 (optional fifth argument).
"""Supplementary table, direct RNA sequencing depth and site yield (referees 1.3 and 3.3).
Built 27 Sep 2026 from existing outputs, revised 18 Nov 2026; nothing is recomputed from raw data here.
Sources (read only):
  Revision_3/analysis/R1.3_depth/depth_table.tsv          FASTQ read and length statistics (June 2026)
  Revision_3/analysis/R1.3_depth/aligned_read_counts.tsv  primary alignments, minimap2 -ax map-ont -k14 -uf (Kaya)
  Supplementary Table 1 sheets (Revision_4/03_Supplementary/Sup_table1.xlsx), extracted to TSV
  Supplementary Table 3 sheets (HMGA1a), extracted to TSV the same way (optional)
"""
import glob, os, sys
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

# Usage: python3 make_dRNA_depth_table.py <SupT1_tsv_dir> <depth_table.tsv> <aligned_read_counts.tsv> <out_dir> [<SupT3_tsv_dir>]
SUP, DEPTH_TSV, ALN_TSV, OUT = sys.argv[1:5]
SUP3 = sys.argv[5] if len(sys.argv) > 5 else None
depth = pd.read_csv(DEPTH_TSV, sep='\t')
aln = pd.read_csv(ALN_TSV, sep='\t')
meta = {  # FASTQ name: (RBP, role, flow cell)
    'SRSF3_DIRCLIP_dRNA_rep1': ('SRSF3', 'dirCLIP IP', 'FAQ69228'),
    'SRSF3_DIRCLIP_dRNA_rep2': ('SRSF3', 'dirCLIP IP', 'FAR36745'),
    'SRSF3-EGFP_BAC_Total_RNA_UV_dRNA_rep1': ('SRSF3', 'Transcriptome, +UV', 'FAR36802'),
    'SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep1': ('SRSF3', 'Transcriptome, -UV', 'FAS46766'),
    'SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2': ('SRSF3', 'Transcriptome, -UV', 'FAR38530'),
    'HMGA1a_DIRCLIP_iv-polyA_dRNA_rep1': ('HMGA1a', 'dirCLIP IP (in vitro polyA)', 'FAO06310'),
    'HMGA1a_DIRCLIP_iv-polyA_dRNA_rep2': ('HMGA1a', 'dirCLIP IP (in vitro polyA)', 'FAO07922'),
    'HMGA1a_DIRCLIP_iv-polyA_dRNA_rep3': ('HMGA1a', 'dirCLIP IP (in vitro polyA)', 'FAQ69580'),
    'WT_Total_RNA_non-UV_iv-polyA_dRNA_rep1': ('WT', 'Transcriptome, -UV (in vitro polyA)', 'FAQ68836'),
}
order = list(meta)
d = depth.set_index('sample').loc[order].reset_index()
a = aln.set_index('sample').loc[order].reset_index()
assert len(d) == 9 and len(a) == 9

raw = pd.DataFrame({'Sample': d['sample'], 'RBP': [meta[s][0] for s in d['sample']], 'Library': [meta[s][1] for s in d['sample']],
                    'Flow_cell': [meta[s][2] for s in d['sample']], 'Number_of_reads': d.reads,
                    'Total_bases': d.total_bases, 'Mean_read_length': d.mean_len, 'Median_read_length': d.median_len,
                    'N50_read_length': d.N50, 'Max_read_length': d.max_len})
bam = pd.DataFrame({'Sample': a['sample'], 'Reference': a.reference, 'Number_of_reads_primary_aligned': a.primary_aligned,
                    'Percent_of_FASTQ_reads': (a.primary_aligned / d.reads.values * 100).round(1)})

def sheet_stats(arm):
    n = pd.read_csv(f'{SUP}/{arm}_SRSF3_Nanocompore.tsv', sep='\t', dtype=str)
    r = pd.read_csv(f'{SUP}/{arm}_SRSF3_Read_Termination.tsv', sep='\t', dtype=str)
    s = pd.read_csv(f'{SUP}/{arm}_SRSF3_Sig_Regions.tsv', sep='\t', dtype=str)
    p = pd.to_numeric(n.GMM_adj_pval, errors='coerce'); l = pd.to_numeric(n.LOR, errors='coerce')
    fin = p.notna() & l.notna()
    rp = pd.to_numeric(r.adj_pval, errors='coerce'); rl = pd.to_numeric(r.L2FC, errors='coerce')
    vc = s.source.value_counts()
    uniq = len(s.drop_duplicates(['ref_id', 'region_center']))
    return {
        'Comparison': f'SRSF3 dirCLIP IP (rep1, rep2) against SRSF3 transcriptome {arm}',
        'Transcripts with a tested position (signal deviation)': int(n.loc[fin, 'ref_id'].nunique()),
        'Positions tested (signal deviation, finite LOR and p)': int(fin.sum()),
        'Significant positions (|LOR| >= 1, adjusted p < 0.005)': int((fin & (p < 0.005) & (l.abs() >= 1)).sum()),
        'Positions in read termination sheet': len(r),
        'Significant read termination positions (log2FC > 2, adjusted p < 0.05; see note 6)': int(((rp < 0.05) & (rl > 2)).sum()),
        'Regions, signal deviation only': int(vc.get('signal_dev', 0)),
        'Regions, read termination only': int(vc.get('read_ends', 0)),
        'Regions, both channels': int(vc.get('read_ends+signal_dev', 0)),
        'Regions, total rows': len(s),
        'Regions, unique (transcript, centre)': uniq,
    }
sites = pd.DataFrame([sheet_stats('+UV'), sheet_stats('-UV')]).set_index('Comparison').T.reset_index().rename(columns={'index': 'Quantity'})

def hmga1a_stats(d3):
    # sheets of Supplementary Table 3 as TSV; the U1 sheet repeats rows of the main sheet and is left out
    nano = [f for f in sorted(glob.glob(f'{d3}/*Nanocompore.tsv')) if 'U1' not in os.path.basename(f)]
    term = sorted(glob.glob(f'{d3}/*Read_Terminat*.tsv'))
    assert len(nano) == 1 and len(term) == 1, (nano, term)
    n = pd.read_csv(nano[0], sep='\t', dtype=str); r = pd.read_csv(term[0], sep='\t', dtype=str)
    p = pd.to_numeric(n.GMM_adj_pval, errors='coerce'); l = pd.to_numeric(n.LOR, errors='coerce')
    fin = p.notna() & l.notna()
    sig = n[fin & (p < 0.05) & (l.abs() >= 1)]
    merged = 0  # consecutive significant positions of a transcript at most 10 nt apart form one site
    for _, g in sig.groupby('ref_id'):
        pos = sorted(int(x) for x in g.pos)
        merged += 1 + sum(1 for a, b in zip(pos, pos[1:]) if b - a > 10)
    rp = pd.to_numeric(r.adj_pval, errors='coerce'); rl = pd.to_numeric(r.L2FC, errors='coerce')
    return pd.DataFrame({'Quantity': [
        'Transcripts with a tested position (signal deviation)',
        'Positions tested (signal deviation, finite LOR and p)',
        'Significant positions (|LOR| >= 1, adjusted p < 0.05)',
        'Transcripts with a significant position',
        'Sites (significant positions merged within 10 nt)',
        'Positions in read termination sheet',
        'Significant read termination positions (log2FC > 2, adjusted p < 0.05)'],
        'HMGA1a dirCLIP IP (rep1 to rep3) against WT transcriptome -UV, in vitro polyadenylated': [
        int(n.loc[fin, 'ref_id'].nunique()), int(fin.sum()), len(sig), int(sig.ref_id.nunique()), merged,
        len(r), int(((rp < 0.05) & (rl > 2)).sum())]})
hmga = hmga1a_stats(SUP3) if SUP3 else None

notes = [
    '1. Read statistics were computed from the nine direct RNA FASTQ files as held on 24 Jun 2026 (fastq.zip), counting every read (pass and fail were concatenated at basecalling, so no q score filter was applied).',
    '2. Basecalling. The FASTQ files were basecalled with Guppy 6.1.2+e0556ff (flow cell FLO-MIN106, kit SQK-RNA002, model 2020-09-07_rna_r9.4.1_minion_256_8f8fc47b), except those held for SRSF3-EGFP_BAC_Total_RNA_non-UV_dRNA_rep2 and SRSF3-EGFP_BAC_Total_RNA_UV_dRNA_rep1, which are Dorado 0.2.1 basecalls (rna002_70bps_hac@v3, October 2023) of the same reads. Read identifiers match the MinKNOW records one to one, so read counts are unaffected; the read lengths and aligned read counts of these two libraries are a Dorado measure.',
    '3. Primary aligned reads: minimap2 -ax map-ont -k14 -uf against Ensembl GRCm39, cDNA reference for SRSF3 libraries and cDNA plus ncRNA reference for HMGA1a and WT (in vitro polyadenylated) libraries; counted with samtools view -c -F 0x904. Computed on the UWA Kaya cluster, June 2026.',
    '4. Site yield is taken directly from Supplementary Table 1 (the published SRSF3 analysis). The adjusted p is the Benjamini Hochberg value output by Nanocompore (column GMM_adj_pval).',
    '5. The +UV comparison is the one reported in the main text (Fig. 1D and 1E). Its 742 regions are 359 signal deviation only, 379 read termination only and 4 found by both channels.',
    '6. The read termination sheet in Supplementary Table 1 lists a subset of tested positions; the significant count shown applies the published rule to the deposited L2FC column, which is not exactly the ratio that was gated (see the parameter table), so it can differ slightly from the positions underlying the regions.',
    '7. The -UV region sheet carries one duplicated region (ENSMUST00000042235.15, centre 1637, listed as read_ends and as signal_dev), hence 1,893 rows and 1,892 unique regions.',
    '8. HMGA1a site yield is taken from Supplementary Table 3: signal deviation positions with adjusted p < 0.05 and |LOR| >= 1, as Nanocompore reports them (first base of the kmer), and read termination positions with adjusted p < 0.05 and log2FC > 2. Sites are the significant positions merged within 10 nt, as for SRSF3.',
]

wb = Workbook()
hdr_fill = PatternFill('solid', fgColor='1F4E78'); hdr_font = Font(color='FFFFFF', bold=True)
def write(ws, df, title):
    ws.append([title]); ws['A1'].font = Font(bold=True)
    ws.append(list(df.columns))
    for c in ws[2]:
        c.fill = hdr_fill; c.font = hdr_font; c.alignment = Alignment(wrap_text=True, vertical='top')
    for row in df.itertuples(index=False):
        ws.append(list(row))
    for i, col in enumerate(df.columns, 1):
        width = max(len(str(col)), *(len(str(v)) for v in df.iloc[:, i - 1])) + 2
        ws.column_dimensions[get_column_letter(i)].width = min(width, 60)
ws = wb.active; ws.title = 'fastq_rawread_stats'
write(ws, raw, 'Direct RNA libraries, read statistics from FASTQ (all reads)')
write(wb.create_sheet('bam_read_stats'), bam, 'Direct RNA libraries, primary alignments to the transcriptome')
write(wb.create_sheet('SRSF3_site_yield'), sites, 'SRSF3 direct RNA dirCLIP, tested positions, significant positions and regions (from Supplementary Table 1)')
if hmga is not None:
    write(wb.create_sheet('HMGA1a_site_yield'), hmga, 'HMGA1a direct RNA dirCLIP, tested positions, significant positions and sites (from Supplementary Table 3)')
wn = wb.create_sheet('notes'); wn.append(['Notes']); wn['A1'].font = Font(bold=True)
for n in notes:
    wn.append([n]); wn.cell(row=wn.max_row, column=1).alignment = Alignment(wrap_text=True, vertical='top')
wn.column_dimensions['A'].width = 140
dst = os.path.join(OUT, 'Supplementary_Table_dRNA_depth_and_site_yield_2026-11-18.xlsx')
wb.save(dst)
print(raw.to_string(index=False)); print(bam.to_string(index=False)); print(sites.to_string(index=False))
if hmga is not None:
    print(hmga.to_string(index=False))
print('saved', dst)
