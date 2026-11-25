# Supplementary Table 1: extraction and checks

- `python3 00_extract_sup_table1.py <Sup_table1.xlsx> <TSV_DIR>` (create `<TSV_DIR>` first): writes
  every sheet of the workbook as a TSV of the same name, streaming, read only on the workbook. The
  other analyses read these TSVs.
- `python3 regions_check.py <TSV_DIR>`: applies the stated rules to the position sheets (signal
  deviation: adjusted p < 0.005 and |LOR| >= 1, merged within 10 nt; read termination: adjusted
  p < 0.05 and log2 fold change > 2, merged within 5 nt; the two channels merged within 5 nt; regions
  reported as centre minus 5 to centre plus 5 with the floored centre) and compares the result with
  the region sheets. It prints the duplicated row, the neighbouring regions of the two channels that
  the sheets leave unmerged, and the counts reproduced. The result is summarised in section 7 of the
  top README.
- `python3 cross_channel_merge_test.py <TSV_DIR>`: cross channel merge distances of 3 to 6 nt, both
  arms, under two read termination gates (log2 fold change above 2, and at or above 2), scored exactly
  and by overlap. Independent rewrite of a coauthor's R test.

Both checks need only the Python standard library.
