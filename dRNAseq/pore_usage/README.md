# Pore usage over the sequencing runs

`make_pore_usage_figure.py` reads the duty time tables of the MinKNOW run reports (one zip per run
in `RUN_REPORTS`, and the unzipped run folder `UV_tsp` of the UV treated transcriptome run under
`KAT_RUN_FILES`), writes the source data (`pore_usage_source_data.csv`, `pore_usage_group_summary.csv`)
and draws the supplementary figure as PDF and PNG. The two CSV files here are those of 27 Sep 2026;
point `OUT` at another folder to keep them.

```bash
RUN_REPORTS=/data/Run_reports KAT_RUN_FILES=/data/Nanopore_run_files OUT=/data/pore_usage_out python3 make_pore_usage_figure.py
```

Needs pandas and matplotlib.
