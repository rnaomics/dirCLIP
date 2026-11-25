# Is the Nanocompore p value column raw or adjusted? Evidence from the source, 5 Oct 2026

Question raised again by a coauthor on 1 Oct 2026. Answer: in every 1.0.x release the column
`GMM_logit_pvalue` of `nanocompore_results.tsv` holds the Benjamini Hochberg adjusted value.

## Files read

| file | origin | lines | md5 |
|---|---|---|---|
| `nanocompore/SampCompDB.py` | PyPI wheel `nanocompore-1.0.4-py3-none-any.whl` | 973 | `0d67450f9242271a166992eb0257045d` |
| `nanocompore/SampCompDB.py` | GitHub `tleonardi/nanocompore`, tag `v1.0.4` | 973 | `0d67450f9242271a166992eb0257045d` |
| `nanocompore/TxComp.py` | PyPI wheel 1.0.4 | 420 | `d7c2b97fc886c1b6a9392765d21b56ff` |
| `nanocompore/TxComp.py` | GitHub tag `v1.0.4` | 420 | `d7c2b97fc886c1b6a9392765d21b56ff` |

The released package and the tagged source are byte identical.

## What the code does (release 1.0.4, line numbers)

1. `TxComp.py` 175 to 229, `gmm_test`: a two component Gaussian mixture (scikit-learn
   `GaussianMixture`, line 16 and `fit_best_gmm`) on scaled intensity and log10 dwell time.
2. `TxComp.py` 294 to 311, `gmm_logit_test`: one pseudo read is added to each condition by cluster
   cell (296, 297), then a logistic regression of cluster on condition
   (`statsmodels.discrete.discrete_model.Logit`, 300 to 305). Its Wald p value is the **raw** p of the
   position, and its coefficient is `Logit_LOR`.
3. `SampCompDB.py` 98 to 101: when the results database is opened, the results table is built by
   `__calculate_results(adjust=True)`. The code comment at line 98 reads "Create results df with
   adjusted p-values".
4. `SampCompDB.py` 169 to 171: with `adjust` set, every p value column is replaced by the output of
   `__multipletests_filter_nan(..., method="fdr_bh")`, under its original name.
5. `SampCompDB.py` 215 to 233: that helper calls `statsmodels.stats.multitest.multipletests`
   (imported at line 23) on the finite values and puts the missing ones back.
6. `SampCompDB.py` 341 to 393, `save_report`: writes the p value columns from that adjusted table
   (374, 379, 380) and takes `Logit_LOR` from the stored model (390). `save_all` (236 to 251) calls it
   to write `nanocompore_results.tsv`, and the `sampcomp` command calls `save_all`
   (`__main__.py` line 218).

Tags `v1.0.1`, `v1.0.2` and `v1.0.3` carry the same two statements (lines 101 and 171 of a 973
line file in each). There is no `v1.0.0` tag, only release candidates and betas.

## Why the development branch reads differently

His link, `nanocompore/comparisons.py` on the master branch, is Nanocompore 2, a rewrite. At master
there is no `SampCompDB.py` and no `TxComp.py`; the test is a chi squared test on a GPU mixture
model, and the output has two columns, `GMM_chi2_pvalue` (raw) and `GMM_chi2_qvalue` (Benjamini
Hochberg; `postprocessing.py`, `_correct_pvalues`; `docs/output.md`). Release 1.0.x, which wrote
`GMM_logit_pvalue` and `Logit_LOR`, has a single column and overwrites it.

His second link, `statsmodels.sandbox.regression.gmm.GMMResults`, is the generalized method of
moments, an unrelated estimator that shares the abbreviation. Nanocompore does not use it.

## Numerical check anyone can repeat

`check_nanocompore_p_is_BH.R` (base R, about 10 s) rebuilds the raw p of every position from
`cluster_counts` and compares. Output on our 2026 rerun, `check_nanocompore_p_is_BH_OUTPUT.txt`:
the rebuilt LOR matches the reported one to 1.2e-13; the reported p equals the raw p at 0.24% of
226,360 positions and the Benjamini Hochberg value of the raw p at 100.00%; the smallest reported
value is the smallest raw p multiplied by 226,360, the number of positions tested. With the raw p
the site rule would give 3,336 positions, not 956.

## Consequences

- Methods: "adjusted p-value <0.005" stays. Figure 1E legend and axis: adjusted p value.
- Parameter table, script 7 row: "Adjusted", not "Raw" (changed in the workflow parameter table on 5 Oct 2026).
- Figure 1F and the read termination row: script 9 gates and plots `adj_pval` of
  `glmGamPoi::test_de` (`pipeline/9`, the `Sig` rule and the volcano axis), so "adjusted" is right there too.
