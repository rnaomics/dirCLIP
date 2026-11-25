# check_nanocompore_p_is_BH.R
# Is the GMM_logit_pvalue column of a Nanocompore 1.0.x results table the raw p value of the
# logit test, or its Benjamini Hochberg adjusted value?
#
# Usage: Rscript check_nanocompore_p_is_BH.R <nanocompore_results.tsv or .tsv.gz>
# Base R only. Needs the columns GMM_logit_pvalue, Logit_LOR and cluster_counts.
#
# How it works. For each position Nanocompore fits a two cluster Gaussian mixture, then a logistic
# regression of cluster on condition (TxComp.py, gmm_logit_test, release 1.0.4 lines 294 to 311),
# after adding one pseudo read to each of the four condition by cluster cells (line 296). With one
# binary predictor that regression is the 2 x 2 table of cluster_counts, so its Wald p value can be
# rebuilt exactly: LOR = log(a*d / (b*c)), SE = sqrt(1/a + 1/b + 1/c + 1/d), p = 2 * pnorm(-|LOR/SE|).
# The rebuilt LOR must match the reported Logit_LOR; that is the proof that the table was rebuilt
# correctly. The rebuilt raw p is then compared with the reported column, before and after BH.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("give the path of a nanocompore_results.tsv file")
f <- args[1]
x <- read.delim(if (grepl("\\.gz$", f)) gzfile(f) else f, stringsAsFactors = FALSE)
stopifnot(all(c("GMM_logit_pvalue", "Logit_LOR", "cluster_counts") %in% names(x)))
x$lor <- suppressWarnings(as.numeric(x$Logit_LOR))
x$rep <- suppressWarnings(as.numeric(x$GMM_logit_pvalue))
n_all <- nrow(x)
# Nanocompore adjusts over every position that has a p value. Where the logistic regression did not
# converge it stores p = 1 and LOR "NC" (TxComp.py line 307); such rows stay in the family with raw p = 1.
x <- x[is.finite(x$rep), ]
ok <- is.finite(x$lor)
cat(sprintf("rows %d, with a p value %d, of which with a finite LOR %d\n", n_all, nrow(x), sum(ok)))

# cluster_counts looks like  label1_1:12/3__label1_2:9/4__label2_1:2/20  (sample:cluster0/cluster1)
parts <- strsplit(x$cluster_counts[ok], "__", fixed = TRUE)
cond_of <- function(p) sub("_[0-9]+$", "", sub(":[0-9]+/[0-9]+$", "", p))
lev <- unique(cond_of(parts[[1]]))
if (length(lev) != 2) stop("expected samples of exactly two conditions, found: ", paste(lev, collapse = ", "))
cat("conditions:", lev[1], "and", lev[2], "\n")
one <- function(p) {
  cond <- cond_of(p)
  ab <- sub("^.*:", "", p)
  a <- as.numeric(sub("/.*$", "", ab)); b <- as.numeric(sub("^.*/", "", ab))
  c(sum(a[cond == lev[1]]), sum(b[cond == lev[1]]), sum(a[cond == lev[2]]), sum(b[cond == lev[2]]))
}
m <- t(vapply(parts, one, numeric(4))) + 1          # the pseudo read in each cell
lor <- log(m[, 1] * m[, 4] / (m[, 2] * m[, 3]))
se <- sqrt(rowSums(1 / m))
raw <- rep(1, nrow(x))
raw[ok] <- 2 * pnorm(-abs(lor / se))
bh <- p.adjust(raw, method = "BH")
x <- x[ok, ]; raw <- raw[ok]; bh <- bh[ok]      # compare where the LOR could be rebuilt

close_to <- function(a, b) abs(a - b) <= 1e-6 * pmax(abs(a), abs(b))
cat(sprintf("largest difference between rebuilt and reported |LOR|: %.2e (must be about zero)\n",
            max(abs(abs(lor) - abs(x$lor)))))
cat(sprintf("reported p equals the raw p at            %6.2f%% of positions\n", 100 * mean(close_to(x$rep, raw))))
cat(sprintf("reported p equals the BH adjusted raw p at %6.2f%% of positions\n", 100 * mean(close_to(x$rep, bh))))
cat(sprintf("smallest reported p %.3g; smallest raw p %.3g; ratio %.0f (positions tested: %d)\n",
            min(x$rep), min(raw), min(x$rep) / min(raw), nrow(x)))
sig <- abs(x$lor) >= 1
cat(sprintf("positions with |LOR| >= 1 and p < 0.005: %d on the reported column, %d on the raw p\n",
            sum(sig & x$rep < 0.005), sum(sig & raw < 0.005)))
