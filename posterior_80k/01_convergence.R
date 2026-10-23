## Convergence of the 3 x 80,000 run.
##   * R-hat two ways: the package's classic Gelman-Rubin (as convergence() does),
##     and the modern rank-normalised split R-hat (Vehtari et al. 2021, posterior pkg)
##   * bulk and tail effective sample size
##   * how the worst R-hat falls as the chains get longer, using the first L
##     iterations of each chain and the same 1/3 burn-in rule
source("posterior_80k/00_load.R")
fit <- readRDS("results_80k/infer_fit.rds")
drop <- c("theta[1]")                                    # fixed, not sampled
pars <- setdiff(colnames(fit$chains[[1]]), drop)

classic_rhat <- function(L) {
  n <- nrow(L[[1]]); M <- sapply(L, colMeans); V <- sapply(L, function(x) apply(x, 2, var))
  B <- n * apply(M, 1, var); W <- rowMeans(V); sqrt(((n - 1) / n * W + B / n) / W)
}
diag_at <- function(len, burn_frac = 1/3) {
  L <- lapply(fit$chains, function(m) m[(floor(len * burn_frac) + 1):len, pars, drop = FALSE])
  arr <- array(unlist(L), dim = c(nrow(L[[1]]), length(pars), length(L)))
  arr <- aperm(arr, c(1, 3, 2)); dimnames(arr) <- list(NULL, NULL, pars)
  da <- as_draws_array(arr)
  data.frame(par = pars, len = len, rhat_classic = classic_rhat(L),
             rhat_split = sapply(pars, function(p) rhat(extract_variable_matrix(da, p))),
             ess_bulk = sapply(pars, function(p) ess_bulk(extract_variable_matrix(da, p))),
             ess_tail = sapply(pars, function(p) ess_tail(extract_variable_matrix(da, p))),
             row.names = NULL)
}
lens <- c(10000, 20000, 40000, 60000, 80000)
D <- do.call(rbind, lapply(lens, diag_at))
saveRDS(D, file.path(OUT, "convergence_by_length.rds"))
full <- D[D$len == 80000, ]; full <- full[order(-full$rhat_split), ]
write.csv(full, file.path(OUT, "convergence_80k.csv"), row.names = FALSE)

cat("=== R-hat by chain length (first L iterations of each chain, first third discarded)\n")
cat(sprintf("%7s  %-28s  %-28s  %s\n", "L", "classic: worst / median / >1.01", "split-rank: worst / median / >1.01", "min bulk ESS"))
for (l in lens) { d <- D[D$len == l, ]
  cat(sprintf("%7d  %6.3f / %6.3f / %2d of %d     %6.3f / %6.3f / %2d of %d     %7.0f (%s)\n", l,
      max(d$rhat_classic), median(d$rhat_classic), sum(d$rhat_classic > 1.01), nrow(d),
      max(d$rhat_split), median(d$rhat_split), sum(d$rhat_split > 1.01), nrow(d),
      min(d$ess_bulk), d$par[which.min(d$ess_bulk)])) }
cat("\n=== 80k, all parameters, worst first\n")
print(within(full[, c("par", "rhat_classic", "rhat_split", "ess_bulk", "ess_tail")],
             { rhat_classic <- round(rhat_classic, 4); rhat_split <- round(rhat_split, 4)
               ess_bulk <- round(ess_bulk); ess_tail <- round(ess_tail) }), row.names = FALSE)

## per-chain medians after burn-in, for the worst parameters: is one chain elsewhere?
worst <- head(full$par, 8)
pc <- sapply(fit$chains, function(m) apply(m[26668:80000, worst, drop = FALSE], 2, median))
cat("\n=== per-chain posterior medians (iterations 26,668-80,000), worst parameters\n")
print(round(pc, 4))
## and in thirds of the kept part, to see drift
cat("\n=== chain-pooled medians in three consecutive windows of the kept draws\n")
win <- list(a = 26668:44444, b = 44445:62222, c = 62223:80000)
print(round(sapply(win, function(w) apply(do.call(rbind, lapply(fit$chains, function(m) m[w, worst, drop = FALSE])), 2, median)), 4))
