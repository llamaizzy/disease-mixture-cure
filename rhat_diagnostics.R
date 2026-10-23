## =====================================================================
##  Modern convergence diagnostics (rank-normalised R-hat, bulk/tail ESS)
##  for a fit you have run.  Moved out of run_model.R, which is sourced by
##  the simulation workers and must not run analysis code on load.
##
##      source("run_model.R")
##      fit <- readRDS("results_80k/infer_fit.rds")     # or any fit object
##      source("rhat_diagnostics.R")                     # prints `diag`
## =====================================================================
if (!exists("fit")) stop("define `fit` first, e.g. fit <- readRDS(\"results_80k/infer_fit.rds\")")
#evaluation of modern Rhat
# install.packages("posterior")   # run once, by hand, if posterior is missing
library(posterior)

## ---- remove first third of EACH chain -------------------------------

burn_frac <- 1/3

L <- lapply(fit$chains, function(m) {
  start <- floor(nrow(m) * burn_frac) + 1
  m[start:nrow(m), , drop = FALSE]
})

## sanity checks
stopifnot(
  length(unique(sapply(L, nrow))) == 1,
  length(unique(sapply(L, ncol))) == 1,
  all(sapply(L, function(x)
    identical(colnames(x), colnames(L[[1]]))))
)

## ---- convert list of chain matrices to posterior draws array --------
n_iter  <- nrow(L[[1]])
n_chain <- length(L)
n_var   <- ncol(L[[1]])
arr <- array(
  NA_real_,
  dim = c(n_iter, n_chain, n_var),
  dimnames = list(
    iteration = NULL,
    chain = paste0("chain", seq_len(n_chain)),
    variable = colnames(L[[1]])
  )
)

for (ch in seq_len(n_chain)) {
  arr[, ch, ] <- L[[ch]]
}
draws <- posterior::as_draws_array(arr)

## ---- modern convergence diagnostics ---------------------------------

diag <- posterior::summarise_draws(
  draws,
  rhat = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail
)
diag
