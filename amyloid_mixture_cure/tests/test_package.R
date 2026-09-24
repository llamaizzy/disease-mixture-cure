## =====================================================================
##  Self-test. Run from the package root:
##
##      Rscript tests/test_package.R
##
##  It checks the things that, when they break, break silently: the bases
##  matching the shipped posterior, the NIMBLE densities agreeing with their
##  plain-R counterparts, the model initialising to a finite log-density, a
##  short MCMC completing, and every figure drawing.
## =====================================================================
source("run_model.R")
ok <- 0L; fail <- character()
chk <- function(name, expr) {
  v <- tryCatch(isTRUE(expr), error = function(e) { attr(v, "msg") <<- conditionMessage(e); FALSE })
  if (v) { ok <<- ok + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else   { fail <<- c(fail, name); cat(sprintf("  FAIL  %s\n", name)) }
}

cat("\n1. data and bases\n")
dat <- load_amyloid_data()
yall <- unlist(lapply(seq_len(dat$n_subj), function(i) dat$y_mat[i, 1:dat$n_obs[i]]))
bA <- make_rate_basis(); vA <- make_var_basis(yall)
chk("1101 subjects, 3393 scans", dat$n_subj == 1101 && sum(dat$n_obs) == 3393)
chk("no single-visit subjects",  min(dat$n_obs) >= 2)
## Centred against the FULL ADNI covariate table, not against these 1101
## subjects, so the column means are near zero but not exactly zero. The
## shipped posterior's intercepts are on that scale -- do not recentre.
chk("covariates are near-centred", max(abs(colMeans(dat$X))) < 0.05)
## The shipped posterior is written in THESE basis coefficients. A different
## K would reindex theta and every figure would be quietly wrong.
chk("rate basis K = 10",         bA$K == 10L)
chk("variance basis K = 8",      vA$K == 8L)
chk("rate basis rows sum to 1",  max(abs(rowSums(bA$B) - 1)) < 1e-8)

cat("\n2. the clock, R against NIMBLE\n")
th <- seq(-6.9, -3.9, length.out = bA$K)
r  <- rate_grid(bA, th)
G_R <- clock_grid(bA, r, SETTINGS$thresh)
cG  <- compileNimble(nfClock)
G_C <- cG(r, bA$step, bA$n, bA$lo, SETTINGS$thresh, SETTINGS$r_min)
chk("clock agrees to 1e-10",     max(abs(G_R - G_C)) < 1e-10)
chk("clock is zero at threshold", abs(sa_interp(bA$y, G_R, SETTINGS$thresh)) < 1e-10)
chk("clock is strictly increasing", all(diff(G_R) > 0))
## G is invertible: going level -> clock -> level must return the level
ys <- seq(bA$lo + .05, bA$hi - .05, length.out = 50)
chk("clock inverts to 1e-3",
    max(abs(sa_inv_monotone(bA$y, G_R, sa_interp(bA$y, G_R, ys)) - ys)) < 1e-3)

cat("\n3. the implied density of x_tilde\n")
cd <- compileNimble(dxtilde)
## Evaluate at the POSTERIOR MEDIAN curve, not at an arbitrary one. Under a
## curve that is too slow every subject needs an infeasible departure age and
## the density is -Inf everywhere, which would make this test pass vacuously.
gp  <- pool_draws(load_fitted_model())
thm <- apply(gp[, sprintf("theta[%d]", 1:bA$K), drop = FALSE], 2, median)
rm_ <- rate_grid(bA, thm); Gm <- clock_grid(bA, rm_, SETTINGS$thresh)
mb  <- median(gp[, "mu_b"])
args <- list(delta = 0.2, atil = 78, mu_b = mb, xg = log(45), sigma_alpha = 0.29,
             alpha_min = 25, alpha_max = 125, G = Gm, rgrid = rm_,
             lo = bA$lo, step = bA$step, ngrid = bA$n, r_min = SETTINGS$r_min)
ld <- function(x) do.call(cd, c(list(x = x), args, list(log = 1)))
chk("finite inside the grid",     is.finite(ld(0.9)))
chk("-Inf above the grid",        ld(bA$hi + 0.1) == -Inf)
## the extension below mu_b must be CONTINUOUS: without it there is an atom
## at mu_b that a slice sampler cannot reach
eps <- c(1e-2, 1e-3, 1e-4)
gap <- sapply(eps, function(e) abs(ld(mb + e) - ld(mb - e)))
chk("continuous across mu_b (gap decays ~10x per decade)",
    all(gap[-1] / gap[-length(gap)] < 0.25))
## it must integrate to (approximately) one over its support
xs <- seq(0.30, bA$hi - 1e-6, length.out = 4000)
chk("integrates to 1 within 1%",
    abs(sum(exp(sapply(xs, ld))) * diff(xs)[1] - 1) < 0.01)

cat("\n4. the model builds and initialises\n")
spec <- build_model(dat, bA, vA)
m <- nimbleModel(spec$code, spec$constants, spec$data, spec$inits,
                 calculate = FALSE, check = FALSE)
lp <- m$calculate()
chk("finite log-density at the inits", is.finite(lp))
chk("theta[1] is pinned as data", m$theta[1] == SETTINGS$theta1)

cat("\n5. a short MCMC runs\n")
fit <- fit_amyloid_model("smoke", n_chains = 2, out_dir = tempdir(), verbose = FALSE)
s <- pool_draws(fit)
chk("two chains returned",        length(fit$chains) == 2L)
chk("all draws finite",           all(is.finite(s)))
chk("mu_b respects its bounds",   all(s[, "mu_b"] > bA$lo & s[, "mu_b"] < SETTINGS$thresh))
chk("sigma_delta respects its bound", all(s[, "sigma_delta"] < SETTINGS$sigma_delta_max))

cat("\n6. the shipped posterior and the figures\n")
pf <- load_fitted_model()
chk("posterior has 1101 subjects",
    length(grep("^x_tilde", colnames(pf$u))) == 1101)
td <- file.path(tempdir(), "figs")
invisible(make_figures(pf, out_dir = td, n_draw = 12))
chk("three figures written",
    length(list.files(td, pattern = "\\.pdf$")) == 3L)

cat(sprintf("\n%d passed, %d failed\n", ok, length(fail)))
if (length(fail)) { cat("failed:", paste(fail, collapse = "; "), "\n"); quit(status = 1) }
