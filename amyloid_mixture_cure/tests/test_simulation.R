## =====================================================================
##  Self-test for the simulator. Run from the package root:
##
##      Rscript tests/test_simulation.R
##
##  check_dat() already guards each dataset on its own. This checks what it
##  cannot: that the simulator is reproducible, that every scenario is inside
##  the prior support, and that the simulated truth has a FINITE log-density
##  under the actual NIMBLE model -- i.e. the data were generated from the
##  model that will be fitted to them, not from something merely similar.
## =====================================================================
source("simulated_data/fit_simulated_data.R")       # sources R/ and the simulator
ok <- 0L; fail <- character()
chk <- function(name, expr) {
  v <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  if (v) { ok <<- ok + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else   { fail <<- c(fail, name); cat(sprintf("  FAIL  %s\n", name)) }
}

cat("\n1. truth and scenarios\n")
truth <- get_baseline_truth()
sc <- make_scenarios(truth)
chk("truth has every population parameter",
    all(lengths(truth[c("phi", "gamma", "psi", "theta", "nu")]) == c(3, 3, 2, 10, 8)))
chk("theta[1] is the pinned value", truth$theta[1] == SETTINGS$theta1)
chk("eight scenarios, all inside the prior support",
    length(sc) == 8L && all(sapply(names(sc), function(n) check_truth_support(sc[[n]], n))))
chk("each scenario changes exactly one thing",
    all(sapply(sc[-1], function(s) sum(!mapply(identical, s, truth)) == 1L)))
cu <- truth_curves(truth); cu2 <- truth_curves(sc$meas_error_hi)
chk("meas_error_hi doubles omega everywhere", max(abs(cu2$omega / cu$omega - 2)) < 1e-8)
chk("an out-of-support truth is refused",
    inherits(try(check_truth_support(modifyList(truth, list(sigma_delta = 3))), silent = TRUE),
             "try-error"))

cat("\n2. one dataset\n")
real <- readRDS("data/adni_amyloid.rds")
sim  <- simulate_one(truth, seed = 1L, verbose = FALSE)
chk("same seed, same dataset", identical(sim$dat, simulate_one(truth, seed = 1L, verbose = FALSE)$dat))
chk("different seed, different dataset",
    !identical(sim$dat$y_mat, simulate_one(truth, seed = 2L, verbose = FALSE)$dat$y_mat))
chk("same fields and types as the real data",
    identical(names(sim$dat), names(real)) &&
      identical(sapply(sim$dat, typeof), sapply(real, typeof)) &&
      identical(dimnames(sim$dat$X), dimnames(real$X)) &&
      ncol(sim$dat$y_mat) == ncol(real$y_mat))
chk("check_dat rejects a corrupted dataset",
    inherits(try(check_dat(within(sim$dat, y_mat[1, 1] <- NA), sim$latent, verbose = FALSE),
                 silent = TRUE), "try-error"))
## noise-free residuals must be scaled t_4: median |resid / omega| = qt(.75, 4)
obs <- col(sim$dat$y_mat) <= sim$dat$n_obs
res <- (sim$dat$y_mat[obs] - sim$latent$A_mat[obs]) /
  sa_interp(cu$vA$y, cu$omega, sim$latent$A_mat[obs])
chk("standardised residuals are t_4", abs(median(abs(res)) - qt(.75, 4)) < 0.05)
chk("accumulators never go down, non-accumulators are flat",
    all(sapply(seq_len(sim$dat$n_subj), function(i) {
      d <- diff(sim$latent$A_mat[i, 1:sim$dat$n_obs[i]])
      if (sim$latent$Z[i] == 1) all(d >= 0) else all(d == 0) })))

cat("\n3. the fitted model accepts the simulated data\n")
## logProb at the inits the fit starts from, and at the TRUTH, with the var
## basis the truth was written in
lp_at_truth <- function(sim) {
  bA <- make_rate_basis(); vA <- true_var_basis()
  spec <- repair_inits(build_model(sim$dat, bA, vA), sim$dat, bA)
  m <- nimbleModel(spec$code, spec$constants, spec$data, spec$inits,
                   calculate = FALSE, check = FALSE)
  lp0 <- m$calculate()
  for (nm in c("theta", "nu", "mu_b", "gamma", "sigma_alpha", "psi", "sigma_delta",
               "phi", "sigma_b")) m[[nm]] <- sim$truth[[nm]]
  for (nm in c("z", "b", "x_tilde")) m[[nm]] <- sim$latent[[nm]]
  c(inits = lp0, truth = m$calculate())
}
lp <- lp_at_truth(sim)
cat(sprintf("        baseline: logProb %.1f at the inits, %.1f at the truth\n", lp[1], lp[2]))
chk("baseline: finite log-density at the inits", is.finite(lp[["inits"]]))
chk("baseline: finite log-density at the truth", is.finite(lp[["truth"]]))
chk("baseline: the truth beats the inits", lp[["truth"]] > lp[["inits"]])
## the scenario that pushes trajectories hardest against the top of the grid
hi <- simulate_one(sc$speed_spread_hi, seed = 1L, verbose = FALSE)
lp <- lp_at_truth(hi)
cat(sprintf("        speed_spread_hi: logProb %.1f at the inits, %.1f at the truth\n", lp[1], lp[2]))
## build_model()'s own inits are -Inf here (subjects at the top of the grid);
## repair_inits() is what makes this scenario fittable
chk("speed_spread_hi: finite log-density at the inits", is.finite(lp[["inits"]]))
chk("speed_spread_hi: finite log-density at the truth", is.finite(lp[["truth"]]))

cat(sprintf("\n%d passed, %d failed\n", ok, length(fail)))
if (length(fail)) { cat("failed:", paste(fail, collapse = "; "), "\n"); quit(status = 1) }
