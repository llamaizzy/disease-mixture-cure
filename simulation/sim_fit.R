## =====================================================================
##  Fitting one simulated cohort on a model compiled ONCE.
##
##  WHERE THE TIME GOES, and how this file avoids most of it. A fit costs
##  ~70 s to compile and 0.117 s per iteration on the full cohort (M3 Pro,
##  one core). Compiling per replicate would add more than a minute to
##  each; instead every worker compiles one model for the tier's design and
##  then, for each replicate, writes the simulated scans into the compiled
##  model's data node, resets the initial values, and reruns the SAME
##  compiled MCMC. Nothing about the sampler changes between replicates:
##  the design (visit ages, scan counts, covariates) is constant, and it is
##  the only thing baked in at compile time.
##
##  The chains start exactly where run_model.R would start them on real
##  data -- from the data-driven inits in build_model(), jittered for
##  chains 2+ -- and NOT from the truth. Starting at the truth would hide
##  burn-in problems and flatter coverage.
## =====================================================================

## Compile once for a design. Returns the uncompiled model (for variable
## names), the compiled model and the compiled MCMC.
sim_compile <- function(dat, bA, vA) {
  spec <- build_model(dat, bA, vA)
  m <- nimbleModel(spec$code, spec$constants, spec$data, spec$inits,
                   calculate = FALSE, check = FALSE)
  stopifnot("design must initialise to a finite log-density" = is.finite(m$calculate()))
  cm  <- compileNimble(m)
  cmc <- compileNimble(buildMCMC(configure_samplers(m, dat, bA, vA)), project = m)
  list(m = m, cm = cm, cmc = cmc, dat = dat, bA = bA, vA = vA)
}

## Feasible initial values for a simulated cohort.
##
## build_model()'s data-driven inits start each subject at their mean scan
## level. Real ADNI subjects never average above 1.53 SUVR, but the fitted
## model puts more mass in the top of the range than the data show (its
## simulated 99th percentile is 1.53 against 1.47 observed), and a t_4 outlier
## can push a two-scan mean higher still. A subject started near the top of
## the level grid has no room left for the visits after their centroid: the
## clock runs off the grid and the likelihood is -Inf before the first
## iteration. That is a property of the INITIAL rate curve, which is flat
## above the threshold, not of the posterior; so the fix is to start those few
## subjects lower and let the sampler find their level. Everything else is
## build_model()'s own rule. This changes where the chain begins, not what it
## targets.
feasible_inits <- function(W, dat_sim, step = 0.1, max_rounds = 15L) {
  vn <- W$m$getVarNames(); n <- dat_sim$n_subj; max_n <- ncol(dat_sim$y_mat)
  for (round in seq_len(max_rounds)) {
    spec <- build_model(dat_sim, W$bA, W$vA)
    for (nm in intersect(names(spec$inits), vn)) W$cm[[nm]] <- spec$inits[[nm]]
    if (is.finite(W$cm$calculate())) return(spec$inits)
    lp  <- vapply(seq_len(n), function(i) W$cm$getLogProb(sprintf("y_obs[%d, 1:%d]", i, max_n)), 0)
    bad <- which(!is.finite(lp))
    if (!length(bad)) break                      # -Inf is somewhere else: fall through to the error
    dat_sim$x_init[bad] <- pmax(dat_sim$x_init[bad] - step, SETTINGS$thresh)
  }
  stop("could not find feasible initial values after ", max_rounds, " rounds")
}

## Fit one simulated cohort. Returns a fit object with the same fields as
## fit_amyloid_model(), so pool_draws(), subject_quantities() and
## crossing_ages() work on it unchanged.
sim_fit <- function(W, dat_sim, n_chains, n_iter, seed) {
  t0 <- Sys.time()
  W$cm$y_obs <- dat_sim$y_mat                 # the ONLY thing that changes
  inits0 <- feasible_inits(W, dat_sim)        # build_model()'s inits, made feasible
  vn <- W$m$getVarNames()
  chains <- vector("list", n_chains); subj <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    set.seed(seed + ch)
    ini <- inits0
    if (ch > 1) {                             # same jitter rule as run_model.R
      for (scl in c(1, .5, .25, .125, .0625)) {
        tryv <- ini
        for (nm in c("theta", "nu", "psi", "gamma"))
          tryv[[nm]] <- ini[[nm]] + scl * 0.05 * rnorm(length(ini[[nm]]))
        tryv$theta[1] <- SETTINGS$theta1
        tryv$mu_b <- ini$mu_b + scl * 0.02 * rnorm(1)
        for (nm in intersect(names(tryv), vn)) W$cm[[nm]] <- tryv[[nm]]
        if (is.finite(W$cm$calculate())) { ini <- tryv; break }
      }
    }
    for (nm in intersect(names(ini), vn)) W$cm[[nm]] <- ini[[nm]]
    lp <- W$cm$calculate()
    if (!is.finite(lp)) stop(sprintf("chain %d: non-finite log-density at the inits", ch))
    W$cmc$run(n_iter, reset = TRUE, progressBar = FALSE)
    chains[[ch]] <- as.matrix(W$cmc$mvSamples)
    subj[[ch]]   <- as.matrix(W$cmc$mvSamples2)
  }
  list(tier = "sim", chains = chains, subj = subj, dat = dat_sim,
       bA = W$bA, vA = W$vA, settings = SETTINGS,
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
