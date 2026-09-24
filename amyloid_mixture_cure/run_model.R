## =====================================================================
##  RUN THE MODEL.  This is the only file you need to call.
##
##      source("run_model.R")
##      fit <- fit_amyloid_model(tier = "smoke")     # ~2 min,  a check
##      fit <- fit_amyloid_model(tier = "struct")    # ~15 min, real cohort
##      fit <- fit_amyloid_model(tier = "infer")     # ~3 h,    for results
##
##  See README.md for what the tiers mean and what comes back.
## =====================================================================


for (f in sort(list.files("R", full.names = TRUE))) source(f)

## ---- data -----------------------------------------------------------
## 1101 ADNI subjects with 2+ amyloid PET scans, 3393 scans. Covariates are
## APOE4 carriage and sex, both CENTRED (so the intercepts read at the
## cohort average, not at a hypothetical non-carrier male).
load_amyloid_data <- function(path = "data/adni_amyloid.rds") readRDS(path)

## ---- the sampler ----------------------------------------------------
##
##  NIMBLE's default samplers will run this model but mix badly, because the
##  parameters that fight each other are not the ones that share a name.
##  Blocking here is BY ROLE, measured on the posterior:
##
##    onset location : mu_b ~ gamma[1] +0.57,  gamma[1] ~ sigma_alpha -0.36
##    rate           : psi[1] ~ sigma_delta +0.32
##    contrasts      : gamma[2:], phi[2:] -- max 0.23, all mixing well
##
##  and gamma[1] is the only one tied to the CURVE: corr(gamma[1], theta[3])
##  = 0.535, against 0.163 for mu_b. So gamma[1] joins the theta block and
##  mu_b does not. Four configurations were compared at 3 x 2000 iterations
##  from STATIONARY starts -- comparing from the inits ranks burn-in escape
##  rather than mixing, and that mistake was made twice before it was caught:
##
##    A  {gamma1, mu_b, theta}      worst R-hat 1.573
##    B  {gamma1, theta}, mu_b solo               1.190   <- used here
##    C  {mu_b, gamma1}, theta solo               2.216
##    D  B + {psi, sigma_delta}                   1.478
##
##  Two further rules, both learned the hard way:
##   * sigma_delta gets a SLICE sampler, never a conjugate one. The conjugate
##     Inverse-Gamma draw is valid only when delta itself is the stochastic
##     node with prior N(0, sigma_delta). Here delta is non-centred and the
##     stochastic node is z ~ N(0,1), whose density contains no sigma at all.
##   * mu_b gets a plain slice sampler, not AF_slice. It sits behind an -Inf
##     wall (the feasibility constraint) that AF_slice's adaptation cannot see.
configure_samplers <- function(model, dat, bA, vA) {
  mon  <- c("theta", "nu", "mu_b", "gamma", "sigma_alpha", "psi",
            "sigma_delta", "sp_th", "phi", "pi_bar", "sigma_b")
  mon2 <- c("x_tilde", "z", "b")
  conf <- configureMCMC(model, monitors = mon, monitors2 = mon2, thin2 = 20)

  ## Remove the defaults for every variable we are about to reconfigure BEFORE
  ## adding any replacement. removeSampler("gamma") drops every sampler whose
  ## target mentions gamma -- including a block already added that happens to
  ## contain gamma[1] -- so interleaving remove and add silently leaves nodes
  ## unsampled. (conf$getUnsampledNodes() is how you catch this.)
  for (v in c("nu", "theta", "gamma", "psi",
              "sigma_delta", "sigma_alpha", "mu_b", "sigma_b")) conf$removeSampler(v)

  ## each random-walk-smoothed vector in its own adaptive-factor block
  conf$addSampler(sprintf("nu[1:%d]", vA$K), "AF_slice")
  conf$addSampler(c("gamma[1]", sprintf("theta[%d]", 2:bA$K)), "AF_slice")
  conf$addSampler(sprintf("gamma[2:%d]", ncol(dat$X) + 1), "AF_slice")
  conf$addSampler(sprintf("psi[1:%d]", ncol(dat$X)), "AF_slice")

  for (v in c("sigma_delta", "sigma_alpha", "mu_b", "sigma_b")) conf$addSampler(v, "slice")

  stopifnot("every node must have a sampler" = length(conf$getUnsampledNodes()) == 0)

  for (i in seq_len(dat$n_subj)) {
    ## b_i is centred with a hyperparameter-dependent truncation; it belongs
    ## to the OTHER branch from the accumulator pair, so it samples alone.
    conf$removeSampler(sprintf("b[%d]", i))
    conf$addSampler(sprintf("b[%d]", i), "slice")
    ## a subject's placement and speed trade off directly against each other
    ## (start late and go fast, or start early and go slow), so they move together
    conf$removeSampler(sprintf("x_tilde[%d]", i))
    conf$removeSampler(sprintf("z[%d]", i))
    conf$addSampler(c(sprintf("x_tilde[%d]", i), sprintf("z[%d]", i)), "RW_block")
  }
  conf
}

## ---- the entry point -------------------------------------------------
fit_amyloid_model <- function(tier = c("smoke", "struct", "infer"),
                              n_chains = 3L, seed = 100L,
                              out_dir = "results", verbose = TRUE) {
  tier <- match.arg(tier)
  cfg <- switch(tier,
    smoke  = list(n_sub = 250L, n_iter =   300L, chunk =  150L),
    struct = list(n_sub =  NA_integer_, n_iter =  800L, chunk = 400L),
    infer  = list(n_sub =  NA_integer_, n_iter = 10000L, chunk = 1000L))
  dir.create(out_dir, showWarnings = FALSE)
  t0 <- Sys.time(); el <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
  say <- function(...) if (verbose) { cat(sprintf(...)); flush.console() }

  dat <- load_amyloid_data()
  if (!is.na(cfg$n_sub)) {                    # subsample SUBJECTS, never visits
    set.seed(7); k <- sort(sample(dat$n_subj, cfg$n_sub))
    dat <- within(dat, {
      n_obs <- n_obs[k]; y_mat <- y_mat[k, , drop = FALSE]
      age_mat <- age_mat[k, , drop = FALSE]; a_tilde <- a_tilde[k]
      x_init <- x_init[k]; X <- X[k, , drop = FALSE]; ids <- ids[k]
      n_subj <- length(k) })
  }
  bA <- make_rate_basis()
  vA <- make_var_basis(unlist(lapply(seq_len(dat$n_subj),
                                     function(i) dat$y_mat[i, 1:dat$n_obs[i]])))
  say("tier %s: %d subjects, %d scans, %d chains x %d iterations\n",
      tier, dat$n_subj, sum(dat$n_obs), n_chains, cfg$n_iter)

  spec <- build_model(dat, bA, vA)
  m <- nimbleModel(spec$code, spec$constants, spec$data, spec$inits,
                   calculate = FALSE, check = FALSE)
  lp <- m$calculate()
  say("logProb at inits: %.2f (finite: %s)\n", lp, is.finite(lp))
  if (!is.finite(lp)) {
    nl  <- m$getLogProb(m$getNodeNames(stochOnly = TRUE))
    bad <- m$getNodeNames(stochOnly = TRUE)[!is.finite(nl)]
    stop(sprintf("%d non-finite nodes at the inits; first few: %s",
                 length(bad), paste(head(bad, 5), collapse = ", ")))
  }
  cm  <- compileNimble(m)
  cmc <- compileNimble(buildMCMC(configure_samplers(m, dat, bA, vA)), project = m)
  say("compiled (%.1f min)\n", el())

  chains <- vector("list", n_chains); subj <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    set.seed(seed + ch)
    ini <- spec$inits
    if (ch > 1) {
      ## jitter the start, shrinking until the state is FEASIBLE. Chains that
      ## all start at the same point make R-hat meaningless.
      for (scl in c(1, .5, .25, .125, .0625)) {
        tryv <- ini
        for (nm in c("theta", "nu", "psi", "gamma"))
          tryv[[nm]] <- ini[[nm]] + scl * 0.05 * rnorm(length(ini[[nm]]))
        tryv$theta[1] <- SETTINGS$theta1          # pinned: must match the data
        tryv$mu_b <- ini$mu_b + scl * 0.02 * rnorm(1)
        for (nm in intersect(names(tryv), m$getVarNames())) cm[[nm]] <- tryv[[nm]]
        if (is.finite(cm$calculate())) { ini <- tryv; break }
      }
    }
    for (nm in intersect(names(ini), m$getVarNames())) cm[[nm]] <- ini[[nm]]
    say("  chain %d start logProb %.2f\n", ch, cm$calculate())
    it <- 0
    while (it < cfg$n_iter) {
      n <- min(cfg$chunk, cfg$n_iter - it)
      cmc$run(n, reset = (it == 0), progressBar = FALSE); it <- it + n
      chains[[ch]] <- as.matrix(cmc$mvSamples)
      subj[[ch]]   <- as.matrix(cmc$mvSamples2)
      saveRDS(chains, file.path(out_dir, sprintf("%s_chains.rds", tier)))
      saveRDS(subj,   file.path(out_dir, sprintf("%s_subj.rds", tier)))
      say("  chain %d: %d/%d  (%.1f min)\n", ch, it, cfg$n_iter, el())
    }
  }
  fit <- list(tier = tier, chains = chains, subj = subj, dat = dat,
              bA = bA, vA = vA, settings = SETTINGS, minutes = el())
  saveRDS(fit, file.path(out_dir, sprintf("%s_fit.rds", tier)))
  say("done in %.1f min -> %s\n", el(), file.path(out_dir, sprintf("%s_fit.rds", tier)))
  invisible(fit)
}

## ---- the posterior from the reported fit ----------------------------
##
##  The results in the report come from 3 chains x 10,000 iterations, which
##  takes about three hours. Those draws ship with the package so the figures
##  can be reproduced immediately. Fit it yourself with
##  fit_amyloid_model("infer") when you want to change something.
##
##  Convergence, honestly stated: R-hat falls monotonically with chain length
##  (worst 1.560 -> 1.141 -> 1.097 -> 1.073 at 5k / 10k / 15k / 20k draws).
##  At the reported length 12 of 32 parameters still exceed 1.01. That is a
##  sampler that needs longer, not a broken one -- roughly 40,000-50,000
##  iterations per chain would clear it.
load_fitted_model <- function(path = "data/fitted_posterior.rds") {
  post <- readRDS(path)
  dat  <- load_amyloid_data()
  bA   <- make_rate_basis()
  vA   <- make_var_basis(unlist(lapply(seq_len(dat$n_subj),
                                       function(i) dat$y_mat[i, 1:dat$n_obs[i]])))
  c(post, list(dat = dat, bA = bA, vA = vA, settings = SETTINGS, tier = "infer"))
}

## ---- the headline numbers -------------------------------------------
## Covariate effects on their natural scales. The AFT coefficients are TIME
## RATIOS and the susceptibility coefficients are ODDS RATIOS; reporting
## either on the log scale is how they get misread.
summarise_fit <- function(fit) {
  s <- pool_draws(fit)
  q <- function(v, f = identity) { x <- f(s[, v]); sprintf("%7.3f [%6.3f, %6.3f]",
        median(x), quantile(x, .025), quantile(x, .975)) }
  cat("\n--- susceptibility: does this subject accumulate at all? (odds ratio)\n")
  cat(sprintf("  APOE4          %s\n", q("phi[2]", exp)))
  cat(sprintf("  female         %s\n", q("phi[3]", exp)))
  cat(sprintf("  pi at average  %s\n", q("pi_bar")))
  cat("\n--- latency: when do they start? (time ratio; < 1 means earlier)\n")
  cat(sprintf("  APOE4          %s\n", q("gamma[2]", exp)))
  cat(sprintf("  female         %s\n", q("gamma[3]", exp)))
  cat("\n--- rate: once started, how fast? (rate ratio)\n")
  cat(sprintf("  APOE4          %s\n", q("psi[1]", exp)))
  cat(sprintf("  female         %s\n", q("psi[2]", exp)))
  cat("\n--- the shared structure\n")
  for (v in c("mu_b", "sigma_b", "sigma_alpha", "sigma_delta"))
    cat(sprintf("  %-14s %s\n", v, q(v)))
  invisible(NULL)
}

## ---- convergence -----------------------------------------------------
##  Split-free Gelman-Rubin R-hat across chains, after discarding the first
##  third of EACH chain. Trimming the pooled matrix instead would throw away
##  most of one chain and none of another -- a mistake that once produced a
##  confidently wrong set of "final" numbers here.
convergence <- function(fit, burn_frac = 1/3, n_show = 8) {
  if (isTRUE(fit$pooled)) stop("R-hat needs the separate chains; use a fit you ran yourself.")
  L <- lapply(fit$chains, function(m) m[(floor(nrow(m) * burn_frac) + 1):nrow(m), , drop = FALSE])
  n <- nrow(L[[1]])
  M <- sapply(L, colMeans); V <- sapply(L, function(x) apply(x, 2, var))
  B <- n * apply(M, 1, var); W <- rowMeans(V)
  r <- sqrt(((n - 1) / n * W + B / n) / W)
  r <- r[is.finite(r)]
  cat(sprintf("%d chains x %d draws: worst R-hat %.3f, %d of %d above 1.01\n",
              length(L), n, max(r), sum(r > 1.01), length(r)))
  print(round(sort(r, decreasing = TRUE)[seq_len(min(n_show, length(r)))], 3))
  invisible(r)
}

