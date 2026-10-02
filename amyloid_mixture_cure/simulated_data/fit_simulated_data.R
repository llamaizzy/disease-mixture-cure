## =====================================================================
##  FIT THE MODEL TO SIMULATED DATA, and score how well it recovers the truth.
##  Run from the package root:
##
##      source("simulated_data/fit_simulated_data.R")
##      fit <- fit_simulated_model(rep = 1, tier = "smoke")    # ~1 min, a check
##      run_simulation_study(n_sim = 5, tier = "infer")        # ~1.5 h, 5 fits in parallel
##      recovery_summary()                                     # the coverage / MSE table
##
##  or, for the whole study in one go:
##
##      Rscript simulated_data/fit_simulated_data.R [tier] [n_sim]
##
##  Everything is written to simulated_data/results/.
## =====================================================================

for (f in sort(list.files("R", full.names = TRUE))) source(f)
source("simulated_data/simulate_data.R")

## ---- data -----------------------------------------------------------
## Simulated cohorts in the format of data/adni_amyloid.rds, written by
## simulate_all() to simulated_data/datasets/. Each file carries the data
## AND the truth that generated it (population parameters and per-subject
## latent quantities), which is what the recovery scores are computed against.
##
##      sim <- load_simulated("baseline", rep = 1)     # $dat, $truth, $latent

## ---- feasible inits --------------------------------------------------
## build_model() starts x_tilde at the subject's mean scan. For a subject
## sitting near the top of the grid, the init curve then carries the later
## visits OFF the grid and the model starts at -Inf. The real cohort never
## gets that high (max 1.585); a simulated one occasionally does. Start those
## subjects lower, so that the last visit is still on the grid with room to
## spare. Like the z inits, this changes where the chain begins, not what it
## targets.
repair_inits <- function(spec, dat, bA, S = SETTINGS, room = 0.10) {
  ini <- spec$inits
  r0  <- rate_grid(bA, ini$theta); G0 <- clock_grid(bA, r0, S$thresh)
  Gb0 <- sa_interp(bA$y, G0, ini$mu_b)
  rb0 <- max(sa_interp(bA$y, r0, ini$mu_b), S$r_min)
  de  <- as.numeric(dat$X %*% ini$psi) + ini$sigma_delta * ini$z
  g0  <- ifelse(ini$x_tilde >= ini$mu_b, sa_interp(bA$y, G0, ini$x_tilde),
                Gb0 + (ini$x_tilde - ini$mu_b) / rb0)
  last <- dat$age_mat[cbind(seq_len(dat$n_subj), dat$n_obs)]
  top  <- G0[bA$n] - room * (G0[bA$n] - Gb0)            # keep `room` of the clock free
  g0_max <- top - exp(de) * (last - dat$a_tilde)
  fix <- which(g0 > g0_max)
  ini$x_tilde[fix] <- ifelse(g0_max[fix] >= Gb0,
                             sa_inv_monotone(bA$y, G0, pmax(g0_max[fix], Gb0)),
                             ini$mu_b + rb0 * (g0_max[fix] - Gb0))
  spec$inits <- ini
  attr(spec, "n_repaired") <- length(fix)
  spec
}

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
##  Same sampler, tiers and checkpointing as fit_amyloid_model(). Two things
##  differ, both so that the fit can be scored against the truth:
##   * the measurement-scale basis uses the TRUTH's knots (true_var_basis),
##     not the quantiles of this dataset. nu is a coefficient vector on a
##     basis; with different knots it is a different parameter and "does the
##     interval cover the true nu" would have no meaning.
##   * the truth and the latent quantities are saved with the fit.
fit_simulated_model <- function(scenario = "baseline", rep = 1L,
                                tier = c("smoke", "struct", "infer"),
                                n_chains = 3L, seed = 100L, n_iter = NULL,
                                data_dir = "simulated_data/datasets",
                                out_dir = "simulated_data/results", verbose = TRUE) {
  tier <- match.arg(tier)
  cfg <- switch(tier,
                smoke  = list(n_sub = 250L, n_iter =   300L, chunk =  150L),
                struct = list(n_sub =  NA_integer_, n_iter =  800L, chunk = 400L),
                infer  = list(n_sub =  NA_integer_, n_iter = 10000L, chunk = 1000L))
  if (!is.null(n_iter)) cfg$n_iter <- as.integer(n_iter)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tag <- sprintf("%s_%s_rep%02d", tier, scenario, rep)
  log_file <- file.path(out_dir, sprintf("%s.log", tag))
  t0 <- Sys.time(); el <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
  say <- function(...) { msg <- sprintf(...); cat(msg, file = log_file, append = TRUE)
    if (verbose) { cat(msg); flush.console() } }
  
  sim <- load_simulated(scenario, rep, data_dir)
  dat <- sim$dat; latent <- sim$latent[c("Z", "alpha", "delta", "z", "b", "x_tilde")]
  if (!is.na(cfg$n_sub) && cfg$n_sub < dat$n_subj) {     # subsample SUBJECTS, never visits
    set.seed(7); k <- sort(sample(dat$n_subj, cfg$n_sub))
    dat <- within(dat, {
      n_obs <- n_obs[k]; y_mat <- y_mat[k, , drop = FALSE]
      age_mat <- age_mat[k, , drop = FALSE]; a_tilde <- a_tilde[k]
      x_init <- x_init[k]; X <- X[k, , drop = FALSE]; ids <- ids[k]
      n_subj <- length(k) })
    latent <- lapply(latent, function(v) v[k])
  }
  bA <- make_rate_basis()
  vA <- true_var_basis()
  say("%s rep %d, tier %s: %d subjects, %d scans, %d chains x %d iterations\n",
      scenario, rep, tier, dat$n_subj, sum(dat$n_obs), n_chains, cfg$n_iter)
  
  spec <- repair_inits(build_model(dat, bA, vA), dat, bA)
  m <- nimbleModel(spec$code, spec$constants, spec$data, spec$inits,
                   calculate = FALSE, check = FALSE)
  lp <- m$calculate()
  say("logProb at inits: %.2f (finite: %s; %d x_tilde inits lowered)\n",
      lp, is.finite(lp), attr(spec, "n_repaired"))
  if (!is.finite(lp)) {
    nl  <- m$getLogProb(m$getNodeNames(stochOnly = TRUE))
    bad <- m$getNodeNames(stochOnly = TRUE)[!is.finite(nl)]
    stop(sprintf("%d non-finite nodes at the inits; first few: %s",
                 length(bad), paste(head(bad, 5), collapse = ", ")))
  }
  cm  <- compileNimble(m)
  cmc <- compileNimble(buildMCMC(configure_samplers(m, dat, bA, vA)), project = m)
  say("compiled (%.1f min)\n", el())
  
  ## An interrupted run keeps its FINISHED chains: each chain is seeded on its
  ## own, so redoing only the unfinished ones gives the same fit.
  chains <- vector("list", n_chains); subj <- vector("list", n_chains)
  partial <- file.path(out_dir, sprintf("%s_partial.rds", tag))
  if (file.exists(partial)) { old <- readRDS(partial)
    for (ch in seq_len(min(n_chains, length(old$chains))))
      if (!is.null(old$chains[[ch]]) && nrow(old$chains[[ch]]) == cfg$n_iter) {
        chains[[ch]] <- old$chains[[ch]]; subj[[ch]] <- old$subj[[ch]] } }
  for (ch in seq_len(n_chains)) {
    if (!is.null(chains[[ch]])) { say("  chain %d: kept from the interrupted run\n", ch); next }
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
      saveRDS(list(chains = chains, subj = subj), partial)
      say("  chain %d: %d/%d  (%.1f min)\n", ch, it, cfg$n_iter, el())
    }
  }
  fit <- list(tier = tier, scenario = scenario, rep = rep, chains = chains, subj = subj,
              dat = dat, truth = sim$truth, latent = latent, data_seed = sim$seed,
              bA = bA, vA = vA, settings = SETTINGS, minutes = el())
  saveRDS(fit, file.path(out_dir, sprintf("%s_fit.rds", tag)))
  unlink(partial)
  say("done in %.1f min -> %s\n", el(), file.path(out_dir, sprintf("%s_fit.rds", tag)))
  invisible(fit)
}

## ---- recovery: one fit against its truth -----------------------------
##  For every parameter that HAS a true value: the posterior median, the 95%
##  credible interval, whether the interval contains the truth, and R-hat.
##
##  Not scored, because the simulation has no true value for them:
##   * lam_sp_th[j], sp_th -- the smoothing hyperparameters. The data are
##     generated from theta directly, not from a draw of the random walk.
##   * theta[1], alpha_min, alpha_max -- fixed, not estimated.
recovery_one <- function(fit, burn_frac = 1/3) {
  tr <- fit$truth
  truth <- c(setNames(tr$phi, sprintf("phi[%d]", 1:3)),
             setNames(tr$gamma, sprintf("gamma[%d]", 1:3)), sigma_alpha = tr$sigma_alpha,
             setNames(tr$psi, sprintf("psi[%d]", 1:2)), sigma_delta = tr$sigma_delta,
             mu_b = tr$mu_b, sigma_b = tr$sigma_b,
             setNames(tr$theta[-1], sprintf("theta[%d]", 2:10)),
             setNames(tr$nu, sprintf("nu[%d]", 1:8)))
  group <- c(rep("Susceptibility", 3), rep("Departure age", 4), rep("Speed", 3),
             rep("Floor", 2), rep("Curves", 17))
  L <- lapply(fit$chains, function(m) m[(floor(nrow(m) * burn_frac) + 1):nrow(m), names(truth), drop = FALSE])
  s <- do.call(rbind, L); n <- nrow(L[[1]])
  q <- apply(s, 2, quantile, c(.025, .5, .975))
  ## same split-free R-hat as convergence()
  M <- sapply(L, colMeans); V <- sapply(L, function(x) apply(x, 2, var))
  rhat <- sqrt(((n - 1) / n * rowMeans(V) + apply(M, 1, var)) / rowMeans(V))
  out <- data.frame(scenario = fit$scenario, rep = fit$rep, group = group, parameter = names(truth),
                    truth = unname(truth), median = q[2, ], lo = q[1, ], hi = q[3, ],
                    covered = as.numeric(truth >= q[1, ] & truth <= q[3, ]),
                    error = q[2, ] - truth, sq_error = (q[2, ] - truth)^2, rhat = rhat)
  
  ## z_i: one per subject, so its row is the AVERAGE over subjects -- `covered`
  ## is the share of subjects whose interval contains their true z_i.
  U  <- do.call(rbind, lapply(fit$subj, function(m) m[(floor(nrow(m) * burn_frac) + 1):nrow(m), , drop = FALSE]))
  zq <- apply(U[, sprintf("z[%d]", seq_len(fit$dat$n_subj)), drop = FALSE], 2, quantile, c(.025, .5, .975))
  zt <- fit$latent$z
  out <- rbind(out, data.frame(scenario = fit$scenario, rep = fit$rep, group = "Speed",
                               parameter = "z[i] (mean over subjects)", truth = NA,
                               median = NA, lo = NA, hi = NA,
                               covered = mean(zt >= zq[1, ] & zt <= zq[3, ]),
                               error = mean(zq[2, ] - zt), sq_error = mean((zq[2, ] - zt)^2), rhat = NA))
  rownames(out) <- NULL
  out
}

## ---- recovery: across simulations ------------------------------------
##  n_covered : in how many of the n_sim fits the 95% interval held the truth
##  bias      : average of (posterior median - truth)
##  mse       : average of (posterior median - truth)^2;  rmse = sqrt(mse)
##  max_rhat  : worst R-hat over the fits. An interval from a chain that has
##              not converged is not a 95% interval, so read coverage with it.
recovery_summary <- function(scenario = "baseline", tier = "infer",
                             out_dir = "simulated_data/results", write = TRUE) {
  files <- list.files(out_dir, sprintf("^%s_%s_rep[0-9]+_fit\\.rds$", tier, scenario), full.names = TRUE)
  if (!length(files)) stop("no fits found for ", tier, " / ", scenario, " in ", out_dir)
  by_sim <- do.call(rbind, lapply(files, function(f) recovery_one(readRDS(f))))
  by_sim$parameter <- factor(by_sim$parameter, unique(by_sim$parameter))
  tab <- do.call(rbind, lapply(split(by_sim, by_sim$parameter), function(d) {
    scalar <- !is.na(d$truth[1])
    data.frame(group = d$group[1], parameter = as.character(d$parameter[1]), truth = d$truth[1],
               n_sim = nrow(d), n_covered = if (scalar) sum(d$covered) else NA,
               coverage = mean(d$covered), bias = mean(d$error), mse = mean(d$sq_error),
               rmse = sqrt(mean(d$sq_error)), mean_ci_width = mean(d$hi - d$lo),
               max_rhat = if (scalar) max(d$rhat) else NA)
  }))
  rownames(tab) <- NULL
  if (write) {
    write.csv(by_sim, file.path(out_dir, sprintf("recovery_%s_%s_by_sim.csv", tier, scenario)), row.names = FALSE)
    write.csv(tab, file.path(out_dir, sprintf("recovery_%s_%s.csv", tier, scenario)), row.names = FALSE)
  }
  sc <- !is.na(tab$truth)
  cat(sprintf("\n%s / %s: %d fits. 95%% intervals held the truth in %d of %d parameter-fits (%.0f%%).\n\n",
              scenario, tier, length(files), sum(tab$n_covered[sc]), sum(tab$n_sim[sc]),
              100 * sum(tab$n_covered[sc]) / sum(tab$n_sim[sc])))
  print(format(tab, digits = 3), row.names = FALSE)
  invisible(list(summary = tab, by_sim = by_sim))
}

## ---- the study -------------------------------------------------------
##  Simulate n_sim replicates of one scenario and fit each. The fits run in
##  PARALLEL, one R process per replicate (NIMBLE compiles into the process's
##  own temp directory, so separate processes cannot collide). A replicate
##  whose fit is already on disk is skipped, so an interrupted study resumes.
run_simulation_study <- function(n_sim = 5L, scenario = "baseline",
                                 tier = c("infer", "struct", "smoke"), n_chains = 3L,
                                 n_iter = NULL, cores = min(n_sim, parallel::detectCores() - 2L),
                                 out_dir = "simulated_data/results") {
  tier <- match.arg(tier)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  simulate_all(n_rep = n_sim, which = scenario, verbose = FALSE)
  done <- file.exists(file.path(out_dir, sprintf("%s_%s_rep%02d_fit.rds", tier, scenario, seq_len(n_sim))))
  todo <- which(!done)
  cat(sprintf("%s / %s: %d replicates, %d to fit on %d cores (logs in %s)\n",
              scenario, tier, n_sim, length(todo), min(cores, length(todo)), out_dir))
  if (length(todo)) {
    cl <- parallel::makeCluster(min(cores, length(todo)))
    on.exit(parallel::stopCluster(cl))
    parallel::clusterCall(cl, function(wd) { setwd(wd)
      source("simulated_data/fit_simulated_data.R"); NULL }, getwd())
    parallel::clusterApplyLB(cl, todo, function(r, scenario, tier, n_chains, n_iter, out_dir) {
      fit_simulated_model(scenario, r, tier, n_chains = n_chains, n_iter = n_iter,
                          out_dir = out_dir, verbose = FALSE); r
    }, scenario, tier, n_chains, n_iter, out_dir)
  }
  recovery_summary(scenario, tier, out_dir)
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

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_simulation_study(n_sim = if (length(args) >= 2) as.integer(args[2]) else 5L,
                       tier  = if (length(args) >= 1) args[1] else "infer")
}
