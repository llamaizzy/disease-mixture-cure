## =====================================================================
##  FIT THE MODEL TO SIMULATED DATA, and score how well it recovers the truth.
##  Run from the package root:
##
##      source("simulated_data/fit_simulated_data.R")
##      simulate_all(n_rep = 10, data_id = "crn-v1")            # the datasets, once
##      fit <- fit_simulated_model(rep = 1, tier = "smoke",     # ~1 min, a check
##                                 run_id = "scratch", data_id = "crn-v1")
##      run_simulation_study(n_sim = 10, tier = "infer",        # 10 fits in parallel
##                           run_id = "2026-10-20_crn", data_id = "crn-v1")
##      recovery_summary(run_id = "2026-10-20_crn")             # the coverage / MSE table
##      plot_sim_fits(run_id = "2026-10-20_crn")                # the fitted figures
##
##  or, for one scenario in one go:
##
##      Rscript simulated_data/fit_simulated_data.R <run_id> <data_id> [tier] [n_sim] [scenario]
##
##  Everything a run produces goes under simulated_data/runs/<run_id>/
##  (fits/, logs/, results/, figures/); see study_paths.R.
## =====================================================================

for (f in sort(list.files("R", full.names = TRUE))) source(f)
source("simulated_data/simulate_data.R")
source("simulated_data/evaluate_simulation.R")
## ---- data -----------------------------------------------------------
## Simulated cohorts in the format of data/adni_amyloid.rds, written by
## simulate_all() to simulated_data/data/<data_id>/. Each file carries the data
## AND the truth that generated it (population parameters and per-subject
## latent quantities), which is what the recovery scores are computed against.
##
##      sim <- load_simulated("baseline", rep = 1, data_id = "v1")   # $dat, $truth, $latent

## ---- subjects who start off the grid -------------------------------
##  At the initial values, a subject's accumulator trajectory can run past the
##  top of the SUVR grid by their last scan (a high starting level x_tilde
##  together with a fast speed). The likelihood is then -Inf and the model
##  cannot be initialised. Those subjects are DROPPED from the fit rather than
##  having their inits moved, and the drop is recorded in fit$dropped: how
##  many, what share of the cohort, their ids, and how many of them were true
##  accumulators. They are removed from `latent` too, so the recovery scores
##  compare the fitted subjects with their own truth.
##
##  room = 0 drops only subjects who are strictly off the grid; room > 0 also
##  drops those within that fraction of the clock below the top.
offgrid_at_inits <- function(spec, dat, bA, S = SETTINGS, room = 0) {
  ini <- spec$inits
  r0  <- rate_grid(bA, ini$theta); G0 <- clock_grid(bA, r0, S$thresh)
  Gb0 <- sa_interp(bA$y, G0, ini$mu_b)                   # clock value at the floor
  rb0 <- max(sa_interp(bA$y, r0, ini$mu_b), S$r_min)
  de  <- as.numeric(dat$X %*% ini$psi) + ini$sigma_delta * ini$z
  g0  <- ifelse(ini$x_tilde >= ini$mu_b, sa_interp(bA$y, G0, ini$x_tilde),
                Gb0 + (ini$x_tilde - ini$mu_b) / rb0)
  last <- dat$age_mat[cbind(seq_len(dat$n_subj), dat$n_obs)]
  top  <- G0[bA$n] - room * (G0[bA$n] - Gb0)            # clock value at the top of the grid
  g0_max <- top - exp(de) * (last - dat$a_tilde)        # highest start that stays on the grid
  which(g0 > g0_max)
}

## keep subjects k (in that order) in a dat list and its matching latent list
subset_subjects <- function(dat, latent, k) {
  dat <- within(dat, {
    n_obs <- n_obs[k]; y_mat <- y_mat[k, , drop = FALSE]
    age_mat <- age_mat[k, , drop = FALSE]; a_tilde <- a_tilde[k]
    x_init <- x_init[k]; X <- X[k, , drop = FALSE]; ids <- ids[k]
    n_subj <- length(k) })
  list(dat = dat, latent = lapply(latent, function(v) v[k]))
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
##  The fit goes to runs/<run_id>/fits/ with a stamp of the code, settings and
##  data that produced it (fit$stamp; see study_paths.R).
TIER_CFG <- list(smoke  = list(n_sub = 250L, n_iter =   300L, chunk =  150L),
                 struct = list(n_sub =  NA_integer_, n_iter =  800L, chunk = 400L),
                 infer  = list(n_sub =  NA_integer_, n_iter = 10000L, chunk = 1000L))

fit_simulated_model <- function(scenario = "baseline", rep = 1L,
                                tier = c("smoke", "struct", "infer"),
                                n_chains = 3L, seed = 100L, n_iter = NULL,
                                run_id, data_id = NULL, verbose = TRUE) {
  tier <- match.arg(tier)
  cfg <- TIER_CFG[[tier]]
  if (!is.null(n_iter)) cfg$n_iter <- as.integer(n_iter)
  data_id <- init_run(run_id, data_id)
  out_dir <- run_dir(run_id, "fits")
  tag <- sprintf("%s_%s_rep%02d", tier, scenario, rep)
  log_file <- run_dir(run_id, "logs", sprintf("%s.log", tag))
  t0 <- Sys.time(); el <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
  say <- function(...) { msg <- sprintf(...); cat(msg, file = log_file, append = TRUE)
    if (verbose) { cat(msg); flush.console() } }

  sim <- load_simulated(scenario, rep, data_id)
  stamp <- fit_stamp(sim, tier, n_chains, cfg$n_iter, seed)
  dat <- sim$dat; latent <- sim$latent[c("Z", "alpha", "delta", "z", "b", "x_tilde")]
  if (!is.na(cfg$n_sub) && cfg$n_sub < dat$n_subj) {     # subsample SUBJECTS, never visits
    set.seed(7); k <- sort(sample(dat$n_subj, cfg$n_sub))
    sub <- subset_subjects(dat, latent, k); dat <- sub$dat; latent <- sub$latent
  }
  bA <- make_rate_basis()
  vA <- true_var_basis()
  say("%s rep %d, tier %s: %d subjects, %d scans, %d chains x %d iterations\n",
      scenario, rep, tier, dat$n_subj, sum(dat$n_obs), n_chains, cfg$n_iter)
  
  spec <- build_model(dat, bA, vA)
  drop <- offgrid_at_inits(spec, dat, bA)                         # strictly off the grid
  near <- setdiff(offgrid_at_inits(spec, dat, bA, room = 0.10), drop) # on it, but within 10% of the top
  dropped <- list(n = length(drop), n_total = dat$n_subj, share = length(drop) / dat$n_subj,
                  ids = dat$ids[drop], n_true_accumulators = sum(latent$Z[drop]),
                  n_near_top = length(near), ids_near_top = dat$ids[near])
  if (length(drop)) {
    sub <- subset_subjects(dat, latent, setdiff(seq_len(dat$n_subj), drop))
    dat <- sub$dat; latent <- sub$latent
    spec <- build_model(dat, bA, vA)
    stopifnot("subjects still off the grid after dropping" =
                length(offgrid_at_inits(spec, dat, bA)) == 0)
  }
  say("dropped %d of %d subjects (%.1f%%) off the grid at the inits (%d true accumulators); fitting %d\n",
      dropped$n, dropped$n_total, 100 * dropped$share, dropped$n_true_accumulators, dat$n_subj)
  say("flagged, kept: %d subjects start within 10%% of the top of the grid\n", dropped$n_near_top)
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
  
  ## An interrupted run keeps its FINISHED chains: each chain is seeded on its
  ## own, so redoing only the unfinished ones gives the same fit.
  chains <- vector("list", n_chains); subj <- vector("list", n_chains)
  partial <- file.path(out_dir, sprintf("%s_partial.rds", tag))
  if (file.exists(partial)) {
    old <- readRDS(partial)
    if (length(stamp_mismatch(old$stamp, stamp))) {     # made by other code: start over
      say("  discarding an interrupted run made with different code or data\n")
      old <- list(chains = list())
    }
    for (ch in seq_len(min(n_chains, length(old$chains))))
      if (NROW(old$chains[[ch]]) == cfg$n_iter) {
        chains[[ch]] <- old$chains[[ch]]; subj[[ch]] <- old$subj[[ch]]
      }
  }
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
      saveRDS(list(chains = chains, subj = subj, stamp = stamp), partial)
      say("  chain %d: %d/%d  (%.1f min)\n", ch, it, cfg$n_iter, el())
    }
  }
  fit <- list(tier = tier, scenario = scenario, rep = rep, chains = chains, subj = subj,
              dat = dat, truth = sim$truth, latent = latent, dropped = dropped, data_seed = sim$seed,
              bA = bA, vA = vA, settings = SETTINGS, minutes = el(),
              run_id = run_id, data_id = data_id, stamp = stamp)
  fit_file <- file.path(out_dir, sprintf("%s_fit.rds", tag))
  saveRDS(fit, fit_file)
  unlink(partial)
  say("done in %.1f min -> %s\n", el(), fit_file)
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
  burn <- function(m) m[(floor(nrow(m) * burn_frac) + 1):nrow(m), , drop = FALSE]
  L <- lapply(fit$chains, function(m) burn(m)[, names(truth), drop = FALSE])
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
  U  <- do.call(rbind, lapply(fit$subj, burn))
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
fit_files <- function(scenario, tier, run_id) {
  d <- run_dir(run_id, "fits")
  files <- list.files(d, sprintf("^%s_%s_rep[0-9]+_fit\\.rds$", tier, scenario), full.names = TRUE)
  if (!length(files)) stop("no fits found for ", tier, " / ", scenario, " in ", d)
  files
}

recovery_summary <- function(scenario = "baseline", tier = "infer", run_id, write = TRUE) {
  files <- fit_files(scenario, tier, run_id)
  out_dir <- run_dir(run_id, "results")
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

## ---- subjects dropped off the grid, across fits -----------------------
##  One row per fit: how many subjects were dropped at the inits, out of how
##  many, how many of them were true accumulators, and how many were kept
##  but flagged for starting within 10% of the top of the grid. Fits made before
##  subjects were dropped (when inits were repaired instead) show NA.
dropped_summary <- function(scenario = "baseline", tier = "infer", run_id, write = TRUE) {
  out_dir <- run_dir(run_id, "results")
  d <- do.call(rbind, lapply(fit_files(scenario, tier, run_id), function(f) {
    fit <- readRDS(f); dr <- fit$dropped
    if (is.null(dr)) dr <- list(n = NA, n_total = NA, n_true_accumulators = NA, n_near_top = NA)
    data.frame(scenario = scenario, rep = fit$rep, n_dropped = dr$n, n_subjects = dr$n_total,
               pct_dropped = 100 * dr$n / dr$n_total, n_true_accumulators = dr$n_true_accumulators,
               n_near_top_kept = dr$n_near_top)
  }))
  if (write) write.csv(d, file.path(out_dir, sprintf("dropped_%s_%s.csv", tier, scenario)), row.names = FALSE)
  cat(sprintf("\n%s / %s: subjects dropped off the grid at the inits, %s total over %d fits (mean %.2f%% per fit); %s more kept but flagged near the top\n",
              scenario, tier, format(sum(d$n_dropped)), nrow(d), mean(d$pct_dropped),
              format(sum(d$n_near_top_kept))))
  print(d, row.names = FALSE)
  invisible(d)
}

## ---- the fitted figures ----------------------------------------------
##  The report figures, drawn from the fits to SIMULATED data, in ONE pdf:
##   * page 1: the rate curve across replicates, where the simulation adds
##     what the real cohort cannot: the TRUE curve, drawn over each band;
##   * then two pages per replicate: fig_aligned() and fig_fit_check() exactly
##     as for the real cohort (R/05_figures.R), so the two can be read side
##     by side.
##  The pages keep their own sizes, so they are drawn separately and joined
##  with qpdf.
plot_sim_fits <- function(scenario = "baseline", tier = "infer", run_id, n_draw = 60) {
  if (!requireNamespace("qpdf", quietly = TRUE)) stop("install.packages(\"qpdf\") to join the pages")
  fig_dir <- run_dir(run_id, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile("simfit_"); dir.create(tmp)
  pages <- file.path(tmp, "rate_curves.pdf")
  bands <- list()
  for (f in fit_files(scenario, tier, run_id)) {
    fit <- readRDS(f); rp <- sprintf("rep%02d", fit$rep)
    cat(sprintf("\n%s: %d of %d subjects are true accumulators\n", rp, sum(fit$latent$Z), fit$dat$n_subj))
    sq <- subject_quantities(fit, n_draw = n_draw)
    pg  <- file.path(tmp, sprintf("%s_%s.pdf", rp, c("aligned", "fit_check")))
    fig_aligned(fit, sq, pg[1])
    fig_fit_check(fit, sq, pg[2])
    pages <- c(pages, pg)
    R  <- exp(fit$bA$B %*% t(pool_draws(fit)[, sprintf("theta[%d]", 1:fit$bA$K)]))
    qd <- quantile(dat_long(fit$dat)$y, c(.02, .98))
    bands[[rp]] <- cbind(.band(R, fit$bA$y), rep = rp, in_data = fit$bA$y >= qd[1] & fit$bA$y <= qd[2])
  }
  B  <- do.call(rbind, bands)
  cu <- truth_curves(fit$truth)
  tr <- data.frame(x = cu$bA$y, y = cu$r)
  rate_panel <- function(p, title, sub)
    p + geom_line(data = tr, aes(x, y, colour = "truth"), inherit.aes = FALSE, linewidth = 0.5) +
    geom_vline(xintercept = to_centiloid(fit$truth$mu_b), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = to_centiloid(SETTINGS$thresh), linetype = "dotted", linewidth = 0.3) +
    scale_colour_manual(values = c("posterior median" = "grey15", truth = "#B03030")) +
    coord_cartesian(ylim = c(0, 0.05)) +
    labs(title = title, subtitle = sub, x = "amyloid level y (Centiloid)", y = "SUVR per year") +
    .thm() + theme(legend.position = "bottom")
  B$x <- to_centiloid(B$x); tr$x <- to_centiloid(tr$x)
  p_all <- rate_panel(ggplot(B, aes(x, mid, group = rep)) +
                        geom_line(aes(colour = "posterior median"), linewidth = 0.3, alpha = 0.7),
                      "A. Rate curve r(y): all replicates",
                      "One posterior median per simulated dataset, against the truth.\nDashed = the true floor mu_b; dotted = threshold.")
  p_rep <- rate_panel(ggplot(B, aes(x, mid)) +
                        geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey25", alpha = 0.2) +
                        geom_ribbon(data = B[!B$in_data & B$x < 0, ], aes(ymin = -Inf, ymax = Inf), fill = "grey60", alpha = 0.2) +
                        geom_ribbon(data = B[!B$in_data & B$x > 0, ], aes(ymin = -Inf, ymax = Inf), fill = "grey60", alpha = 0.2) +
                        geom_line(aes(colour = "posterior median"), linewidth = 0.4) + facet_wrap(~rep, nrow = 1),
                      "B. Each replicate, with its 95% band",
                      "Grey = outside the range of that dataset, where the curve is prior-driven.")
  pdf(pages[1], width = 15, height = 4.3); print(p_all + p_rep + plot_layout(widths = c(1, 3.6))); dev.off()
  out <- file.path(fig_dir, sprintf("fit_%s_%s.pdf", tier, scenario))
  qpdf::pdf_combine(pages, out)
  cat(sprintf("\nwrote %s (%d pages)\n", out, length(pages)))
  invisible(B)
}

## ---- the study -------------------------------------------------------
##  Simulate n_sim replicates of one scenario and fit each. The fits run in
##  PARALLEL, one R process per replicate (NIMBLE compiles into the process's
##  own temp directory, so separate processes cannot collide).
##
##  A replicate whose fit is already in the run is skipped ONLY if its stamp
##  matches what this call would produce (same model code, settings, data,
##  tier, chains, iterations, seed), so an interrupted study resumes. A fit
##  that does not match stops the study: a run never mixes code versions.
##  Start a new run_id instead, or delete that fit on purpose.
run_simulation_study <- function(n_sim = 10L, scenario = "baseline",
                                 tier = c("infer", "struct", "smoke"), n_chains = 3L,
                                 n_iter = NULL, cores = min(n_sim, parallel::detectCores() - 2L),
                                 seed = 100L, run_id, data_id = NULL, notes = "") {
  tier <- match.arg(tier)
  data_id <- init_run(run_id, data_id, notes)
  simulate_all(n_rep = n_sim, which = scenario, data_id = data_id, verbose = FALSE)
  n_it <- if (is.null(n_iter)) TIER_CFG[[tier]]$n_iter else as.integer(n_iter)
  files <- run_dir(run_id, "fits", sprintf("%s_%s_rep%02d_fit.rds", tier, scenario, seq_len(n_sim)))
  done <- file.exists(files)
  for (r in which(done)) {
    want <- fit_stamp(load_simulated(scenario, r, data_id), tier, n_chains, n_it, seed)
    bad  <- stamp_mismatch(readRDS(files[r])$stamp, want)
    if (length(bad))
      stop(sprintf("%s in run '%s' was made differently (%s). Use a new run_id, or delete that fit on purpose.",
                   basename(files[r]), run_id, paste(bad, collapse = ", ")))
  }
  todo <- which(!done)
  cat(sprintf("%s / %s, run %s on data %s: %d replicates, %d to fit on %d cores\n",
              scenario, tier, run_id, data_id, n_sim, length(todo), min(cores, length(todo))))
  if (length(todo)) {
    cl <- parallel::makeCluster(min(cores, length(todo)))
    on.exit(parallel::stopCluster(cl))
    parallel::clusterCall(cl, function(wd) { setwd(wd)
      source("simulated_data/fit_simulated_data.R"); NULL }, getwd())
    parallel::clusterApplyLB(cl, todo, function(r, scenario, tier, n_chains, n_iter, seed, run_id) {
      fit_simulated_model(scenario, r, tier, n_chains = n_chains, n_iter = n_iter, seed = seed,
                          run_id = run_id, verbose = FALSE); r
    }, scenario, tier, n_chains, n_iter, seed, run_id)
  }
  dropped_summary(scenario, tier, run_id)
  recovery_summary(scenario, tier, run_id)
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

## ---- traceplots --------------------------------------------------------
##  One panel per parameter, one line per chain, over the FULL chain; the
##  discarded first `burn_frac` of each chain is shaded. A dashed line marks
##  the true value when the fit is of simulated data. R-hat in each panel
##  title is the same split-free version as convergence(), after burn-in.
##
##      fit <- readRDS(run_dir("2026-10-07_v1", "fits", "infer_baseline_rep05_fit.rds"))
##      traceplot_fit(fit)                                   # the key parameters
##      traceplot_fit(fit, pars = sprintf("theta[%d]", 2:10))
##      traceplot_fit(fit, file = run_dir("2026-10-07_v1", "figures", "trace_baseline_rep05.pdf"))
traceplot_fit <- function(fit, pars = c("sigma_alpha", "sigma_delta", "sigma_b", "mu_b",
                                        "gamma[1]", "gamma[2]", "gamma[3]",
                                        "phi[1]", "phi[2]", "phi[3]", "psi[1]", "psi[2]"),
                          burn_frac = 1/3, max_points = 2000, file = NULL,
                          width = 10, height = 7) {
  pars <- intersect(pars, colnames(fit$chains[[1]]))
  n <- nrow(fit$chains[[1]]); keep <- unique(round(seq(1, n, length.out = min(n, max_points))))
  burn <- floor(n * burn_frac)
  r <- convergence_quiet(fit, burn_frac)
  lab <- setNames(sprintf("%s  (R-hat %.3f)", pars, r[pars]), pars)
  d <- do.call(rbind, lapply(seq_along(fit$chains), function(k)
    data.frame(chain = factor(k), iter = rep(keep, length(pars)),
               par = rep(pars, each = length(keep)),
               value = as.numeric(fit$chains[[k]][keep, pars, drop = FALSE]))))
  d$par <- factor(lab[d$par], lab)
  tv <- c(setNames(fit$truth$phi, sprintf("phi[%d]", seq_along(fit$truth$phi))),
          setNames(fit$truth$gamma, sprintf("gamma[%d]", seq_along(fit$truth$gamma))),
          setNames(fit$truth$psi, sprintf("psi[%d]", seq_along(fit$truth$psi))),
          setNames(fit$truth$theta, sprintf("theta[%d]", seq_along(fit$truth$theta))),
          setNames(fit$truth$nu, sprintf("nu[%d]", seq_along(fit$truth$nu))),
          unlist(fit$truth[c("sigma_alpha", "sigma_delta", "sigma_b", "mu_b")]))
  tr <- if (length(tv)) data.frame(par = factor(lab[intersect(pars, names(tv))], lab),
                                   value = tv[intersect(pars, names(tv))]) else NULL
  p <- ggplot(d, aes(iter, value, colour = chain)) +
    annotate("rect", xmin = -Inf, xmax = burn, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
    geom_line(linewidth = 0.25, alpha = 0.8) +
    facet_wrap(~ par, scales = "free_y") +
    labs(title = sprintf("%s rep %s (%s): %d chains x %d iterations",
                         fit$scenario, fit$rep, fit$tier, length(fit$chains), n),
         subtitle = "shaded: discarded as burn-in; dashed: true value",
         x = "iteration", y = NULL) +
    sim_theme()
  if (!is.null(tr) && nrow(tr)) p <- p + geom_hline(data = tr, aes(yintercept = value),
                                                    linetype = "dashed", linewidth = 0.4)
  if (!is.null(file)) ggsave(file, p, width = width, height = height)
  p
}

## R-hat for every monitored parameter, without the printout
convergence_quiet <- function(fit, burn_frac = 1/3) {
  out <- NULL
  capture.output(out <- convergence(fit, burn_frac = burn_frac, n_show = 0))
  out
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("usage: Rscript simulated_data/fit_simulated_data.R <run_id> <data_id> [tier] [n_sim] [scenario]")
  run_id   <- args[1]; data_id <- args[2]
  tier     <- if (length(args) >= 3) args[3] else "infer"
  n_sim    <- if (length(args) >= 4) as.integer(args[4]) else 10L
  scenario <- if (length(args) >= 5) args[5] else "baseline"
  run_simulation_study(n_sim = n_sim, scenario = scenario, tier = tier, run_id = run_id, data_id = data_id)
  plot_sim_fits(scenario = scenario, tier = tier, run_id = run_id)
  evaluate_scenario(scenario, tier = tier, run_id = run_id)
}
