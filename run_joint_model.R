## =====================================================================
##  The joint amyloid-tau model: fit it.
##
##      source("run_joint_model.R")
##      fit <- fit_joint_model("smoke")     # 250 subjects, 300 iterations
##      fit <- fit_joint_model("struct")    # everyone, 800 iterations
##      fit <- fit_joint_model("infer")     # the real thing, checkpointed
##
##  Everything the amyloid-only runner does (run_model.R) is reused: data,
##  bases, the perturbed-feasible-start scheme, chunked checkpointing. What
##  differs is the model (R/07_joint_model.R) and the sampler blocking.
## =====================================================================
suppressMessages(source("run_model.R"))

## ---- the sampler ----------------------------------------------------
##
##  Blocking follows the amyloid runner, per biomarker:
##   * each random-walk-smoothed vector in its own AF_slice block, with the
##     AFT intercept joining its biomarker's rate-curve block (they move
##     together through the clock scale);
##   * contrasts (gamma, psi, phi) in their own small blocks;
##   * scalars behind -Inf walls (floors, scales, u_rho) on plain slice;
##   * per subject, ONE random-walk block over (x_tilde, z, x_tilde_T, z_T)
##     for subjects with tau scans -- an amyloid proposal moves the implied
##     tau departure age through the exposure clock, so the two pairs are
##     coupled in the posterior -- and two biomarker pairs for subjects
##     without tau, whose tau coordinates are auxiliary. SETTINGS_JOINT$
##     subj_block = "pairs" uses pairs for everyone (the note asks that the
##     two be compared on effective samples per unit time).
##
##  Three changes after the 3 x 25k fit (joint posterior 60k/convergence_report.md):
##   1. beta and sigma_delta_T join the tau rate-curve block. beta's
##      conditional sd was 0.03 against a marginal sd of 0.30, with
##      corr(beta, sigma_delta_T) = 0.63 and corr(beta, theta_T[4:6]) = -0.5:
##      a univariate slice on beta walked the ridge at 1/10 of its width.
##   2. mu_b, sigma_alpha and phi[1] share an AF_slice block (corr -0.84 /
##      0.71 / -0.57 among gamma[1], sigma_alpha, mu_b, phi[1]). CAUTION:
##      run_model.R found mu_b in an AF_slice block with gamma[1] or theta
##      mixed WORSE (its feasibility wall is invisible to the adaptation).
##      This block avoids those partners; it is untested beyond compilation,
##      so SETTINGS_JOINT$amy_block = FALSE restores the plain slice samplers.
##   3. Subject RW blocks start from a proposal covariance on the right
##      scales (x_tilde ~ 0.05 SUVR, z ~ 1, x_tilde_T ~ 0.15) instead of the
##      identity, which NIMBLE warned about on every subject.
configure_joint_samplers <- function(model, dat, bA, vA, bT, vT, SJ = SETTINGS_JOINT) {
  mon  <- c("theta", "nu", "mu_b", "gamma", "sigma_alpha", "psi", "sigma_delta",
            "sp_th", "phi", "pi_bar", "sigma_b",
            "theta_T", "nu_T", "mu_b_T", "gamma_T", "sigma_alpha_T", "psi_T",
            "sigma_delta_T", "sp_th_T", "phi_T", "q_bar", "sigma_b_T",
            "beta", "rho_alpha")
  mon2 <- c("x_tilde", "z", "b", "x_tilde_T", "z_T", "b_T")
  conf <- configureMCMC(model, monitors = mon, monitors2 = mon2, thin2 = 20)

  ## remove every default we replace BEFORE adding (see run_model.R)
  amy_block <- !identical(SJ$amy_block, FALSE)
  for (v in c("nu", "theta", "gamma", "psi", "sigma_delta", "sigma_alpha",
              "mu_b", "sigma_b", "nu_T", "theta_T", "gamma_T", "psi_T",
              "sigma_delta_T", "sigma_alpha_T", "mu_b_T", "sigma_b_T",
              "phi_T", "beta", "u_rho", if (amy_block) "phi")) conf$removeSampler(v)

  p <- ncol(dat$X)
  conf$addSampler(sprintf("nu[1:%d]", vA$K), "AF_slice")
  conf$addSampler(c("gamma[1]", sprintf("theta[%d]", 2:bA$K)), "AF_slice")
  conf$addSampler(sprintf("gamma[2:%d]", p + 1), "AF_slice")
  conf$addSampler(sprintf("psi[1:%d]", p), "AF_slice")
  if (amy_block) {                                       # change 2
    conf$addSampler(c("mu_b", "sigma_alpha", "phi[1]"), "AF_slice")
    conf$addSampler(sprintf("phi[2:%d]", p + 1), "AF_slice")
  } else {
    for (v in c("sigma_alpha", "mu_b")) conf$addSampler(v, "slice")
  }

  conf$addSampler(sprintf("nu_T[1:%d]", vT$K), "AF_slice")
  conf$addSampler(c("gamma_T[1]", sprintf("theta_T[%d]", 1:bT$K),
                    "beta", "sigma_delta_T"), "AF_slice")   # change 1
  conf$addSampler(sprintf("gamma_T[2:%d]", p + 1), "AF_slice")
  conf$addSampler(sprintf("psi_T[1:%d]", p), "AF_slice")
  conf$addSampler(sprintf("phi_T[1:%d]", p + 1), "AF_slice")

  for (v in c("sigma_delta", "sigma_b", "sigma_alpha_T", "mu_b_T", "sigma_b_T",
              "u_rho")) conf$addSampler(v, "slice")

  stopifnot("every node must have a sampler" = length(conf$getUnsampledNodes()) == 0)

  ## change 3: proposal covariances on the coordinates' own scales
  sc <- SJ$subj_scales
  if (is.null(sc)) sc <- c(x_tilde = 0.05, z = 1, x_tilde_T = 0.15, z_T = 1)
  covA <- diag(sc[c("x_tilde", "z")]^2); covT <- diag(sc[c("x_tilde_T", "z_T")]^2)
  cov4 <- diag(sc[c("x_tilde", "z", "x_tilde_T", "z_T")]^2)

  four <- identical(SJ$subj_block, "four")
  for (i in seq_len(dat$n_subj)) {
    for (v in c("b", "b_T")) {
      conf$removeSampler(sprintf("%s[%d]", v, i))
      conf$addSampler(sprintf("%s[%d]", v, i), "slice")
    }
    pairA <- c(sprintf("x_tilde[%d]", i), sprintf("z[%d]", i))
    pairT <- c(sprintf("x_tilde_T[%d]", i), sprintf("z_T[%d]", i))
    for (v in c(pairA, pairT)) conf$removeSampler(v)
    if (four && dat$has_tau[i]) {
      conf$addSampler(c(pairA, pairT), "RW_block", control = list(propCov = cov4))
    } else {
      conf$addSampler(pairA, "RW_block", control = list(propCov = covA))
      conf$addSampler(pairT, "RW_block", control = list(propCov = covT))
    }
  }
  conf
}

## ---- one chain, start to finish ---------------------------------------
##
##  Builds and compiles the model, jitters the start for chains > 1 (shrinking
##  the jitter until the state is feasible; the coupling and correlation are
##  jittered too, so chains do not all start from the uncoupled model), runs
##  the chunks and checkpoints after each. Self-contained so that it can run
##  in a worker process. Returns list(chains = <matrix>, subj = <matrix>).
run_joint_chain <- function(ch, job, verbose = TRUE) {
  spec <- job$spec; cfg <- job$cfg; tier <- job$tier; out_dir <- job$out_dir
  t0 <- Sys.time(); el <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
  logf <- file.path(out_dir, sprintf("%s_chain%d.log", tier, ch))
  say <- function(...) {
    msg <- sprintf(...); cat(msg, file = logf, append = TRUE)
    if (verbose) { cat(msg); flush.console() }
  }
  cat("", file = logf)
  m <- nimbleModel(spec$code, spec$constants, spec$data, spec$inits,
                   calculate = FALSE, check = FALSE)
  lp <- m$calculate()
  say("chain %d: logProb at inits %.2f (finite: %s)\n", ch, lp, is.finite(lp))
  if (!is.finite(lp)) {
    nl  <- sapply(m$getNodeNames(stochOnly = TRUE), function(n) m$getLogProb(n))
    bad <- m$getNodeNames(stochOnly = TRUE)[!is.finite(nl)]
    stop(sprintf("%d non-finite nodes at the inits; first few: %s",
                 length(bad), paste(head(bad, 8), collapse = ", ")))
  }
  cm  <- compileNimble(m)
  cmc <- compileNimble(buildMCMC(configure_joint_samplers(m, job$dat, job$bA, job$vA,
                                                          job$bT, job$vT)), project = m)
  say("chain %d: compiled (%.1f min)\n", ch, el())

  set.seed(job$seed + ch)
  ini <- spec$inits
  if (ch > 1) {
    for (scl in c(1, .5, .25, .125, .0625)) {
      tryv <- ini
      for (nm in c("theta", "nu", "psi", "gamma", "theta_T", "nu_T", "psi_T",
                   "gamma_T", "phi_T"))
        tryv[[nm]] <- ini[[nm]] + scl * job$jitter * rnorm(length(ini[[nm]]))
      tryv$theta[1] <- SETTINGS$theta1
      tryv$mu_b   <- ini$mu_b   + scl * 0.02 * rnorm(1)
      tryv$mu_b_T <- ini$mu_b_T + scl * 0.02 * rnorm(1)
      tryv$beta   <- abs(scl * 0.2 * rnorm(1))
      tryv$u_rho  <- 0.5 + scl * 0.1 * rnorm(1)
      for (nm in intersect(names(tryv), m$getVarNames())) cm[[nm]] <- tryv[[nm]]
      if (is.finite(cm$calculate())) { ini <- tryv; break }
    }
  }
  for (nm in intersect(names(ini), m$getVarNames())) cm[[nm]] <- ini[[nm]]
  say("chain %d: start logProb %.2f\n", ch, cm$calculate())
  it <- 0; chains <- NULL; subj <- NULL
  while (it < cfg$n_iter) {
    n <- min(cfg$chunk, cfg$n_iter - it)
    cmc$run(n, reset = (it == 0), progressBar = FALSE); it <- it + n
    chains <- as.matrix(cmc$mvSamples); subj <- as.matrix(cmc$mvSamples2)
    saveRDS(chains, file.path(out_dir, sprintf("%s_chain%d_chains.rds", tier, ch)))
    saveRDS(subj,   file.path(out_dir, sprintf("%s_chain%d_subj.rds", tier, ch)))
    say("chain %d: %d/%d  (%.1f min)\n", ch, it, cfg$n_iter, el())
  }
  list(chains = chains, subj = subj)
}

## ---- the entry point -------------------------------------------------
##  n_iter and chunk override the tier's defaults. Chains run IN PARALLEL,
##  one R process each (n_cores; default min(n_chains, cores - 1)); each
##  process compiles its own copy of the model (about a minute) and writes
##  <tier>_chain<k>_chains.rds / _subj.rds after every chunk, so a run that
##  is killed keeps what every chain has done. The combined files
##  <tier>_chains.rds, <tier>_subj.rds, <tier>_fit.rds are written at the end.
##  n_cores = 1 runs the chains one after another in this process.
fit_joint_model <- function(tier = c("smoke", "struct", "infer"),
                            n_chains = 3L, seed = 100L,
                            out_dir = "results_joint", verbose = TRUE,
                            jitter = 0.05, n_iter = NULL, chunk = NULL,
                            dat = NULL, n_cores = NULL) {
  ## dat: a dataset in the format of data/adni_joint.rds (e.g. a simulated
  ## one from simulation/sim_joint.R); the real data when NULL.
  tier <- match.arg(tier)
  cfg <- switch(tier,
    smoke  = list(n_sub = 250L, n_iter =   300L, chunk =  150L),
    struct = list(n_sub =  NA_integer_, n_iter =  800L, chunk = 400L),
    infer  = list(n_sub =  NA_integer_, n_iter = 20000L, chunk = 1000L))
  if (!is.null(n_iter)) cfg$n_iter <- as.integer(n_iter)
  if (!is.null(chunk))  cfg$chunk  <- as.integer(chunk)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time(); el <- function() as.numeric(difftime(Sys.time(), t0, units = "mins"))
  say <- function(...) if (verbose) { cat(sprintf(...)); flush.console() }

  if (is.null(dat)) dat <- load_amyloid_data()
  if (!is.na(cfg$n_sub) && cfg$n_sub < dat$n_subj) {
    set.seed(7); k <- sort(sample(dat$n_subj, cfg$n_sub))
    dat <- subset_subjects(dat, k)
  }
  scans <- function(y, n) unlist(lapply(seq_along(n), function(i) y[i, seq_len(n[i])]))
  bA <- make_rate_basis(SETTINGS)
  vA <- make_var_basis(scans(dat$y_mat, dat$n_obs), SETTINGS)
  bT <- make_rate_basis(SETTINGS_TAU)
  vT <- make_var_basis(scans(dat$y_tau_mat, dat$n_obs_tau), SETTINGS_TAU)
  say("tier %s: %d subjects, %d amyloid + %d tau scans (%d subjects with tau), %d chains x %d iterations\n",
      tier, dat$n_subj, sum(dat$n_obs), sum(dat$n_obs_tau), dat$n_subj_tau,
      n_chains, cfg$n_iter)

  spec <- build_joint_model(dat, bA, vA, bT, vT)
  job  <- list(spec = spec, dat = dat, bA = bA, vA = vA, bT = bT, vT = vT,
               cfg = cfg, seed = seed, jitter = jitter, out_dir = out_dir, tier = tier)

  if (is.null(n_cores)) n_cores <- min(n_chains, max(1L, parallel::detectCores() - 1L))
  n_cores <- min(n_cores, n_chains)
  if (n_cores > 1L && n_chains > 1L) {
    ## one R process per chain, each compiling its own copy of the model and
    ## writing its own checkpoint files; the parent prints each worker's log
    say(paste0("running %d chains on %d workers. This console stays SILENT until all\n",
               "chains finish; each worker logs to %s/%s_chain<k>.log (compile takes\n",
               "several minutes per worker, then one line per chunk). Watch with\n",
               "    tail -f \"%s/%s_chain1.log\"\n"),
        n_chains, n_cores, out_dir, tier, normalizePath(out_dir), tier)
    cl <- parallel::makeCluster(n_cores, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    wd <- getwd()
    parallel::clusterCall(cl, function(wd) { setwd(wd); suppressMessages(source("run_joint_model.R")); NULL }, wd)
    res <- parallel::parLapply(cl, seq_len(n_chains), run_joint_chain, job = job, verbose = FALSE)
  } else {
    res <- lapply(seq_len(n_chains), run_joint_chain, job = job, verbose = verbose)
  }
  chains <- lapply(res, `[[`, "chains"); subj <- lapply(res, `[[`, "subj")
  saveRDS(chains, file.path(out_dir, sprintf("%s_chains.rds", tier)))
  saveRDS(subj,   file.path(out_dir, sprintf("%s_subj.rds", tier)))
  fit <- list(tier = tier, chains = chains, subj = subj, dat = dat,
              bA = bA, vA = vA, bT = bT, vT = vT,
              settings = list(A = SETTINGS, T = SETTINGS_TAU, J = SETTINGS_JOINT),
              minutes = el())
  saveRDS(fit, file.path(out_dir, sprintf("%s_fit.rds", tier)))
  say("done in %.1f min -> %s\n", el(), file.path(out_dir, sprintf("%s_fit.rds", tier)))
  invisible(fit)
}

## ---- posterior state probabilities ---------------------------------
##
##  r_{s,i} = P(S_i = s | data, latent draw), averaged over draws. Computed
##  by re-evaluating the three branches of dsubjJoint in R at each saved
##  subject draw, so it uses exactly the likelihood the sampler used.
joint_state_probs <- function(fit, max_draws = 200L) {
  dat <- fit$dat; S <- fit$settings
  ch  <- do.call(rbind, fit$chains); sb <- do.call(rbind, fit$subj)
  ## subject draws are thinned 20x relative to the population draws
  idx <- round(seq(1, nrow(sb), length.out = min(max_draws, nrow(sb))))
  pop <- ch[pmin(idx * 20, nrow(ch)), , drop = FALSE]
  get <- function(M, v, k) M[, sprintf("%s[%d]", v, k)]
  Xg  <- cbind(1, dat$X)
  R <- matrix(0, dat$n_subj, 3, dimnames = list(NULL, c("00", "10", "11")))
  for (d in seq_along(idx)) {
    th <- sapply(seq_len(fit$bA$K), function(k) get(pop[d, , drop = FALSE], "theta", k))
    thT <- sapply(seq_len(fit$bT$K), function(k) get(pop[d, , drop = FALSE], "theta_T", k))
    nu  <- sapply(seq_len(fit$vA$K), function(k) get(pop[d, , drop = FALSE], "nu", k))
    nuT <- sapply(seq_len(fit$vT$K), function(k) get(pop[d, , drop = FALSE], "nu_T", k))
    rA <- rate_grid(fit$bA, th); GA <- clock_grid(fit$bA, rA, S$A$thresh)
    rT <- rate_grid(fit$bT, thT); GT <- clock_grid(fit$bT, rT, S$T$thresh)
    beta <- pop[d, "beta"]
    mA <- exp(beta * (plogis((fit$bA$y - S$A$thresh) / S$J$d_coup) - 0.5))
    EA <- sa_cumtrapz(mA / pmax(rA, S$A$r_min), fit$bA$step)
    omA <- as.numeric(exp(fit$vA$B %*% nu)); omT <- as.numeric(exp(fit$vT$B %*% nuT))
    mu_b <- pop[d, "mu_b"]; mu_bT <- pop[d, "mu_b_T"]
    phi  <- sapply(1:3, function(k) get(pop[d, , drop = FALSE], "phi", k))
    phiT <- sapply(1:3, function(k) get(pop[d, , drop = FALSE], "phi_T", k))
    psi  <- sapply(1:2, function(k) get(pop[d, , drop = FALSE], "psi", k))
    psiT <- sapply(1:2, function(k) get(pop[d, , drop = FALSE], "psi_T", k))
    sdl <- pop[d, "sigma_delta"]; sdlT <- pop[d, "sigma_delta_T"]
    for (i in seq_len(dat$n_subj)) {
      s <- sb[idx[d], ]
      xA <- s[sprintf("x_tilde[%d]", i)]; zA <- s[sprintf("z[%d]", i)]; bAi <- s[sprintf("b[%d]", i)]
      xT <- s[sprintf("x_tilde_T[%d]", i)]; zT <- s[sprintf("z_T[%d]", i)]; bTi <- s[sprintf("b_T[%d]", i)]
      dA <- sum(dat$X[i, ] * psi) + sdl * zA; dT <- sum(dat$X[i, ] * psiT) + sdlT * zT
      path <- joint_paths(dat, i, xA, dA, xT, dT, mu_b, mu_bT,
                          fit$bA, rA, GA, EA, fit$bT, rT, GT, S, beta)
      ya <- dat$y_mat[i, seq_len(dat$n_obs[i])]; aa <- path$A
      lA1 <- sum(dt_log(ya, aa, sa_interp(fit$vA$y, omA, aa), S$A$t_df))
      lA0 <- sum(dt_log(ya, bAi, sa_interp(fit$vA$y, omA, bAi), S$A$t_df))
      lT0 <- 0; lT1 <- 0
      if (dat$n_obs_tau[i] > 0) {
        yt <- dat$y_tau_mat[i, seq_len(dat$n_obs_tau[i])]; tt <- path$Tau
        lT1 <- sum(dt_log(yt, tt, sa_interp(fit$vT$y, omT, tt), S$T$t_df))
        lT0 <- sum(dt_log(yt, bTi, sa_interp(fit$vT$y, omT, bTi), S$T$t_df))
      }
      lpi <- plogis(sum(Xg[i, ] * phi), log.p = TRUE); l1pi <- plogis(-sum(Xg[i, ] * phi), log.p = TRUE)
      lq  <- plogis(sum(Xg[i, ] * phiT), log.p = TRUE); l1q <- plogis(-sum(Xg[i, ] * phiT), log.p = TRUE)
      l <- c(l1pi + lA0 + lT0, lpi + l1q + lA1 + lT0, lpi + lq + lA1 + lT1)
      l[!is.finite(l)] <- -Inf
      R[i, ] <- R[i, ] + exp(l - max(l)) / sum(exp(l - max(l)))
    }
  }
  R / length(idx)
}

dt_log <- function(y, mu, om, df) dt((y - mu) / om, df, log = TRUE) - log(om)

## Latent paths for one subject at their own scan ages, in plain R, by the
## same rule as the compiled densities (nfExposure / nfInvExposure).
joint_paths <- function(dat, i, xA, dA, xT, dT, mu_b, mu_bT,
                        bA, rA, GA, EA, bT, rT, GT, S, beta) {
  cpl <- function(y) exp(beta * (plogis((y - S$A$thresh) / S$J$d_coup) - 0.5))
  Gb <- sa_interp(bA$y, GA, mu_b); rb <- max(sa_interp(bA$y, rA, mu_b), S$A$r_min)
  g0 <- if (xA >= mu_b) sa_interp(bA$y, GA, xA) else Gb + (xA - mu_b) / rb
  alphaA <- dat$a_tilde[i] - exp(-dA) * (g0 - Gb)
  Eb <- sa_interp(bA$y, EA, mu_b)
  level_at <- function(u) {
    gt <- pmax(Gb, g0 + exp(dA) * (u - dat$a_tilde[i]))
    sa_inv_monotone(bA$y, GA, gt)
  }
  expo <- function(u) ifelse(u <= alphaA, cpl(mu_b) * (u - alphaA),
                             exp(-dA) * (sa_interp(bA$y, EA, level_at(u)) - Eb))
  A <- level_at(dat$age_mat[i, seq_len(dat$n_obs[i])])
  Tau <- numeric(0)
  if (dat$n_obs_tau[i] > 0) {
    GbT <- sa_interp(bT$y, GT, mu_bT); rbT <- max(sa_interp(bT$y, rT, mu_bT), S$T$r_min)
    g0T <- if (xT >= mu_bT) sa_interp(bT$y, GT, xT) else GbT + (xT - mu_bT) / rbT
    at <- dat$age_tau_mat[i, seq_len(dat$n_obs_tau[i])]
    gt <- pmax(GbT, g0T + exp(dT) * (expo(at) - expo(dat$a_tilde_tau[i])))
    Tau <- sa_inv_monotone(bT$y, GT, gt)
  }
  list(A = A, Tau = Tau, alphaA = alphaA)
}
