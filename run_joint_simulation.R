## =====================================================================
##  SIMULATION STUDY FOR THE JOINT MODEL: does using tau help recover
##  AMYLOID?  This is the only file you need to call.
##
##      source("run_joint_simulation.R")
##      run_joint_simulation("pilot")   # ~6 min: 250 subjects, 1 rep, 200 iters -- wiring
##      run_joint_simulation("small")   # ~1-2 h: full design, 3 reps, 1 chain x 4000, in parallel
##      summarise_joint_sim("small")    # tables from whatever has finished
##
##  DESIGN (agreed 2026-10-23)
##   * Truth: the JOINT model at the posterior means of the production fit
##     (posterior 75k NEO modified: beta = 1.76, rho = 0.39), i.e. the
##     coupling is real. Simulated on the REAL design:
##     every subject's amyloid and tau scan ages, counts and covariates as
##     observed; only the measurements are drawn (simulation/sim_joint.R).
##   * Each simulated cohort is fitted TWICE: by the amyloid-only model
##     (R/04_model.R, which never sees the tau scans) and by the joint model
##     (R/07_joint_model.R). Both start from their own data-driven inits,
##     not from the truth.
##   * Scored on the AMYLOID side, identically for both fits:
##       subject level  posterior-mean latent amyloid level vs the true
##                      latent level, at the subject's own scan ages and at
##                      5 years past the last scan (the extrapolation the
##                      model exists for); departure-age error for true
##                      accumulators who departed before their last scan;
##                      AUC and accuracy of P(amyloid accumulator)
##       population     cure fractions by APOE4, pi_bar, mu_b, sigma_b,
##                      sigma_alpha, sigma_delta, gamma[1], covariate effects,
##                      transit budget, log rate at 10 levels -- bias, RMSE,
##                      95% coverage (Morris, White & Crowther 2019)
##     plus, for the joint fit only, recovery of beta, rho, q_bar and the
##     tau population parameters.
##   * Every replicate is checkpointed as results_joint/sim/<tier>/rep_<r>.rds
##     (scores only, never chains); a killed run resumes.
## =====================================================================
if (!file.exists("run_joint_model.R") || !dir.exists("simulation"))
  stop("run this from the amyloid_mixture_cure project folder: setwd(\"/path/to/amyloid_mixture_cure\") first")
suppressMessages(source("run_joint_model.R"))
for (f in sort(list.files("simulation", pattern = "\\.R$", full.names = TRUE))) source(f)

JOINT_SIM_TIERS <- list(
  pilot  = list(n_sub = 250L,        n_rep =   1L, n_iter =   200L, n_chains = 1L),  # ~6 min: wiring
  small  = list(n_sub = NA_integer_, n_rep =   3L, n_iter =  4000L, n_chains = 1L),  # ~1 h:   direction only
  medium = list(n_sub = NA_integer_, n_rep =  20L, n_iter =  8000L, n_chains = 1L),  # ~5 h:   bias/RMSE with MC error
  full   = list(n_sub = NA_integer_, n_rep = 100L, n_iter = 10000L, n_chains = 1L))  # ~1 day: coverage

## ---- truth -----------------------------------------------------------
## Posterior means of the NEO-T fit, with everything the joint model needs.
joint_sim_truth <- function(fit_path = "posterior 75k NEO modified/infer_fit.rds") {
  fit <- readRDS(fit_path)
  tr <- joint_truth_from_fit(fit)
  tr$source <- fit_path
  tr
}

## ---- the true latent amyloid level, for scoring -------------------------
## For accumulators the clamped clock trajectory; for non-accumulators the
## flat level b. `ages` is a list, one numeric vector per subject.
true_amyloid_levels <- function(truth, latent, ages, bA, S = SETTINGS) {
  r <- rate_grid(bA, truth$theta); G <- clock_grid(bA, r, S$thresh, S$r_min)
  Gb <- sa_interp(bA$y, G, truth$mu_b)
  lapply(seq_along(ages), function(i) {
    if (latent$ZA[i] == 1) {
      gt <- pmin(Gb + exp(latent$deltaA[i]) * pmax(ages[[i]] - latent$alphaA[i], 0), max(G))
      sa_inv_monotone(bA$y, G, gt)
    } else rep(latent$bA[i], length(ages[[i]]))
  })
}

## ---- per-draw amyloid predictions from EITHER fit -----------------------
##
##  For each of n_draw posterior draws and each subject: the accumulator
##  path A_i(a) at the requested ages, the flat level b_i, and the within-draw
##  posterior probability r_i that the subject is an amyloid accumulator.
##  The prediction scored is the mixture mean r A(a) + (1 - r) b.
##
##  r_i is where the two models differ: the amyloid-only fit forms it from
##  the two amyloid branches; the joint fit from the three states, so the
##  tau scans enter through the state weights and the tau likelihood.
amyloid_prediction_draws <- function(fit, ages, n_draw = 60, seed = 2) {
  set.seed(seed)
  joint <- !is.null(fit$bT)
  S  <- if (joint) fit$settings else list(A = fit$settings, T = NULL, J = NULL)
  SA <- S$A; dat <- fit$dat; n <- dat$n_subj; bA <- fit$bA; vA <- fit$vA
  pp <- pool_paired(fit); g <- pp$g; u <- pp$u
  d <- sort(sample.int(nrow(g), min(n_draw, nrow(g)))); g <- g[d, , drop = FALSE]; u <- u[d, , drop = FALSE]
  K <- nrow(g); Xg <- cbind(1, dat$X)
  PRED <- lapply(seq_len(n), function(i) matrix(NA_real_, K, length(ages[[i]])))
  RESP <- ALPHA <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) {
    th <- g[k, sprintf("theta[%d]", 1:bA$K)]; rA <- rate_grid(bA, th); GA <- clock_grid(bA, rA, SA$thresh, SA$r_min)
    mu_b <- g[k, "mu_b"]; Gb <- sa_interp(bA$y, GA, mu_b); rb <- max(sa_interp(bA$y, rA, mu_b), SA$r_min)
    omA <- as.numeric(exp(vA$B %*% g[k, sprintf("nu[%d]", 1:vA$K)]))
    psi <- g[k, sprintf("psi[%d]", 1:ncol(dat$X))]; sdl <- g[k, "sigma_delta"]
    lpi <- as.numeric(Xg %*% g[k, sprintf("phi[%d]", 1:3)])
    xA <- u[k, sprintf("x_tilde[%d]", 1:n)]; zA <- u[k, sprintf("z[%d]", 1:n)]; bA_i <- u[k, sprintf("b[%d]", 1:n)]
    dA <- as.numeric(dat$X %*% psi) + sdl * zA
    g0 <- ifelse(xA >= mu_b, sa_interp(bA$y, GA, xA), Gb + (xA - mu_b) / rb)
    ALPHA[, k] <- dat$a_tilde - exp(-dA) * (g0 - Gb)
    if (joint) {
      bT <- fit$bT; vT <- fit$vT; ST <- S$T
      thT <- g[k, sprintf("theta_T[%d]", 1:bT$K)]; rT <- rate_grid(bT, thT); GT <- clock_grid(bT, rT, ST$thresh, ST$r_min)
      beta <- g[k, "beta"]; EA <- sa_cumtrapz(exp(beta * (plogis((bA$y - SA$thresh) / S$J$d_coup) - 0.5)) / pmax(rA, SA$r_min), bA$step)
      omT <- as.numeric(exp(vT$B %*% g[k, sprintf("nu_T[%d]", 1:vT$K)])); mu_bT <- g[k, "mu_b_T"]
      psiT <- g[k, sprintf("psi_T[%d]", 1:ncol(dat$X))]; sdlT <- g[k, "sigma_delta_T"]
      lq <- as.numeric(Xg %*% g[k, sprintf("phi_T[%d]", 1:3)])
      xT <- u[k, sprintf("x_tilde_T[%d]", 1:n)]; zT <- u[k, sprintf("z_T[%d]", 1:n)]; bT_i <- u[k, sprintf("b_T[%d]", 1:n)]
      dT <- as.numeric(dat$X %*% psiT) + sdlT * zT
    }
    for (i in seq_len(n)) {
      level_at <- function(a) sa_inv_monotone(bA$y, GA, pmin(pmax(g0[i] + exp(dA[i]) * (a - dat$a_tilde[i]), Gb), max(GA)))
      ## amyloid branches at the scan ages
      ya <- dat$y_mat[i, seq_len(dat$n_obs[i])]; aa <- dat$age_mat[i, seq_len(dat$n_obs[i])]
      mu <- level_at(aa)
      lA1 <- sum(dt_log(ya, mu, sa_interp(vA$y, omA, pmin(pmax(mu, vA$lo), vA$hi)), SA$t_df))
      lA0 <- sum(dt_log(ya, bA_i[i], sa_interp(vA$y, omA, pmin(pmax(bA_i[i], vA$lo), vA$hi)), SA$t_df))
      if (!joint) {
        aS <- -log1p(exp(-lpi[i])) + lA1; aN <- -log1p(exp(lpi[i])) + lA0
        r <- 1 / (1 + exp(aN - aS))
      } else {
        lT0 <- 0; lT1 <- 0
        if (dat$n_obs_tau[i] > 0) {
          path <- joint_paths(dat, i, xA[i], dA[i], xT[i], dT[i], mu_b, mu_bT, bA, rA, GA, EA, bT, rT, GT, S, beta)
          yt <- dat$y_tau_mat[i, seq_len(dat$n_obs_tau[i])]; tt <- path$Tau
          lT1 <- sum(dt_log(yt, tt, sa_interp(vT$y, omT, pmin(pmax(tt, vT$lo), vT$hi)), ST$t_df))
          lT0 <- sum(dt_log(yt, bT_i[i], sa_interp(vT$y, omT, pmin(pmax(bT_i[i], vT$lo), vT$hi)), ST$t_df))
        }
        l <- c(-log1p(exp(lpi[i])) + lA0 + lT0,
               -log1p(exp(-lpi[i])) - log1p(exp(lq[i])) + lA1 + lT0,
               -log1p(exp(-lpi[i])) - log1p(exp(-lq[i])) + lA1 + lT1)
        l[!is.finite(l)] <- -Inf; w <- exp(l - max(l)); w <- w / sum(w)
        r <- w[2] + w[3]
      }
      RESP[i, k] <- r
      PRED[[i]][k, ] <- r * level_at(ages[[i]]) + (1 - r) * bA_i[i]
    }
  }
  list(pred = PRED, resp = RESP, alpha = ALPHA, g = g)
}

## ---- score one fit against the truth --------------------------------------
score_amyloid_recovery <- function(fit, sim, truth, model, n_draw = 60) {
  dat <- fit$dat; n <- dat$n_subj; bA <- fit$bA; lat <- sim$latent
  SA <- if (!is.null(fit$bT)) fit$settings$A else fit$settings
  ## evaluation ages: the scan ages, then last scan + 5 years
  ages <- lapply(seq_len(n), function(i) { a <- dat$age_mat[i, seq_len(dat$n_obs[i])]; c(a, max(a) + 5) })
  tru <- true_amyloid_levels(truth, lat, ages, bA, SA)
  pd  <- amyloid_prediction_draws(fit, ages, n_draw = n_draw)
  acc <- lat$ZA == 1
  pm  <- lapply(pd$pred, colMeans)
  lo  <- lapply(pd$pred, function(M) apply(M, 2, quantile, .05)); hi <- lapply(pd$pred, function(M) apply(M, 2, quantile, .95))
  err_scan <- unlist(lapply(seq_len(n), function(i) { m <- length(ages[[i]]); pm[[i]][-m] - tru[[i]][-m] }))
  acc_scan <- unlist(lapply(seq_len(n), function(i) rep(acc[i], length(ages[[i]]) - 1)))
  err_ext  <- vapply(seq_len(n), function(i) { m <- length(ages[[i]]); pm[[i]][m] - tru[[i]][m] }, 0)
  cov_ext  <- vapply(seq_len(n), function(i) { m <- length(ages[[i]]); lo[[i]][m] <= tru[[i]][m] && tru[[i]][m] <= hi[[i]][m] }, TRUE)
  rmse <- function(e) sqrt(mean(e^2))
  ## departure age: true accumulators who departed before their last scan
  last <- vapply(seq_len(n), function(i) max(dat$age_mat[i, seq_len(dat$n_obs[i])]), 0)
  obs_dep <- acc & lat$alphaA < last
  a_med <- apply(pd$alpha, 1, median); a_lo <- apply(pd$alpha, 1, quantile, .05); a_hi <- apply(pd$alpha, 1, quantile, .95)
  p_med <- rowMeans(pd$resp)
  subject <- c(
    rmse_scan_all = rmse(err_scan), rmse_scan_acc = rmse(err_scan[acc_scan]), rmse_scan_non = rmse(err_scan[!acc_scan]),
    rmse_ext5_all = rmse(err_ext), rmse_ext5_acc = rmse(err_ext[acc]), rmse_ext5_non = rmse(err_ext[!acc]),
    bias_ext5_acc = mean(err_ext[acc]), cover90_ext5_all = mean(cov_ext), cover90_ext5_acc = mean(cov_ext[acc]),
    alpha_bias = mean((a_med - lat$alphaA)[obs_dep]), alpha_rmse = rmse((a_med - lat$alphaA)[obs_dep]),
    alpha_cover90 = mean((a_lo <= lat$alphaA & lat$alphaA <= a_hi)[obs_dep]), n_obs_dep = sum(obs_dep),
    auc_susc = auc(p_med, acc), accuracy_susc = mean((p_med >= 0.5) == acc),
    p_susc_true_acc = mean(p_med[acc]), p_susc_true_non = mean(p_med[!acc]))
  ## population estimands, the same transformations as the amyloid study
  s <- pool_draws(fit); ed <- estimand_draws(s, bA, SA); xc <- dat$x_centre
  Pi <- function(w) as.numeric(plogis(s[, sprintf("phi[%d]", 1:3)] %*% w))
  ed <- cbind(ed, cure_apoe0 = 1 - Pi(c(1, 0 - xc[1], 0)), cure_apoe1 = 1 - Pi(c(1, 1 - xc[1], 0)))
  tv <- c(sim_estimands(truth, bA, SA), cure_apoe0 = 1 - plogis(sum(c(1, 0 - xc[1], 0) * truth$phi)),
          cure_apoe1 = 1 - plogis(sum(c(1, 1 - xc[1], 0) * truth$phi)))[colnames(ed)]
  q <- apply(ed, 2, quantile, c(.5, .025, .975))
  global <- data.frame(model = model, estimand = colnames(ed), truth = tv, est = q[1, ], lo = q[2, ], hi = q[3, ], row.names = NULL)
  if (!is.null(fit$bT)) {                       # joint-only extras
    jt <- c(beta = truth$beta, rho_alpha = truth$rho_alpha, q_bar = plogis(truth$phi_T[1]),
            mu_b_T = truth$mu_b_T, sigma_b_T = truth$sigma_b_T, gamma_T1 = truth$gamma_T[1],
            sigma_alpha_T = truth$sigma_alpha_T, sigma_delta_T = truth$sigma_delta_T, psi_T1 = truth$psi_T[1])
    col <- c(beta = "beta", rho_alpha = "rho_alpha", q_bar = "q_bar", mu_b_T = "mu_b_T", sigma_b_T = "sigma_b_T",
             gamma_T1 = "gamma_T[1]", sigma_alpha_T = "sigma_alpha_T", sigma_delta_T = "sigma_delta_T", psi_T1 = "psi_T[1]")
    qj <- sapply(col, function(v) quantile(s[, v], c(.5, .025, .975)))
    global <- rbind(global, data.frame(model = model, estimand = names(jt), truth = unname(jt),
                                       est = qj[1, ], lo = qj[2, ], hi = qj[3, ], row.names = NULL))
  }
  list(model = model, subject = subject, global = global, minutes = fit$minutes,
       n_iter = nrow(fit$chains[[1]]), n_chains = length(fit$chains))
}

## ---- one replicate: simulate, fit twice, score twice ----------------------
joint_sim_replicate <- function(r, truth, dat, bA, vA, bT, vT, cfg, seed, out_dir) {
  f <- file.path(out_dir, sprintf("rep_%02d.rds", r))
  if (file.exists(f)) return(invisible(readRDS(f)))
  t0 <- Sys.time()
  sim <- simulate_joint(dat, truth, bA, bT, vA, vT, seed = seed + r)
  res <- list(rep = r, seed = seed + r, realised = c(frac_ZA = mean(sim$latent$ZA), frac_ZT = mean(sim$latent$ZT),
                                                    off_grid = sim$off_grid))
  ## amyloid-only: the existing simulation machinery (compile once, fit)
  res$amyloid <- tryCatch({
    W <- sim_compile(dat, bA, vA)
    fitA <- sim_fit(W, sim$dat, n_chains = cfg$n_chains, n_iter = cfg$n_iter, seed = seed + 100 * r)
    score_amyloid_recovery(fitA, sim, truth, "amyloid_only")
  }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
  ## joint
  res$joint <- tryCatch({
    fitJ <- fit_joint_model("struct", n_chains = cfg$n_chains, n_iter = cfg$n_iter, chunk = cfg$n_iter,
                            seed = seed + 100 * r, dat = sim$dat, n_cores = 1L, verbose = FALSE,
                            out_dir = file.path(out_dir, sprintf("rep_%02d_jointfit", r)))
    score_amyloid_recovery(fitJ, sim, truth, "joint")
  }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
  res$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(res, f)
  invisible(res)
}

## ---- the entry point -------------------------------------------------------
run_joint_simulation <- function(tier = c("pilot", "small"), n_rep = NULL, n_iter = NULL,
                                 seed = 2026L, parallel = TRUE, dir = "results_joint/sim",
                                 truth = NULL) {
  tier <- match.arg(tier); cfg <- JOINT_SIM_TIERS[[tier]]
  if (!is.null(n_rep))  cfg$n_rep  <- as.integer(n_rep)
  if (!is.null(n_iter)) cfg$n_iter <- as.integer(n_iter)
  out_dir <- file.path(dir, tier); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(truth)) truth <- joint_sim_truth()
  saveRDS(truth, file.path(out_dir, "truth.rds"))
  dat <- load_amyloid_data()
  if (!is.na(cfg$n_sub) && cfg$n_sub < dat$n_subj) { set.seed(7); dat <- subset_subjects(dat, sort(sample(dat$n_subj, cfg$n_sub))) }
  scans <- function(y, n) unlist(lapply(seq_along(n), function(i) y[i, seq_len(n[i])]))
  bA <- make_rate_basis(SETTINGS); vA <- make_var_basis(scans(dat$y_mat, dat$n_obs), SETTINGS)
  bT <- make_rate_basis(SETTINGS_TAU); vT <- make_var_basis(scans(dat$y_tau_mat, dat$n_obs_tau), SETTINGS_TAU)
  todo <- Filter(function(r) !file.exists(file.path(out_dir, sprintf("rep_%02d.rds", r))), seq_len(cfg$n_rep))
  cat(sprintf("tier %s: %d subjects, truth beta = %.2f; %d replicates (%d to do), %d chain(s) x %d iterations each, two fits per replicate\n",
              tier, dat$n_subj, truth$beta, cfg$n_rep, length(todo), cfg$n_chains, cfg$n_iter))
  if (!length(todo)) return(invisible(summarise_joint_sim(tier, dir)))
  job <- list(truth = truth, dat = dat, bA = bA, vA = vA, bT = bT, vT = vT, cfg = cfg, seed = seed, out_dir = out_dir)
  worker <- function(r, job) joint_sim_replicate(r, job$truth, job$dat, job$bA, job$vA, job$bT, job$vT, job$cfg, job$seed, job$out_dir)
  if (parallel && length(todo) > 1) {
    nc <- min(length(todo), max(1L, parallel::detectCores() - 1L))
    cat(sprintf("running %d replicates on %d workers; this console is silent until they finish. Progress: ls %s\n", length(todo), nc, out_dir))
    cl <- parallel::makeCluster(nc, outfile = ""); on.exit(parallel::stopCluster(cl), add = TRUE)
    wd <- getwd(); parallel::clusterCall(cl, function(wd) { setwd(wd); suppressMessages(source("run_joint_simulation.R")); NULL }, wd)
    parallel::parLapply(cl, todo, worker, job = job)
  } else for (r in todo) { cat(sprintf("replicate %d ...\n", r)); worker(r, job) }
  summarise_joint_sim(tier, dir)
}

## ---- summary: amyloid-only vs joint ------------------------------------------
summarise_joint_sim <- function(tier = "small", dir = "results_joint/sim") {
  out_dir <- file.path(dir, tier)
  reps <- lapply(list.files(out_dir, pattern = "^rep_\\d+\\.rds$", full.names = TRUE), readRDS)
  if (!length(reps)) stop("no finished replicates in ", out_dir)
  ok <- function(x) !isTRUE(x$error)
  errs <- unlist(lapply(reps, function(r) c(if (!ok(r$amyloid)) sprintf("rep %d amyloid: %s", r$rep, r$amyloid$message),
                                             if (!ok(r$joint))   sprintf("rep %d joint: %s",   r$rep, r$joint$message))))
  if (length(errs)) cat("errors:\n ", paste(errs, collapse = "\n  "), "\n")
  ## subject-level metrics, mean (sd) over replicates, side by side
  S <- do.call(rbind, lapply(reps, function(r) rbind(
    if (ok(r$amyloid)) data.frame(rep = r$rep, model = "amyloid_only", t(r$amyloid$subject)),
    if (ok(r$joint))   data.frame(rep = r$rep, model = "joint", t(r$joint$subject)))))
  met <- setdiff(colnames(S), c("rep", "model"))
  agg <- function(m) sapply(c("amyloid_only", "joint"), function(mod) { v <- S[S$model == mod, m]; c(mean = mean(v), sd = sd(v), n = length(v)) })
  tab <- do.call(rbind, lapply(met, function(m) { a <- agg(m)
    data.frame(metric = m, amyloid_only = a["mean", "amyloid_only"], joint = a["mean", "joint"],
               diff_joint_minus_amy = a["mean", "joint"] - a["mean", "amyloid_only"],
               sd_amyloid = a["sd", "amyloid_only"], sd_joint = a["sd", "joint"], n_rep = min(a["n", ])) }))
  ## population estimands
  G <- do.call(rbind, lapply(reps, function(r) rbind(if (ok(r$amyloid)) cbind(rep = r$rep, r$amyloid$global),
                                                      if (ok(r$joint)) cbind(rep = r$rep, r$joint$global))))
  G$err <- G$est - G$truth; G$cov <- G$lo <= G$truth & G$truth <= G$hi
  P <- do.call(rbind, lapply(split(G, list(G$model, G$estimand), drop = TRUE), function(d)
    data.frame(model = d$model[1], estimand = d$estimand[1], truth = d$truth[1], mean_est = mean(d$est),
               bias = mean(d$err), rmse = sqrt(mean(d$err^2)), coverage = mean(d$cov), width = mean(d$hi - d$lo), n_rep = nrow(d))))
  rownames(P) <- NULL
  write.csv(tab, file.path(out_dir, "subject_comparison.csv"), row.names = FALSE)
  write.csv(P,   file.path(out_dir, "population_performance.csv"), row.names = FALSE)
  ## print
  cat(sprintf("\n== tier %s: %d replicates; truth beta = %.2f ==\n", tier, length(reps), readRDS(file.path(out_dir, "truth.rds"))$beta))
  cat("\n-- AMYLOID recovery, subject level (mean over replicates; lower RMSE / higher AUC and coverage is better)\n")
  cat(sprintf("  %-20s %12s %12s %10s\n", "metric", "amyloid-only", "joint", "joint-amy"))
  for (i in seq_len(nrow(tab))) cat(sprintf("  %-20s %12.4f %12.4f %+10.4f\n", tab$metric[i], tab$amyloid_only[i], tab$joint[i], tab$diff_joint_minus_amy[i]))
  show <- c("cure_apoe0", "cure_apoe1", "pi_bar", "OR_apoe4", "TR_apoe4", "RR_apoe4", "mu_b", "sigma_b", "sigma_alpha", "sigma_delta", "gamma1", "transit")
  cat("\n-- AMYLOID population estimands: truth | amyloid-only mean (RMSE, cover) | joint mean (RMSE, cover)\n")
  for (e in show) { a <- P[P$model == "amyloid_only" & P$estimand == e, ]; j <- P[P$model == "joint" & P$estimand == e, ]
    if (nrow(a) && nrow(j)) cat(sprintf("  %-12s %7.3f | %7.3f (%.3f, %.2f) | %7.3f (%.3f, %.2f)\n", e, a$truth, a$mean_est, a$rmse, a$coverage, j$mean_est, j$rmse, j$coverage)) }
  cat(sprintf("\n-- rate curve log r(y): RMSE over the 10 levels   amyloid-only %.3f   joint %.3f\n",
      sqrt(mean(P$rmse[P$model == "amyloid_only" & grepl("^logr_", P$estimand)]^2)),
      sqrt(mean(P$rmse[P$model == "joint" & grepl("^logr_", P$estimand)]^2))))
  jx <- P[P$model == "joint" & P$estimand %in% c("beta", "rho_alpha", "q_bar", "mu_b_T", "gamma_T1", "sigma_delta_T", "sigma_alpha_T"), ]
  if (nrow(jx)) { cat("\n-- joint-only parameters: truth | mean (RMSE, cover)\n")
    for (i in seq_len(nrow(jx))) cat(sprintf("  %-14s %7.3f | %7.3f (%.3f, %.2f)\n", jx$estimand[i], jx$truth[i], jx$mean_est[i], jx$rmse[i], jx$coverage[i])) }
  mins <- sapply(reps, function(r) r$minutes); cat(sprintf("\n%.1f min per replicate (both fits)\n", mean(mins)))
  invisible(list(subject = tab, population = P, reps = reps))
}
