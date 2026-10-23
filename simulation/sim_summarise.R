## =====================================================================
##  Scoring. Two levels:
##
##    score_replicate()  one fit against its truth -> a small list saved
##                       to disk (never the chains: 200 replicates of
##                       chains would be gigabytes and answer nothing)
##    summarise_sim()    all saved replicates of a tier -> the performance
##                       tables and the figures
##
##  Performance measures follow Morris, White & Crowther (2019): bias,
##  empirical SE, RMSE, coverage of the 95% interval, mean interval width,
##  each with its Monte Carlo standard error, so that "coverage 0.92 from 40
##  replicates" is read with its +-0.04 attached.
## =====================================================================

## Gelman-Rubin R-hat per column, the same rule as convergence() but quiet.
rhat_quiet <- function(chains, burn_frac = 1/3) {
  L <- lapply(chains, function(m) m[(floor(nrow(m) * burn_frac) + 1):nrow(m), , drop = FALSE])
  n <- nrow(L[[1]])
  M <- sapply(L, colMeans); V <- sapply(L, function(x) apply(x, 2, var))
  B <- n * apply(M, 1, var); W <- rowMeans(V)
  r <- sqrt(((n - 1) / n * W + B / n) / W); r[!is.finite(r)] <- NA
  r
}

## Area under the ROC curve of a score against a binary truth (rank form).
auc <- function(score, truth) {
  if (length(unique(truth)) < 2) return(NA_real_)
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## Per-draw estimand vectors from the pooled global draws; the same
## transformations as sim_estimands() so truth and estimate line up.
estimand_draws <- function(s, bA, S = SETTINGS, y_eval = seq(0.6, 1.5, by = 0.1)) {
  TH <- s[, sprintf("theta[%d]", 1:bA$K), drop = FALSE]
  Gb <- vapply(seq_len(nrow(s)), function(k) {
    r <- rate_grid(bA, TH[k, ]); sa_interp(bA$y, clock_grid(bA, r, S$thresh), s[k, "mu_b"]) }, 0)
  LR <- t(vapply(seq_len(nrow(s)), function(k) sa_interp(bA$y, as.numeric(bA$B %*% TH[k, ]), y_eval),
                 numeric(length(y_eval))))
  colnames(LR) <- sprintf("logr_%.1f", y_eval)
  cbind(OR_apoe4 = exp(s[, "phi[2]"]), OR_female = exp(s[, "phi[3]"]), pi_bar = plogis(s[, "phi[1]"]),
        TR_apoe4 = exp(s[, "gamma[2]"]), TR_female = exp(s[, "gamma[3]"]),
        RR_apoe4 = exp(s[, "psi[1]"]), RR_female = exp(s[, "psi[2]"]),
        mu_b = s[, "mu_b"], sigma_b = s[, "sigma_b"], sigma_alpha = s[, "sigma_alpha"],
        sigma_delta = s[, "sigma_delta"], gamma1 = s[, "gamma[1]"], transit = -Gb, LR)
}

## Score one replicate. `sim` is simulate_cohort()'s output, `truth` the
## scenario truth.
score_replicate <- function(fit, sim, truth, n_draw_subj = 100) {
  bA <- fit$bA; S <- fit$settings; dat <- fit$dat; xc <- dat$x_centre
  s  <- pool_draws(fit)
  ed <- estimand_draws(s, bA)
  ## the cure fraction at the two APOE4 profiles, as fig_survival() defines it
  Pi <- function(w) as.numeric(plogis(s[, sprintf("phi[%d]", 1:3)] %*% w))   # a vector, or cbind drops the name
  ed <- cbind(ed, cure_apoe0 = 1 - Pi(c(1, 0 - xc[1], 0)), cure_apoe1 = 1 - Pi(c(1, 1 - xc[1], 0)))
  q  <- apply(ed, 2, quantile, c(.5, .025, .975))

  ## subject level, WITHIN draw, using the package's own routine
  sq <- subject_quantities(fit, n_draw = n_draw_subj)
  cr <- capture.output(ca <- crossing_ages(fit, sq))
  tr_sub <- sim$subj; acc <- tr_sub$susc
  qi <- function(M) t(apply(M, 1, quantile, c(.5, .025, .975)))
  Aq <- qi(sq$alpha); Dq <- qi(sq$delta)
  by_state <- function(err, cov, st) c(
    bias = mean(err[st]), mae = mean(abs(err[st])), rmse = sqrt(mean(err[st]^2)),
    cover = mean(cov[st]), n = sum(st))
  a_err <- Aq[, 1] - tr_sub$alpha; a_cov <- Aq[, 2] <= tr_sub$alpha & tr_sub$alpha <= Aq[, 3]
  d_err <- Dq[, 1] - tr_sub$delta; d_cov <- Dq[, 2] <= tr_sub$delta & tr_sub$delta <= Dq[, 3]
  P <- tr_sub$state == "P"; L <- tr_sub$state == "L"
  est_acc <- sq$p_med >= 0.5
  ## the pitfall from README section 7: medians correlate, draws do not
  cm_ <- if (sum(est_acc) > 3) cor(sq$alpha_med[est_acc], sq$delta_med[est_acc]) else NA
  cw_ <- mean(vapply(seq_len(ncol(sq$alpha)), function(k) {
    a <- sq$p_susc[, k] >= 0.5; if (sum(a) > 3) cor(sq$alpha[a, k], sq$delta[a, k]) else NA }, 0), na.rm = TRUE)

  truth_vec <- c(sim_estimands(truth, bA), sim$realised[c("cure_apoe0", "cure_apoe1")])
  truth_vec <- truth_vec[colnames(ed)]
  glob <- data.frame(estimand = colnames(ed), truth = truth_vec,
                     est = q[1, ], lo = q[2, ], hi = q[3, ], row.names = NULL)
  cross <- data.frame(estimand = c("cross_apoe0", "cross_apoe1"),
                      truth = unname(sim$realised[c("cross_apoe0", "cross_apoe1")]),
                      est = ca$ages[, "median"], lo = ca$ages[, "lo"], hi = ca$ages[, "hi"], row.names = NULL)
  rh <- rhat_quiet(fit$chains)
  list(
    global   = rbind(glob, cross),
    rhat     = rh,
    subject  = list(
      auc_susc   = auc(sq$p_med, acc),
      acc_susc   = mean(est_acc == acc),
      p_by_state = tapply(sq$p_med, tr_sub$state, mean),
      alpha_P = by_state(a_err, a_cov, P), alpha_L = by_state(a_err, a_cov, L),
      delta_P = by_state(d_err, d_cov, P), delta_L = by_state(d_err, d_cov, L),
      cor_medians = cm_, cor_within = cw_),
    realised = sim$realised,
    draws    = s[seq(1, nrow(s), length.out = min(400, nrow(s))), , drop = FALSE],
    minutes  = fit$minutes, n_iter = nrow(fit$chains[[1]]), n_chains = length(fit$chains))
}

## ---- aggregate a tier -------------------------------------------------
load_replicates <- function(tier, dir = "results/sim") {
  fs <- list.files(file.path(dir, tier), pattern = "^rep_.*\\.rds$", full.names = TRUE, recursive = TRUE)
  out <- lapply(fs, readRDS)
  names(out) <- basename(fs)
  out[!vapply(out, function(r) isTRUE(r$error), TRUE)]
}

performance_table <- function(reps) {
  g <- do.call(rbind, lapply(seq_along(reps), function(k)
    cbind(rep = reps[[k]]$rep, scenario = reps[[k]]$scenario, reps[[k]]$global)))
  g$err <- g$est - g$truth; g$cov <- g$lo <= g$truth & g$truth <= g$hi; g$width <- g$hi - g$lo
  agg <- function(d) {
    R <- nrow(d); cov <- mean(d$cov)
    data.frame(n_rep = R, truth = mean(d$truth), mean_est = mean(d$est),
               bias = mean(d$err), bias_mcse = sd(d$err) / sqrt(R),
               rel_bias = mean(d$err) / mean(d$truth),
               emp_se = sd(d$est), rmse = sqrt(mean(d$err^2)),
               coverage = cov, cov_mcse = sqrt(cov * (1 - cov) / R),
               mean_width = mean(d$width))
  }
  sp <- split(g, list(g$scenario, g$estimand), drop = TRUE)
  out <- do.call(rbind, lapply(sp, agg))
  out <- cbind(scenario = sub("\\..*$", "", rownames(out)),
               estimand = sub("^[^.]*\\.", "", rownames(out)), out)
  rownames(out) <- NULL
  out[order(out$scenario, match(out$estimand, unique(g$estimand))), ]
}

subject_table <- function(reps) {
  rows <- lapply(reps, function(r) { s <- r$subject
    data.frame(scenario = r$scenario, auc = s$auc_susc, accuracy = s$acc_susc,
               p_C = s$p_by_state["C"], p_L = s$p_by_state["L"], p_P = s$p_by_state["P"],
               alpha_P_bias = s$alpha_P["bias"], alpha_P_rmse = s$alpha_P["rmse"], alpha_P_cover = s$alpha_P["cover"],
               alpha_L_bias = s$alpha_L["bias"], alpha_L_rmse = s$alpha_L["rmse"], alpha_L_cover = s$alpha_L["cover"],
               delta_P_bias = s$delta_P["bias"], delta_P_rmse = s$delta_P["rmse"], delta_P_cover = s$delta_P["cover"],
               cor_medians = s$cor_medians, cor_within = s$cor_within,
               worst_rhat = max(r$rhat, na.rm = TRUE), n_above_1.05 = sum(r$rhat > 1.05, na.rm = TRUE),
               minutes = r$minutes, row.names = NULL) })
  d <- do.call(rbind, rows)
  agg <- aggregate(d[, -1], list(scenario = d$scenario), function(v) mean(v, na.rm = TRUE))
  list(per_rep = d, mean = agg)
}

## ---- figures ----------------------------------------------------------
summarise_sim <- function(tier, dir = "results/sim", out_dir = file.path("figures", "sim", tier)) {
  suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
  gs <- function(...) suppressWarnings(ggsave(...))     # one replicate makes empty MC-error bars
  reps <- load_replicates(tier, dir)
  if (!length(reps)) stop("no finished replicates for tier '", tier, "'")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  perf <- performance_table(reps); subj <- subject_table(reps)
  write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)
  write.csv(subj$mean, file.path(out_dir, "subject_level.csv"), row.names = FALSE)
  write.csv(subj$per_rep, file.path(out_dir, "subject_level_per_rep.csv"), row.names = FALSE)

  head_est <- c("OR_apoe4", "OR_female", "pi_bar", "TR_apoe4", "TR_female", "RR_apoe4", "RR_female",
                "mu_b", "sigma_b", "sigma_alpha", "sigma_delta", "transit",
                "cure_apoe0", "cure_apoe1", "cross_apoe0", "cross_apoe1")
  g <- do.call(rbind, lapply(seq_along(reps), function(k)
    cbind(rep = reps[[k]]$rep, scenario = reps[[k]]$scenario, reps[[k]]$global)))
  g$estimand <- factor(g$estimand, levels = unique(g$estimand))
  gh <- g[g$estimand %in% head_est, ]
  th <- theme_minimal(base_size = 8.5) + theme(plot.title = element_text(face = "bold", size = 9),
                                               plot.subtitle = element_text(size = 6.8, colour = "grey30"))

  ## 1. every replicate's posterior median and 95% interval against the truth
  gh$cov <- gh$lo <= gh$truth & gh$truth <= gh$hi
  tru <- unique(gh[, c("scenario", "estimand", "truth")])
  p1 <- ggplot(gh, aes(x = rep, y = est, colour = cov)) +
    geom_hline(data = tru, aes(yintercept = truth), linetype = "dashed", colour = "firebrick",
               inherit.aes = FALSE) +
    geom_linerange(aes(ymin = lo, ymax = hi), alpha = .5, linewidth = .3) + geom_point(size = .6) +
    scale_colour_manual(values = c(`TRUE` = "grey30", `FALSE` = "#B03030"), guide = "none") +
    facet_grid(estimand ~ scenario, scales = "free_y") +
    labs(title = "Posterior medians and 95% intervals, every replicate",
         subtitle = "Dashed red = truth. Red intervals miss it.", x = "replicate", y = NULL) + th
  gs(file.path(out_dir, "fig_sim1_intervals.pdf"), p1, width = 3 + 2.4 * length(unique(gh$scenario)),
         height = 1.1 * length(head_est) + 1, limitsize = FALSE)

  ## 2. coverage and standardised bias with Monte Carlo error
  ph <- perf[perf$estimand %in% head_est, ]
  ph$estimand <- factor(ph$estimand, levels = head_est)
  p2a <- ggplot(ph, aes(x = estimand, y = coverage, colour = scenario)) +
    geom_hline(yintercept = .95, linetype = "dotted") +
    geom_pointrange(aes(ymin = pmax(coverage - 1.96 * cov_mcse, 0), ymax = pmin(coverage + 1.96 * cov_mcse, 1)),
                    position = position_dodge(width = .6), size = .25) +
    coord_cartesian(ylim = c(0.5, 1)) + labs(title = "A. Coverage of the 95% interval", x = NULL, y = "coverage") +
    th + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  p2b <- ggplot(ph, aes(x = estimand, y = bias / emp_se, colour = scenario)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_pointrange(aes(ymin = (bias - 1.96 * bias_mcse) / emp_se, ymax = (bias + 1.96 * bias_mcse) / emp_se),
                    position = position_dodge(width = .6), size = .25) +
    labs(title = "B. Bias in units of the empirical SE", x = NULL, y = "bias / emp. SE") +
    th + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  gs(file.path(out_dir, "fig_sim2_coverage_bias.pdf"), p2a / p2b, width = 9, height = 6.5)

  ## 3. the rate curve, pointwise: truth, the spread of the estimates, coverage
  gr <- g[grepl("^logr_", g$estimand), ]
  gr$y <- as.numeric(sub("logr_", "", gr$estimand)); gr$cov <- gr$lo <= gr$truth & gr$truth <= gr$hi
  cr <- do.call(rbind, lapply(split(gr, list(gr$scenario, gr$y), drop = TRUE), function(d)
    data.frame(scenario = d$scenario[1], y = d$y[1], truth = d$truth[1], med = median(d$est),
               lo = quantile(d$est, .025), hi = quantile(d$est, .975), coverage = mean(d$cov))))
  p3a <- ggplot(cr, aes(x = to_centiloid(y))) +
    geom_ribbon(aes(ymin = exp(lo), ymax = exp(hi), fill = scenario), alpha = .2) +
    geom_line(aes(y = exp(med), colour = scenario)) +
    geom_line(aes(y = exp(truth)), colour = "black", linetype = "dashed") +
    labs(title = "A. Rate curve r(y): truth (dashed) and the 2.5-97.5% range of posterior medians",
         x = "amyloid level (Centiloid)", y = "SUVR per year") + th
  p3b <- ggplot(cr, aes(x = to_centiloid(y), y = coverage, colour = scenario)) +
    geom_hline(yintercept = .95, linetype = "dotted") + geom_line() + geom_point(size = .8) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(title = "B. Pointwise coverage of log r(y)", x = "amyloid level (Centiloid)", y = "coverage") + th
  gs(file.path(out_dir, "fig_sim3_rate_curve.pdf"), p3a / p3b, width = 8, height = 6.5)

  ## 4. subject level: membership recovery and the alpha/delta pitfall
  d <- subj$per_rep
  pl <- rbind(data.frame(scenario = d$scenario, state = "C (cured)", p = d$p_C),
              data.frame(scenario = d$scenario, state = "L (latent)", p = d$p_L),
              data.frame(scenario = d$scenario, state = "P (progressing)", p = d$p_P))
  p4a <- ggplot(pl, aes(x = state, y = p, colour = scenario)) + geom_boxplot(outlier.size = .5) +
    labs(title = "A. Mean posterior P(susceptible), by TRUE state",
         subtitle = "C should be low, P high; L is genuinely ambiguous and should sit between.", x = NULL, y = "P(susceptible)") + th
  p4b <- ggplot(d, aes(x = scenario, y = auc)) + geom_boxplot(outlier.size = .5) +
    labs(title = "B. AUC: P(susceptible) against true membership", x = NULL, y = "AUC") + th
  pc <- rbind(data.frame(scenario = d$scenario, kind = "corr of posterior medians", r = d$cor_medians),
              data.frame(scenario = d$scenario, kind = "corr within draw", r = d$cor_within))
  p4c <- ggplot(pc, aes(x = scenario, y = r, colour = kind)) + geom_hline(yintercept = 0, linetype = "dotted") +
    geom_boxplot(outlier.size = .5) +
    labs(title = "C. corr(alpha_i, delta_i) among estimated accumulators",
         subtitle = "Truth: alpha and delta are independent. Medians correlate anyway; draws should not.", x = NULL, y = "correlation") + th
  pe <- rbind(data.frame(scenario = d$scenario, state = "P", rmse = d$alpha_P_rmse, cover = d$alpha_P_cover),
              data.frame(scenario = d$scenario, state = "L", rmse = d$alpha_L_rmse, cover = d$alpha_L_cover))
  p4d <- ggplot(pe, aes(x = state, y = cover, colour = scenario)) + geom_hline(yintercept = .95, linetype = "dotted") +
    geom_boxplot(outlier.size = .5) + coord_cartesian(ylim = c(0, 1)) +
    labs(title = "D. Coverage of the subject-level 95% interval for the departure age", x = "true state", y = "coverage") + th
  gs(file.path(out_dir, "fig_sim4_subjects.pdf"), (p4a | p4b) / (p4c | p4d), width = 11, height = 7)

  cat(sprintf("\n== tier '%s': %d replicates scored; tables and 4 figures -> %s/\n", tier, length(reps), out_dir))
  print_sim_summary(perf, subj$mean)
  invisible(list(performance = perf, subject = subj))
}

print_sim_summary <- function(perf, subj_mean) {
  show <- c("OR_apoe4", "TR_apoe4", "RR_apoe4", "pi_bar", "mu_b", "sigma_alpha", "sigma_delta",
            "transit", "cure_apoe0", "cure_apoe1", "cross_apoe0", "cross_apoe1")
  for (sc in unique(perf$scenario)) {
    p <- perf[perf$scenario == sc & perf$estimand %in% show, ]
    cat(sprintf("\n--- scenario %s (%d replicates)\n", sc, p$n_rep[1]))
    cat(sprintf("  %-12s %8s %8s %8s %8s %9s\n", "estimand", "truth", "mean", "rel.bias", "RMSE", "coverage"))
    for (i in seq_len(nrow(p))) cat(sprintf("  %-12s %8.3f %8.3f %7.1f%% %8.3f %5.2f +-%.2f\n",
      p$estimand[i], p$truth[i], p$mean_est[i], 100 * p$rel_bias[i], p$rmse[i], p$coverage[i], p$cov_mcse[i]))
    s <- subj_mean[subj_mean$scenario == sc, ]
    cat(sprintf("  subjects: AUC %.3f, accuracy %.3f | P(susc) by true state C %.2f  L %.2f  P %.2f\n",
                s$auc, s$accuracy, s$p_C, s$p_L, s$p_P))
    cat(sprintf("  alpha (state P): bias %+.2f y, RMSE %.2f y, coverage %.2f | corr medians %+.2f vs within-draw %+.2f\n",
                s$alpha_P_bias, s$alpha_P_rmse, s$alpha_P_cover, s$cor_medians, s$cor_within))
    cat(sprintf("  mean worst R-hat %.3f, mean %.1f params above 1.05, %.1f min per replicate\n",
                s$worst_rhat, s$n_above_1.05, s$minutes))
  }
}
