## Analysis of the completed "full" simulation tier (fitted scenario,
## 8 replicates, 3 chains x 20,000 iterations, 1,101 subjects).
## Run from the project root:  Rscript sim_results_full/analyse.R
source("run_simulation.R")
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
OUT <- "sim_results_full"
reps <- load_replicates("full"); R <- length(reps)
reps <- reps[order(sapply(reps, `[[`, "rep"))]
bA <- make_rate_basis(); truth <- sim_truth("fitted")
xc <- load_amyloid_data()$x_centre
th <- theme_minimal(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.5, colour = "grey30"), panel.grid.minor = element_blank())
cat(sprintf("%d replicates; settings: %s\n", R, paste(unique(sapply(reps, function(r) paste(r$cfg$n_chains, "x", r$cfg$n_iter))), collapse = ",")))

## ---- 1. per replicate, per estimand: estimate, interval, posterior SD, z, rank
G <- do.call(rbind, lapply(reps, function(r) {
  ed <- estimand_draws(r$draws, bA)
  Pi <- function(w) as.numeric(plogis(r$draws[, sprintf("phi[%d]", 1:3)] %*% w))
  ed <- cbind(ed, cure_apoe0 = 1 - Pi(c(1, -xc[1], 0)), cure_apoe1 = 1 - Pi(c(1, 1 - xc[1], 0)))
  g <- r$global
  g$post_sd <- ifelse(g$estimand %in% colnames(ed), apply(ed, 2, sd)[g$estimand], (g$hi - g$lo) / 3.92)
  g$rank <- ifelse(g$estimand %in% colnames(ed),
                   sapply(seq_len(nrow(g)), function(i) if (g$estimand[i] %in% colnames(ed)) mean(ed[, g$estimand[i]] < g$truth[i]) else NA), NA)
  g$rep <- r$rep; g$worst_rhat <- max(r$rhat, na.rm = TRUE); g }))
G$err <- G$est - G$truth; G$z <- G$err / G$post_sd; G$cover <- G$lo <= G$truth & G$truth <= G$hi
write.csv(G, file.path(OUT, "per_replicate_estimands.csv"), row.names = FALSE)

head_est <- c("OR_apoe4", "OR_female", "pi_bar", "TR_apoe4", "TR_female", "RR_apoe4", "RR_female",
              "mu_b", "sigma_b", "sigma_alpha", "sigma_delta", "transit", "cure_apoe0", "cure_apoe1",
              "cross_apoe0", "cross_apoe1")
curve_est <- sprintf("logr_%.1f", seq(.6, 1.5, .1))
perf <- do.call(rbind, lapply(split(G, G$estimand), function(d) data.frame(
  estimand = d$estimand[1], truth = mean(d$truth), mean_est = mean(d$est), bias = mean(d$err),
  bias_mcse = sd(d$err) / sqrt(nrow(d)), bias_over_mcse = mean(d$err) / (sd(d$err) / sqrt(nrow(d))),
  rel_bias = mean(d$err) / mean(d$truth), emp_se = sd(d$est), mean_post_sd = mean(d$post_sd),
  sd_ratio = mean(d$post_sd) / sd(d$est), rmse = sqrt(mean(d$err^2)), covered = sum(d$cover), n = nrow(d),
  mean_z = mean(d$z), sd_z = sd(d$z))))
perf <- perf[match(c(head_est, "gamma1", curve_est), perf$estimand), ]
write.csv(perf, file.path(OUT, "performance_detailed.csv"), row.names = FALSE)
op <- options(width = 200)
cat("\n=== performance (8 replicates)\n"); print(perf, digits = 3, row.names = FALSE)

## ---- 2. pooled coverage and calibration ---------------------------------
H <- G[G$estimand %in% head_est, ]; C <- G[G$estimand %in% curve_est, ]
pc <- function(d) c(intervals = nrow(d), covered = sum(d$cover), rate = mean(d$cover),
                    p_at_most = pbinom(sum(d$cover), nrow(d), .95))
cat("\n=== pooled coverage\n"); print(rbind(headline = pc(H), rate_curve = pc(C), all = pc(G)), digits = 3)
cat("\nmisses:\n"); print(G[!G$cover, c("estimand", "rep", "truth", "est", "lo", "hi", "worst_rhat")], digits = 4, row.names = FALSE)
zs <- H$z[is.finite(H$z)]
cat(sprintf("\nheadline z-scores: mean %.2f, SD %.2f (calibrated: 0 and 1)\n", mean(zs), sd(zs)))
zs2 <- H$z[is.finite(H$z) & !H$estimand %in% c("OR_apoe4", "cure_apoe0", "cure_apoe1", "OR_female")]
cat(sprintf("headline z excluding the susceptibility contrasts: mean %.2f, SD %.2f\n", mean(zs2), sd(zs2)))

## ---- 3. prior shrinkage of the covariate contrasts, on the log scale ------
## Normal-normal approximation: with prior N(0, s0^2) and posterior SD s, the
## prior's weight is w = s^2 / s0^2, and the posterior median is pulled toward
## 0 by about w * truth.  Compare with the observed bias.
con <- data.frame(par = c("phi[2]", "phi[3]", "gamma[2]", "gamma[3]", "psi[1]", "psi[2]"),
                  label = c("OR APOE4", "OR female", "TR APOE4", "TR female", "RR APOE4", "RR female"),
                  prior_sd = c(.5, .5, .25, .25, .25, .25),
                  truth = c(truth$phi[2:3], truth$gamma[2:3], truth$psi))
con <- do.call(rbind, lapply(seq_len(nrow(con)), function(k) {
  p <- con$par[k]
  est <- sapply(reps, function(r) median(r$draws[, p])); s <- sapply(reps, function(r) sd(r$draws[, p]))
  w <- mean(s)^2 / con$prior_sd[k]^2
  cbind(con[k, ], mean_est = mean(est), bias = mean(est) - con$truth[k], bias_mcse = sd(est) / sqrt(R),
        post_sd = mean(s), prior_weight = w, predicted_bias = -w * con$truth[k],
        truth_in_prior_sd = con$truth[k] / con$prior_sd[k]) }))
write.csv(con, file.path(OUT, "prior_shrinkage.csv"), row.names = FALSE)
cat("\n=== shrinkage of contrasts (log scale)\n"); print(con, digits = 3, row.names = FALSE)

## ---- 4. subject level and convergence ------------------------------------
S <- subject_table(reps)$per_rep; S$rep <- sapply(reps, `[[`, "rep")
write.csv(S, file.path(OUT, "subject_level_per_rep.csv"), row.names = FALSE)
cat("\n=== subject level, mean [min, max] over replicates\n")
for (v in c("auc", "accuracy", "p_C", "p_L", "p_P", "alpha_P_bias", "alpha_P_rmse", "alpha_P_cover",
            "alpha_L_bias", "alpha_L_rmse", "alpha_L_cover", "delta_P_bias", "delta_P_rmse", "delta_P_cover",
            "cor_medians", "cor_within", "worst_rhat"))
  cat(sprintf("  %-14s %7.3f [%7.3f, %7.3f]\n", v, mean(S[[v]]), min(S[[v]]), max(S[[v]])))
real <- sapply(reps, function(r) r$realised)
cat("\nrealised cohort composition, mean over replicates:\n"); print(round(rowMeans(real), 3))
rh <- do.call(rbind, lapply(reps, function(r) data.frame(rep = r$rep, par = names(r$rhat), rhat = r$rhat)))
cat("\nR-hat above 1.02:\n"); print(rh[rh$rhat > 1.02 & !is.na(rh$rhat), ], row.names = FALSE, digits = 4)
saveRDS(list(G = G, perf = perf, con = con, S = S, rh = rh, real = real), file.path(OUT, "analysis.rds"))

## ---- figures ---------------------------------------------------------------
## A. bias with 95% Monte Carlo interval, in units of the empirical SE
ph <- perf[perf$estimand %in% c(head_est, curve_est), ]
ph$estimand <- factor(ph$estimand, levels = rev(c(head_est, curve_est)))
ph$flag <- abs(ph$bias_over_mcse) > 2
f1 <- ggplot(ph, aes(bias / emp_se, estimand, colour = flag)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_errorbar(aes(xmin = (bias - 1.96 * bias_mcse) / emp_se, xmax = (bias + 1.96 * bias_mcse) / emp_se), width = 0, linewidth = .5, orientation = "y") +
  geom_point(size = 1.6) + scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#B03030"), guide = "none") +
  labs(title = "A. Bias, in units of the estimate's spread across replicates",
       subtitle = "Bars = 95% Monte Carlo interval from 8 replicates. Red = bias more than 2 Monte Carlo SEs from 0.",
       x = "bias / empirical SE", y = NULL) + th
cv <- aggregate(cover ~ estimand, G[G$estimand %in% c(head_est, curve_est), ], sum)
cv$estimand <- factor(cv$estimand, levels = rev(c(head_est, curve_est)))
f2 <- ggplot(cv, aes(cover, estimand)) + geom_vline(xintercept = 7.6, linetype = "dotted") +
  geom_segment(aes(x = 0, xend = cover, yend = estimand), colour = "grey70") + geom_point(size = 1.6, colour = "#2C6FA8") +
  scale_x_continuous(breaks = 0:8, limits = c(0, 8)) +
  labs(title = "B. Intervals covering the truth", subtitle = "Out of 8. Dotted = 7.6, the 95% expectation.", x = "replicates covered", y = NULL) +
  th + theme(axis.text.y = element_blank())
ggsave(file.path(OUT, "fig1_bias_coverage.pdf"), f1 + f2 + plot_layout(widths = c(2.2, 1)), width = 9, height = 6.2)

## B. the susceptibility odds ratio: replicate intervals, and predicted vs observed shrinkage
or <- G[G$estimand == "OR_apoe4", ]
g1 <- ggplot(or, aes(factor(rep), est)) + geom_hline(yintercept = or$truth[1], linetype = "dashed", colour = "firebrick") +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey40") + geom_point(colour = "#2C6FA8", size = 1.8) +
  labs(title = "A. APOE4 odds ratio, every replicate", subtitle = "Dashed = truth, 5.32. Point = posterior median.",
       x = "replicate", y = "odds ratio") + th
g2 <- ggplot(con, aes(predicted_bias, bias, label = label)) + geom_abline(linetype = "dotted") +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_errorbar(aes(ymin = bias - 1.96 * bias_mcse, ymax = bias + 1.96 * bias_mcse), width = 0, colour = "grey50") +
  geom_point(colour = "#B03030", size = 1.8) +
  ## label positions set by hand (order: OR APOE4, OR female, TR APOE4, TR female, RR APOE4, RR female)
  geom_segment(aes(xend = predicted_bias + c(.02, .02, .035, .035, .035, .035) - .003,
                   yend = c(bias[1], bias[2], -.045, -.020, .045, .020)), colour = "grey65", linewidth = .25) +
  geom_text(aes(x = predicted_bias + c(.02, .02, .035, .035, .035, .035),
                y = c(bias[1], bias[2], -.045, -.020, .045, .020)), size = 2.6, hjust = 0) +
  coord_cartesian(xlim = c(-.35, .1)) +
  labs(title = "B. Pull toward zero: predicted from the prior vs observed",
       subtitle = "Log scale. Dotted = perfect agreement. Bars = 95% Monte Carlo interval.",
       x = "predicted bias from prior weight", y = "observed bias") + th
ggsave(file.path(OUT, "fig2_shrinkage.pdf"), g1 | g2, width = 9, height = 3.6)

## C. calibration: z-scores and SD ratio
zz <- H[is.finite(H$z), ]; zz$group <- ifelse(zz$estimand %in% c("OR_apoe4", "OR_female", "cure_apoe0", "cure_apoe1"),
                                              "susceptibility contrasts", "everything else")
c1 <- ggplot(zz, aes(z, fill = group)) + geom_histogram(binwidth = .4, boundary = 0, colour = "white") +
  stat_function(fun = function(x) dnorm(x) * nrow(zz) * .4, inherit.aes = FALSE, linetype = "dashed") +
  scale_fill_manual(values = c(`everything else` = "grey60", `susceptibility contrasts` = "#B03030"), name = NULL) +
  labs(title = "A. Standardised errors, headline quantities", subtitle = "(estimate - truth) / posterior SD. Dashed = N(0, 1), what perfect calibration gives.",
       x = "z", y = "count") + th + theme(legend.position = c(.2, .85))
sr <- perf[perf$estimand %in% head_est, ]; sr$estimand <- factor(sr$estimand, levels = rev(head_est))
c2 <- ggplot(sr, aes(sd_ratio, estimand)) + geom_vline(xintercept = 1, linetype = "dotted") +
  geom_point(colour = "#2C6FA8", size = 1.6) + scale_x_log10() +
  labs(title = "B. Posterior SD / spread of estimates", subtitle = "Near 1 = intervals sized right. Note: 8 replicates make each ratio noisy.",
       x = "ratio (log scale)", y = NULL) + th
ggsave(file.path(OUT, "fig3_calibration.pdf"), c1 | c2, width = 9, height = 3.8)
cat("\nfigures written\n")
