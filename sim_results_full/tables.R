## LaTeX table rows for the report, generated from the analysis outputs.
OUT <- "sim_results_full"
P <- read.csv(file.path(OUT, "performance_detailed.csv")); K <- read.csv(file.path(OUT, "prior_shrinkage.csv"))
S <- read.csv(file.path(OUT, "subject_level_per_rep.csv"))
f <- function(x, d) formatC(x, format = "f", digits = d)
lab <- c(OR_apoe4 = "OR, APOE4", OR_female = "OR, female", pi_bar = "$\\pi$ at average", TR_apoe4 = "TR, APOE4",
  TR_female = "TR, female", RR_apoe4 = "RR, APOE4", RR_female = "RR, female", mu_b = "floor $\\mu_b$",
  sigma_b = "$\\sigma_b$", sigma_alpha = "$\\sigma_\\alpha$", sigma_delta = "$\\sigma_\\delta$", transit = "transit, years",
  cure_apoe0 = "never accumulate, APOE4$-$", cure_apoe1 = "never accumulate, APOE4$+$",
  cross_apoe0 = "crossing age, APOE4$-$", cross_apoe1 = "crossing age, APOE4$+$")
dg <- c(OR_apoe4 = 2, OR_female = 2, pi_bar = 3, TR_apoe4 = 3, TR_female = 3, RR_apoe4 = 3, RR_female = 3, mu_b = 4,
  sigma_b = 4, sigma_alpha = 3, sigma_delta = 3, transit = 2, cure_apoe0 = 3, cure_apoe1 = 3, cross_apoe0 = 1, cross_apoe1 = 1)
rows <- vapply(names(lab), function(e) { r <- P[P$estimand == e, ]
  flag <- if (abs(r$bias_over_mcse) > 2) "\\textbf{" else "{"
  sprintf("%s & %s & %s & %s%s} & %s & %+.1f & %d/8 & %.2f \\\\", lab[e], f(r$truth, dg[e]), f(r$mean_est, dg[e]),
          flag, f(r$bias, dg[e]), f(r$bias_mcse, dg[e]), 100 * r$rel_bias, r$covered, r$sd_ratio) }, "")
writeLines(rows, file.path(OUT, "tab_perf.tex"))
kr <- vapply(seq_len(nrow(K)), function(i) with(K[i, ], sprintf("%s & %.2f & %.1f & %.3f & %.3f & %+.3f & %+.3f (%.3f) \\\\",
  label, prior_sd, truth_in_prior_sd, post_sd, prior_weight, predicted_bias, bias, bias_mcse)), "")
writeLines(kr, file.path(OUT, "tab_shrink.tex"))
sv <- function(v, d, pct = FALSE) { x <- S[[v]] * if (pct) 100 else 1; sprintf("%s [%s, %s]", f(mean(x), d), f(min(x), d), f(max(x), d)) }
sr <- c(sprintf("AUC, P(accumulator) vs true membership & %s \\\\", sv("auc", 3)),
  sprintf("Accuracy at a 0.5 cut-off & %s \\\\", sv("accuracy", 3)),
  sprintf("Mean P(accumulator), true state C & %s \\\\", sv("p_C", 2)),
  sprintf("Mean P(accumulator), true state L & %s \\\\", sv("p_L", 2)),
  sprintf("Mean P(accumulator), true state P & %s \\\\", sv("p_P", 2)),
  "\\midrule",
  sprintf("Departure age, state P: bias, years & %s \\\\", sv("alpha_P_bias", 1)),
  sprintf("Departure age, state P: RMSE, years & %s \\\\", sv("alpha_P_rmse", 1)),
  sprintf("Departure age, state P: coverage, \\%% & %s \\\\", sv("alpha_P_cover", 0, TRUE)),
  sprintf("Departure age, state L: bias, years & %s \\\\", sv("alpha_L_bias", 1)),
  sprintf("Departure age, state L: RMSE, years & %s \\\\", sv("alpha_L_rmse", 1)),
  sprintf("Departure age, state L: coverage, \\%% & %s \\\\", sv("alpha_L_cover", 0, TRUE)),
  sprintf("Log speed $\\delta$, state P: RMSE & %s \\\\", sv("delta_P_rmse", 2)),
  sprintf("Log speed $\\delta$, state P: coverage, \\%% & %s \\\\", sv("delta_P_cover", 0, TRUE)),
  "\\midrule",
  sprintf("corr($\\alpha_i, \\delta_i$) of posterior medians & %s \\\\", sv("cor_medians", 2)),
  sprintf("corr($\\alpha_i, \\delta_i$) within draw & %s \\\\", sv("cor_within", 2)))
writeLines(sr, file.path(OUT, "tab_subj.tex"))
cat(rows, kr, sr, sep = "\n")
