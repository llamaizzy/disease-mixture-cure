## LaTeX table rows for the report, generated from the saved results so no
## number is transcribed by hand.
OUT <- "posterior_80k"
Q <- read.csv(file.path(OUT, "estimands_80k_vs_shipped.csv")); T <- read.csv(file.path(OUT, "sim_truth_shipped_vs_80k.csv"))
C <- read.csv(file.path(OUT, "convergence_80k.csv"))
f <- function(x, d) formatC(x, format = "f", digits = d)
ci <- function(m, l, h, d) sprintf("%s [%s, %s]", f(m, d), f(l, d), f(h, d))
row <- function(lab, e, d, mult = 1) { r <- Q[Q$estimand == e, ]
  sprintf("%s & %s & %s \\\\", lab, ci(mult * r$med_80k, mult * r$lo_80k, mult * r$hi_80k, d),
          ci(mult * r$med_shipped, mult * r$lo_shipped, mult * r$hi_shipped, d)) }
res <- c(
  "\\multicolumn{3}{@{}l}{\\emph{Whether: odds ratio}} \\\\",
  row("APOE4", "OR_apoe4", 2), row("female", "OR_female", 2), row("$\\pi$ at cohort average", "pi_bar", 3),
  "\\multicolumn{3}{@{}l}{\\emph{When: time ratio, below 1 is earlier}} \\\\",
  row("APOE4", "TR_apoe4", 3), row("female", "TR_female", 3),
  "\\multicolumn{3}{@{}l}{\\emph{How fast: rate ratio}} \\\\",
  row("APOE4", "RR_apoe4", 3), row("female", "RR_female", 3),
  "\\multicolumn{3}{@{}l}{\\emph{Shared structure and derived quantities}} \\\\",
  row("floor $\\mu_b$, SUVR", "mu_b", 3), row("$\\sigma_b$", "sigma_b", 4), row("$\\sigma_\\alpha$", "sigma_alpha", 3),
  row("$\\sigma_\\delta$", "sigma_delta", 3), row("transit budget, years", "transit", 1),
  row("never accumulate, APOE4$-$, \\%", "cure_apoe0", 1, 100), row("never accumulate, APOE4$+$, \\%", "cure_apoe1", 1, 100),
  row("crossing age, APOE4$-$", "cross_apoe0", 1), row("crossing age, APOE4$+$", "cross_apoe1", 1))
writeLines(res, file.path(OUT, "tab_results.tex"))

lab <- c(OR_apoe4 = "OR, APOE4", OR_female = "OR, female", pi_bar = "$\\pi$ at average", TR_apoe4 = "TR, APOE4",
         TR_female = "TR, female", RR_apoe4 = "RR, APOE4", RR_female = "RR, female", mu_b = "$\\mu_b$",
         sigma_alpha = "$\\sigma_\\alpha$", sigma_delta = "$\\sigma_\\delta$", transit = "transit, years",
         cure_apoe0 = "never accumulate, APOE4$-$", cure_apoe1 = "never accumulate, APOE4$+$", logr_1.2 = "$\\log r$ at 1.2 SUVR")
dg <- c(OR_apoe4 = 2, OR_female = 2, pi_bar = 3, TR_apoe4 = 3, TR_female = 3, RR_apoe4 = 3, RR_female = 3, mu_b = 4,
        sigma_alpha = 3, sigma_delta = 3, transit = 2, cure_apoe0 = 3, cure_apoe1 = 3, logr_1.2 = 3)
tt <- vapply(names(lab), function(e) { r <- T[T$estimand == e, ]
  sprintf("%s & %s & %s & %+.2f \\\\", lab[e], f(r$truth_shipped, dg[e]), f(r$truth_80k, dg[e]), r$change_sd) }, "")
writeLines(tt, file.path(OUT, "tab_truth.tex"))
cat(res, sep = "\n"); cat(tt, sep = "\n")
cat(sprintf("\nconvergence: worst split %.4f (%s), worst classic %.4f, min bulk ESS %.0f (%s), min tail ESS %.0f (%s)\n",
    max(C$rhat_split), C$par[which.max(C$rhat_split)], max(C$rhat_classic), min(C$ess_bulk), C$par[which.min(C$ess_bulk)],
    min(C$ess_tail), C$par[which.min(C$ess_tail)]))
cat(sprintf("max |change_sd| in sim truth: %.2f (%s)\n", max(abs(T$change_sd)), T$estimand[which.max(abs(T$change_sd))]))
