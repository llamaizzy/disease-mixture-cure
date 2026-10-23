## Figures for the short report, plus the Monte Carlo precision of each run's medians.
source("posterior_80k/00_load.R")
fit80 <- readRDS("results_80k/infer_fit.rds"); shp <- load_fitted_model()
D <- readRDS(file.path(OUT, "convergence_by_length.rds")); S <- readRDS(file.path(OUT, "summaries.rds"))
P <- S$P; Q <- S$Q
th <- theme_minimal(base_size = 9) + theme(plot.title = element_text(face = "bold", size = 9.5),
  plot.subtitle = element_text(size = 7.5, colour = "grey30"), panel.grid.minor = element_blank())

## ---- Monte Carlo SE of the median, in units of the 80k posterior SD -------
## 80k: from bulk ESS. shipped: 252 pooled draws whose chain order is not stored,
## so treat them as one ordered sequence (autocorrelation-aware) AND as independent;
## the truth lies between, and both are reported.
g80 <- pool_draws(fit80); g10 <- shp$g; pars <- P$par
full <- D[D$len == 80000, ]
mc <- data.frame(par = pars,
  mcse80_sd = sapply(pars, function(p) mcse_median(do.call(cbind, lapply(fit80$chains, function(m) m[26668:80000, p])))) / P$sd_80k,
  mcse10_seq_sd = sapply(pars, function(p) mcse_median(g10[, p])) / P$sd_80k,
  mcse10_iid_sd = sapply(pars, function(p) sqrt(pi / 2) * sd(g10[, p]) / sqrt(nrow(g10))) / P$sd_80k, row.names = NULL)
mc$shift_sd <- P$shift_sd; mc$shift_over_mcse <- abs(mc$shift_sd) / pmax(mc$mcse10_seq_sd, mc$mcse10_iid_sd)
write.csv(mc, file.path(OUT, "mcse_medians.csv"), row.names = FALSE)
cat(sprintf("MCSE of median / posterior SD -- 80k: max %.3f | shipped, sequence: median %.3f max %.3f | shipped, iid: median %.3f\n",
    max(mc$mcse80_sd), median(mc$mcse10_seq_sd), max(mc$mcse10_seq_sd), median(mc$mcse10_iid_sd)))
cat(sprintf("|shift| / shipped MCSE: median %.2f, max %.2f (%s); %d of %d below 2\n",
    median(mc$shift_over_mcse), max(mc$shift_over_mcse), mc$par[which.max(mc$shift_over_mcse)],
    sum(mc$shift_over_mcse < 2), nrow(mc)))

## ---- Figure 1: R-hat and ESS by chain length -----------------------------
agg <- do.call(rbind, lapply(split(D, D$len), function(d) data.frame(len = d$len[1],
  worst = max(d$rhat_split), median = median(d$rhat_split), min_ess = min(d$ess_bulk))))
pa <- ggplot(agg, aes(len / 1000)) +
  geom_hline(yintercept = 1.01, linetype = "dashed", colour = "firebrick") +
  geom_line(aes(y = worst, colour = "worst parameter")) + geom_point(aes(y = worst, colour = "worst parameter")) +
  geom_line(aes(y = median, colour = "median parameter")) + geom_point(aes(y = median, colour = "median parameter")) +
  scale_colour_manual(values = c(`worst parameter` = "#B03030", `median parameter` = "#2C6FA8"), name = NULL) +
  scale_x_continuous(breaks = agg$len / 1000) +
  labs(title = "A. Split R-hat as the chains get longer", subtitle = "First L iterations of each chain, first third discarded. Dashed line = 1.01.",
       x = "iterations per chain used (thousands)", y = "rank-normalised split R-hat") + th + theme(legend.position = c(.72, .8))
pb <- ggplot(agg, aes(len / 1000, min_ess)) + geom_hline(yintercept = 400, linetype = "dashed", colour = "firebrick") +
  geom_line(colour = "grey20") + geom_point(colour = "grey20") + scale_x_continuous(breaks = agg$len / 1000) +
  labs(title = "B. Smallest bulk effective sample size", subtitle = "Across all parameters. Dashed line = 400, the usual minimum.",
       x = "iterations per chain used (thousands)", y = "bulk ESS") + th
ggsave(file.path(OUT, "fig1_convergence.pdf"), pa | pb, width = 8.5, height = 3.2)

## ---- Figure 2: traces of the four slowest parameters ----------------------
slow <- c("gamma[1]", "mu_b", "phi[1]", "theta[3]")
tr <- do.call(rbind, lapply(seq_along(fit80$chains), function(c_) { m <- fit80$chains[[c_]]; k <- seq(1, nrow(m), by = 40)
  do.call(rbind, lapply(slow, function(p) data.frame(iter = k, val = m[k, p], par = p, chain = factor(c_)))) }))
tr$par <- factor(tr$par, levels = slow)
p2 <- ggplot(tr, aes(iter / 1000, val, colour = chain)) +
  annotate("rect", xmin = 0, xmax = 80 / 3, ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_line(linewidth = .2, alpha = .8) + facet_wrap(~ par, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c("#2C6FA8", "#C07020", "#3A8F5C")) +
  labs(title = "Trace plots of the four parameters with the smallest effective sample size",
       subtitle = "Every 40th draw. Grey = first third of each chain, discarded as burn-in. The chains overlap and wander over the same range.",
       x = "iteration (thousands)", y = NULL) + th + theme(legend.position = "bottom")
ggsave(file.path(OUT, "fig2_traces.pdf"), p2, width = 8.5, height = 4.6)

## ---- Figure 3: how far the shipped values sit from the 80k posterior -----
keep_q <- c("OR_apoe4", "OR_female", "pi_bar", "TR_apoe4", "TR_female", "RR_apoe4", "RR_female",
            "mu_b", "sigma_b", "sigma_alpha", "sigma_delta", "transit", "cure_apoe0", "cure_apoe1",
            "cross_apoe0", "cross_apoe1", sprintf("logr_%.1f", seq(.6, 1.5, .1)))
qq <- Q[Q$estimand %in% keep_q, ]; qq$estimand <- factor(qq$estimand, levels = rev(keep_q))
band <- median(pmax(mc$mcse10_seq_sd, mc$mcse10_iid_sd))
p3 <- ggplot(qq, aes(shift_sd, estimand)) +
  annotate("rect", xmin = -2 * band, xmax = 2 * band, ymin = -Inf, ymax = Inf, fill = "#EEF4FA") +
  geom_vline(xintercept = 0, colour = "grey40") + geom_point(size = 1.6, colour = "#2C6FA8") +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(title = "Shipped posterior medians, measured against the 80k posterior",
       subtitle = sprintf("x = (shipped median - 80k median) / 80k posterior SD. Shaded = +-2 typical Monte Carlo SEs of a\nmedian from the shipped 252 draws (+-%.2f SD): shifts inside it are indistinguishable from sampling noise.", 2 * band),
       x = "shift, in 80k posterior standard deviations", y = NULL) + th
ggsave(file.path(OUT, "fig3_shift.pdf"), p3, width = 6.5, height = 5.2)
cat("figures written\n")
