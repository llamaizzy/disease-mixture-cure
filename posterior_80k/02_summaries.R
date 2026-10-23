## Posterior inference on the 3 x 80,000 run, and comparison with the shipped
## posterior (3 x 10,000, 252 thinned draws) that the simulation study uses as truth.
source("posterior_80k/00_load.R")
for (f in c("simulation/sim_truth.R", "simulation/sim_summarise.R")) source(f)
fit80 <- readRDS("results_80k/infer_fit.rds")
shp   <- load_fitted_model()
g80 <- pool_draws(fit80)          # first third of each chain discarded -> 160,002 draws
g10 <- pool_draws(shp)            # already pooled: 252 draws
bA <- fit80$bA; xc <- fit80$dat$x_centre
cat(sprintf("draws: 80k run %d, shipped %d\n", nrow(g80), nrow(g10)))

## ---- 1. every sampled parameter: medians, intervals, shift in posterior-SD units
pars <- setdiff(intersect(colnames(g80), colnames(g10)), "theta[1]")
P <- data.frame(par = pars,
  med_80k = apply(g80[, pars], 2, median), sd_80k = apply(g80[, pars], 2, sd),
  lo_80k = apply(g80[, pars], 2, quantile, .025), hi_80k = apply(g80[, pars], 2, quantile, .975),
  med_shipped = apply(g10[, pars], 2, median),
  lo_shipped = apply(g10[, pars], 2, quantile, .025), hi_shipped = apply(g10[, pars], 2, quantile, .975),
  row.names = NULL)
P$shift_sd <- (P$med_shipped - P$med_80k) / P$sd_80k
P$width_ratio <- (P$hi_shipped - P$lo_shipped) / (P$hi_80k - P$lo_80k)
P$inside <- P$med_shipped >= P$lo_80k & P$med_shipped <= P$hi_80k
write.csv(P, file.path(OUT, "parameters_80k_vs_shipped.csv"), row.names = FALSE)

## ---- 2. reported quantities, posterior summaries computed per draw
est_draws <- function(g) {
  set.seed(1); k <- if (nrow(g) > 20000) sort(sample.int(nrow(g), 20000)) else seq_len(nrow(g))
  e <- estimand_draws(g[k, ], bA)
  Pi <- function(w) as.numeric(plogis(g[k, sprintf("phi[%d]", 1:3)] %*% w))
  cbind(e, cure_apoe0 = 1 - Pi(c(1, 0 - xc[1], 0)), cure_apoe1 = 1 - Pi(c(1, 1 - xc[1], 0)))
}
E80 <- est_draws(g80); E10 <- est_draws(g10)
qs <- function(E) t(apply(E, 2, quantile, c(.5, .025, .975)))
q80 <- qs(E80); q10 <- qs(E10)
Q <- data.frame(estimand = rownames(q80), med_80k = q80[, 1], lo_80k = q80[, 2], hi_80k = q80[, 3],
                med_shipped = q10[, 1], lo_shipped = q10[, 2], hi_shipped = q10[, 3],
                sd_80k = apply(E80, 2, sd), row.names = NULL)
Q$shift_sd <- (Q$med_shipped - Q$med_80k) / Q$sd_80k

## crossing ages need the per-subject draws, computed within draw by the package routine
cat("subject quantities (this takes a few minutes) ...\n")
sq80 <- subject_quantities(fit80, n_draw = 200, seed = 2)
sq10 <- subject_quantities(shp,   n_draw = 200, seed = 2)
co80 <- capture.output(ca80 <- crossing_ages(fit80, sq80))
co10 <- capture.output(ca10 <- crossing_ages(shp, sq10))
add <- function(nm, a80, a10, sdv) data.frame(estimand = nm, med_80k = a80[1], lo_80k = a80[2], hi_80k = a80[3],
  med_shipped = a10[1], lo_shipped = a10[2], hi_shipped = a10[3], sd_80k = sdv, shift_sd = (a10[1] - a80[1]) / sdv, row.names = NULL)
sd_cross <- function(sq, fit) {     # per-draw mean crossing ages, to get a posterior SD
  K <- ncol(sq$alpha); g <- sq$g
  Tb <- vapply(seq_len(K), function(k) -sa_interp(bA$y, clock_grid(bA, rate_grid(bA, g[k, sprintf("theta[%d]", 1:bA$K)]), SETTINGS$thresh), g[k, "mu_b"]), 0)
  CR <- sq$alpha + sweep(exp(-sq$delta), 2, Tb, "*"); sus <- sq$p_susc >= .5
  apoe <- fit$dat$X[, 1] + fit$dat$x_centre[1]
  sapply(c(0, 1), function(a) sd(vapply(seq_len(K), function(k) mean(CR[which(apoe == a & sus[, k]), k]), 0)))
}
sdc <- sd_cross(sq80, fit80)
Q <- rbind(Q, add("cross_apoe0", ca80$ages[1, ], ca10$ages[1, ], sdc[1]),
              add("cross_apoe1", ca80$ages[2, ], ca10$ages[2, ], sdc[2]))
## number of subjects classed as accumulators (P >= 0.5), a figure the report quotes
Q <- rbind(Q, data.frame(estimand = "n_accumulators", med_80k = sum(sq80$p_med >= .5), lo_80k = NA, hi_80k = NA,
  med_shipped = sum(sq10$p_med >= .5), lo_shipped = NA, hi_shipped = NA, sd_80k = NA, shift_sd = NA))
write.csv(Q, file.path(OUT, "estimands_80k_vs_shipped.csv"), row.names = FALSE)

## ---- 3. the simulation truth each posterior would imply (plug-in medians)
t10 <- sim_truth("fitted", post = shp$g); t80 <- sim_truth("fitted", post = g80)
T <- data.frame(estimand = names(sim_estimands(t10, bA)),
                truth_shipped = sim_estimands(t10, bA), truth_80k = sim_estimands(t80, bA), row.names = NULL)
T <- rbind(T, data.frame(estimand = c("cure_apoe0", "cure_apoe1"),
  truth_shipped = c(1 - plogis(sum(c(1, -xc[1], 0) * t10$phi)), 1 - plogis(sum(c(1, 1 - xc[1], 0) * t10$phi))),
  truth_80k     = c(1 - plogis(sum(c(1, -xc[1], 0) * t80$phi)), 1 - plogis(sum(c(1, 1 - xc[1], 0) * t80$phi)))))
T$change <- T$truth_80k - T$truth_shipped
sdmap <- setNames(Q$sd_80k, Q$estimand); T$change_sd <- T$change / sdmap[T$estimand]
write.csv(T, file.path(OUT, "sim_truth_shipped_vs_80k.csv"), row.names = FALSE)
saveRDS(list(P = P, Q = Q, T = T, t10 = t10, t80 = t80, sq80_pmed = sq80$p_med, sq10_pmed = sq10$p_med),
        file.path(OUT, "summaries.rds"))

op <- options(width = 160)
cat("\n=== sampled parameters: shipped median vs 80k posterior (shift in 80k posterior SDs)\n")
print(within(P[order(-abs(P$shift_sd)), c("par", "med_80k", "lo_80k", "hi_80k", "med_shipped", "shift_sd", "width_ratio", "inside")],
             { med_80k <- signif(med_80k, 4); lo_80k <- signif(lo_80k, 4); hi_80k <- signif(hi_80k, 4)
               med_shipped <- signif(med_shipped, 4); shift_sd <- round(shift_sd, 2); width_ratio <- round(width_ratio, 2) }), row.names = FALSE)
cat("\n=== reported quantities\n")
print(within(Q[, c("estimand", "med_80k", "lo_80k", "hi_80k", "med_shipped", "lo_shipped", "hi_shipped", "shift_sd")],
             { for (v in c("med_80k", "lo_80k", "hi_80k", "med_shipped", "lo_shipped", "hi_shipped")) assign(v, signif(get(v), 4)); shift_sd <- round(shift_sd, 2) }), row.names = FALSE)
cat("\n=== simulation truth: shipped vs what the 80k run implies\n")
print(within(T, { truth_shipped <- signif(truth_shipped, 4); truth_80k <- signif(truth_80k, 4); change <- signif(change, 3); change_sd <- round(change_sd, 2) }), row.names = FALSE)
cat("\n80k crossing ages:\n"); cat(co80, sep = "\n"); cat("\nshipped crossing ages:\n"); cat(co10, sep = "\n")
