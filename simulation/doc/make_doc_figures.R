## Figures for simulation_guide.pdf. Run from the project root:
##   Rscript simulation/doc/make_doc_figures.R
## Figure 1 uses the real data-generating code (simulate_cohort) with the
## "fitted" truth. Figure 2 is a toy illustration, not model output.
## Loads the model's R/ files directly, not run_model.R, so the figures do not
## depend on anything else run_model.R does when sourced.
for (f in sort(list.files("R", full.names = TRUE))) suppressMessages(source(f))
load_amyloid_data <- function(path = "data/adni_amyloid.rds") readRDS(path)
for (f in c("simulation/sim_truth.R", "simulation/sim_generate.R")) source(f)
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
out <- "simulation/doc"
th <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 9.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey30"),
        panel.grid.minor = element_blank())

## ---- Figure 1: what a simulated cohort looks like ---------------------
bA <- make_rate_basis(); vA <- sim_var_basis(); dat <- sim_design(NA)
truth <- sim_truth("fitted")
sim <- simulate_cohort(dat, truth, bA, vA, seed = 10001)
r <- rate_grid(bA, truth$theta); G <- clock_grid(bA, r, SETTINGS$thresh)
Gb <- sa_interp(bA$y, G, truth$mu_b)
life <- function(i, ages) {
  if (!sim$subj$susc[i]) return(rep(sim$subj$b[i], length(ages)))
  gt <- pmin(pmax(Gb + exp(sim$subj$delta[i]) * (ages - sim$subj$alpha[i]), Gb), max(G))
  sa_inv_monotone(bA$y, G, gt)
}
st <- sim$subj$state; nobs <- dat$n_obs
pick <- function(s, cond = TRUE) { set.seed(3); w <- which(st == s & nobs >= 4 & cond); w[sample.int(length(w), 2)] }
ids <- c(pick("C"), pick("L", sim$subj$alpha < 95), pick("P", sim$subj$alpha > 45))
lab <- c(C = "C: never accumulates", L = "L: will accumulate, has not started", P = "P: accumulating")
ages <- seq(40, 100, by = 0.25)
curves <- do.call(rbind, lapply(ids, function(i)
  data.frame(id = i, state = unname(lab[st[i]]), age = ages, cl = to_centiloid(life(i, ages)))))
scans <- do.call(rbind, lapply(ids, function(i) { no <- nobs[i]
  data.frame(id = i, state = unname(lab[st[i]]), age = dat$age_mat[i, 1:no], cl = to_centiloid(sim$dat$y_mat[i, 1:no])) }))
wins <- do.call(rbind, lapply(ids, function(i) { no <- nobs[i]
  data.frame(id = i, state = unname(lab[st[i]]), lo = min(dat$age_mat[i, 1:no]), hi = max(dat$age_mat[i, 1:no])) }))
alph <- data.frame(id = ids, state = unname(lab[st[ids]]), a = ifelse(sim$subj$susc[ids], sim$subj$alpha[ids], NA))
cols <- setNames(c("#2C6FA8", "#B03030"), as.character(ids[c(1, 2)]))
pal <- c("#2C6FA8", "#C07020")
curves$who <- factor(ave(curves$id, curves$state, FUN = function(v) match(v, unique(v))))
scans$who  <- factor(ave(scans$id,  scans$state,  FUN = function(v) match(v, unique(v))))
wins$who   <- factor(ave(wins$id,   wins$state,   FUN = function(v) match(v, unique(v))))
alph$who   <- factor(ave(alph$id,   alph$state,   FUN = function(v) match(v, unique(v))))
for (d in c("curves", "scans", "wins", "alph")) assign(d, within(get(d), state <- factor(state, levels = lab)))
p1 <- ggplot() +
  geom_rect(data = wins, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = who), alpha = .10) +
  geom_hline(yintercept = to_centiloid(SETTINGS$thresh), linetype = "dotted") +
  geom_line(data = curves, aes(age, cl, colour = who), linewidth = .6) +
  geom_point(data = scans, aes(age, cl, colour = who), size = 1.3) +
  geom_point(data = alph[!is.na(alph$a), ], aes(x = a, y = to_centiloid(truth$mu_b), colour = who),
             shape = 4, size = 2.2, stroke = .9) +
  facet_wrap(~ state, nrow = 1) +
  scale_colour_manual(values = pal, guide = "none") + scale_fill_manual(values = pal, guide = "none") +
  coord_cartesian(xlim = c(40, 100), ylim = c(-30, 150)) +
  labs(title = "Six subjects from one simulated cohort, two per hidden state",
       subtitle = paste("Line = the TRUE lifelong amyloid curve the simulation drew.  Dots = the simulated scans, at the subject's real ADNI visit ages.",
                        "Shaded band = the follow-up window.  X = the true departure age.  Dotted line = positivity threshold (24 Centiloid).", sep = "\n"),
       x = "age (years)", y = "amyloid (Centiloid)") + th
ggsave(file.path(out, "fig_doc1_cohort.pdf"), p1, width = 9, height = 3.4)

## ---- Figure 2: what bias and coverage mean (toy) ----------------------
## seed chosen so the "good" case misses once or twice, as a 95% interval should
set.seed(5)
truth0 <- 0.73; R <- 30
mk <- function(name, shift, se) {
  est <- truth0 + shift + rnorm(R, 0, se)
  data.frame(case = name, rep = 1:R, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
}
toy <- rbind(mk("Good: unbiased, honest intervals", 0, 0.03),
             mk("Biased: intervals centred in the wrong place", 0.05, 0.03),
             mk("Overconfident: intervals too narrow", 0, 0.03))
toy$lo[toy$case == "Overconfident: intervals too narrow"] <- with(toy[toy$case == "Overconfident: intervals too narrow", ], est - 1.96 * 0.012)
toy$hi[toy$case == "Overconfident: intervals too narrow"] <- with(toy[toy$case == "Overconfident: intervals too narrow", ], est + 1.96 * 0.012)
toy$hit <- toy$lo <= truth0 & truth0 <= toy$hi
toy$case <- factor(toy$case, levels = unique(toy$case))
cov_lab <- aggregate(hit ~ case, toy, function(v) sprintf("covers the truth %d of %d times", sum(v), length(v)))
p2 <- ggplot(toy, aes(rep, est, colour = hit)) +
  geom_hline(yintercept = truth0, linetype = "dashed", colour = "black") +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .45) + geom_point(size = .9) +
  geom_text(data = cov_lab, aes(x = 15.5, y = 0.612, label = hit), inherit.aes = FALSE, size = 2.7) +
  facet_wrap(~ case, nrow = 1) + coord_cartesian(ylim = c(0.60, 0.90)) +
  scale_colour_manual(values = c(`TRUE` = "grey35", `FALSE` = "#B03030"), guide = "none") +
  labs(title = "Toy illustration: 30 replicates of one estimand, three ways an estimator can behave",
       subtitle = "Dashed line = the truth (here 0.73, the APOE4 time ratio). Each vertical bar is one replicate's 95% interval; red bars miss the truth.\nNumbers are invented to illustrate the ideas; they are not output of the model.",
       x = "replicate", y = "estimate") + th
ggsave(file.path(out, "fig_doc2_bias_coverage.pdf"), p2, width = 9, height = 3.2)

cat(sprintf("fig 1 subjects: %s | states %s | alphas %s\n", paste(ids, collapse = ","),
            paste(st[ids], collapse = ","), paste(round(sim$subj$alpha[ids], 1), collapse = ",")))
print(table(toy$case, toy$hit))
