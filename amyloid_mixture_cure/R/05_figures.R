## =====================================================================
##  The report figures.
##
##      figs <- make_figures(fit)          # writes to figures/
##
##  Every figure is drawn from POSTERIOR DRAWS, not from a point estimate,
##  and that distinction is load-bearing in this model: several quantities
##  of interest (the departure age alpha_i, the speed delta_i) are
##  coordinated by the estimator, so a scatter of posterior medians shows
##  an association that does not exist within any single draw. Where a
##  figure uses medians it is a DISPLAY choice and is labelled as such.
## =====================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

## SUVR -> Centiloid, the scale published amyloid results are reported on.
## A LEVEL transforms with both terms; a DIFFERENCE of levels scales only.
to_centiloid <- function(y) y * 169.27 - 102.35

## Pool the chains and drop the first third of each as burn-in. Trimming the
## POOLED matrix instead would discard most of one chain and none of another.
pool_draws <- function(fit, burn_frac = 1/3) {
  if (isTRUE(fit$pooled)) return(fit$g)          # a shipped posterior: already done
  do.call(rbind, lapply(fit$chains, function(m)
    m[(floor(nrow(m) * burn_frac) + 1):nrow(m), , drop = FALSE]))
}

## Pool the per-subject draws alongside the global ones, keeping them ALIGNED:
## monitors2 is thinned, so the two matrices index different iterations.
pool_paired <- function(fit, thin2 = 20, burn_frac = 1/3) {
  if (isTRUE(fit$pooled)) return(list(g = fit$g, u = fit$u))
  gl <- list(); su <- list()
  for (c_ in seq_along(fit$chains)) {
    ns <- nrow(fit$subj[[c_]]); idx <- seq_len(ns) * thin2
    keep <- idx <= nrow(fit$chains[[c_]]); b0 <- floor(sum(keep) * burn_frac) + 1
    gl[[c_]] <- fit$chains[[c_]][idx[keep], , drop = FALSE][b0:sum(keep), , drop = FALSE]
    su[[c_]] <- fit$subj[[c_]][which(keep), , drop = FALSE][b0:sum(keep), , drop = FALSE]
  }
  list(g = do.call(rbind, gl), u = do.call(rbind, su))
}

.band <- function(M, x) {
  q <- apply(M, 1, quantile, c(.025, .5, .975), na.rm = TRUE)
  data.frame(x = x, lo = q[1, ], mid = q[2, ], hi = q[3, ])
}
.thm <- function() theme_minimal(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.8, colour = "grey30"),
        legend.title = element_blank(), legend.text = element_text(size = 6.5))

## ---------------------------------------------------------------------
##  FIGURE 1. The rate curve, and the survival curves with and without the
##  cure component. The point of the pair C/D is that a standard survival
##  model forces everyone to depart eventually; the cure model does not.
## ---------------------------------------------------------------------
fig_survival <- function(fit, out = "figures/fig1_survival.pdf", n_draw = 200, seed = 7) {
  set.seed(seed)
  bA <- fit$bA; S <- fit$settings; thr <- S$thresh; xc <- fit$dat$x_centre
  s <- pool_draws(fit); d <- sort(sample.int(nrow(s), min(n_draw, nrow(s))))
  ages <- seq(S$alpha_min, 100, by = 1)
  yobs <- unlist(lapply(seq_len(fit$dat$n_subj),
                        function(i) fit$dat$y_mat[i, 1:fit$dat$n_obs[i]]))
  prof <- list("APOE4 negative" = c(1, 0 - xc[1], 0),
               "APOE4 positive" = c(1, 1 - xc[1], 0))
  cols <- c("APOE4 negative" = "#2C6FA8", "APOE4 positive" = "#B03030")

  ## A -- the rate curve r_A(y). Grey = outside the range of the data, where
  ## the curve is prior-driven and should not be read.
  ys <- seq(bA$lo, bA$hi, length.out = 200)
  R  <- exp(predict(bA$B, ys) %*% t(s[d, sprintf("theta[%d]", 1:bA$K), drop = FALSE]))
  qd <- quantile(yobs, c(.02, .98))
  pA <- ggplot(.band(R, ys), aes(to_centiloid(x), mid)) +
    annotate("rect", xmin = -Inf, xmax = to_centiloid(qd[1]), ymin = -Inf, ymax = Inf, fill = "grey90") +
    annotate("rect", xmin = to_centiloid(qd[2]), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey90") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey25", alpha = .2) +
    geom_line(colour = "grey15") +
    geom_vline(xintercept = to_centiloid(median(s[, "mu_b"])), linetype = "dashed", colour = "firebrick") +
    geom_vline(xintercept = to_centiloid(thr), linetype = "dotted") +
    labs(title = "A. Rate curve r(y)",
         subtitle = "How fast amyloid accumulates at each LEVEL.\nDashed red = the floor mu_b; dotted = threshold;\ngrey = outside the range of the data.",
         x = "amyloid level y (Centiloid)", y = "SUVR per year") + .thm()

  ## the AFT survivor among susceptibles, and susceptibility itself
  npg <- ncol(fit$dat$X) + 1
  Sc <- function(k, w, a) {
    lam <- exp(sum(w * s[k, sprintf("gamma[%d]", 1:npg)])); kap <- 1 / s[k, "sigma_alpha"]
    Fm  <- 1 - exp(-((S$alpha_max - S$alpha_min) / lam)^kap)
    1 - (1 - exp(-(pmax(a - S$alpha_min, 0) / lam)^kap)) / Fm
  }
  Pi <- function(k, w) 1 / (1 + exp(-sum(w * s[k, sprintf("phi[%d]", 1:npg)])))

  ## C -- proper survival: conditional on being susceptible, falls to zero
  dC <- do.call(rbind, lapply(names(prof), function(pn)
    cbind(.band(vapply(d, function(k) Sc(k, prof[[pn]], ages), numeric(length(ages))), ages), g = pn)))
  pC <- ggplot(dC, aes(x, mid, colour = g, fill = g)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) + geom_line() +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "B. WITHOUT the cure component",
         subtitle = "Proper survival S(a) = P(not yet departed | susceptible).\nEveryone departs eventually, so it falls to zero.",
         x = "age", y = "probability not yet departed") +
    .thm() + theme(legend.position = c(.25, .25))

  ## D -- improper survival: the population curve, which plateaus at 1 - pi
  dD <- do.call(rbind, lapply(names(prof), function(pn)
    cbind(.band(vapply(d, function(k) 1 - Pi(k, prof[[pn]]) * (1 - Sc(k, prof[[pn]], ages)),
                       numeric(length(ages))), ages), g = pn)))
  pl <- sapply(names(prof), function(pn) median(1 - vapply(d, function(k) Pi(k, prof[[pn]]), 0)))
  pD <- ggplot(dD, aes(x, mid, colour = g, fill = g)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, colour = NA) + geom_line() +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    geom_hline(yintercept = pl, linetype = "dotted", colour = cols[names(pl)]) +
    annotate("text", x = 30, y = pl[1] + .04, label = sprintf("1 - pi = %.2f", pl[1]),
             size = 2.2, colour = cols[1], hjust = 0) +
    annotate("text", x = 30, y = pl[2] + .04, label = sprintf("1 - pi = %.2f", pl[2]),
             size = 2.2, colour = cols[2], hjust = 0) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "C. WITH the cure component",
         subtitle = "Improper survival 1 - pi F(a). It does NOT fall to zero:\nit plateaus at the cure fraction 1 - pi.",
         x = "age", y = "probability never yet positive") +
    .thm() + theme(legend.position = c(.25, .28))

  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  pdf(out, width = 12, height = 4.3); print(pA | pC | pD); dev.off()
  cat(sprintf("wrote %s   (cure fraction: APOE4- %.3f, APOE4+ %.3f)\n", out, pl[1], pl[2]))
  invisible(list(cure = pl))
}

## ---------------------------------------------------------------------
##  Per-subject quantities, computed WITHIN each draw. Returns the departure
##  age, the speed, the fitted path, and the posterior probability that the
##  subject is an accumulator.
## ---------------------------------------------------------------------
subject_quantities <- function(fit, n_draw = 60, seed = 2) {
  set.seed(seed)
  bA <- fit$bA; vA <- fit$vA; S <- fit$settings; dat <- fit$dat; n <- dat$n_subj
  pp <- pool_paired(fit); g <- pp$g; u <- pp$u
  d <- sort(sample.int(nrow(g), min(n_draw, nrow(g))))
  g <- g[d, , drop = FALSE]; u <- u[d, , drop = FALSE]; K <- nrow(g)
  XT <- u[, sprintf("x_tilde[%d]", 1:n), drop = FALSE]
  ZZ <- u[, sprintf("z[%d]", 1:n), drop = FALSE]
  BB <- u[, sprintf("b[%d]", 1:n), drop = FALSE]
  Xg <- cbind(1, dat$X); npg <- ncol(Xg)
  AL <- DE <- W <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) {
    th  <- g[k, sprintf("theta[%d]", 1:bA$K)]
    rgr <- rate_grid(bA, th); G <- clock_grid(bA, rgr, S$thresh)
    om  <- as.numeric(exp(vA$B %*% g[k, sprintf("nu[%d]", 1:vA$K)]))
    mb  <- g[k, "mu_b"]; Gb <- sa_interp(bA$y, G, mb)
    rb  <- max(sa_interp(bA$y, rgr, mb), S$r_min)
    de  <- as.numeric(dat$X %*% g[k, grep("^psi", colnames(g))]) + g[k, "sigma_delta"] * ZZ[k, ]
    DE[, k] <- de
    g0  <- ifelse(XT[k, ] >= mb, sa_interp(bA$y, G, XT[k, ]), Gb + (XT[k, ] - mb) / rb)
    AL[, k] <- dat$a_tilde - exp(-de) * (g0 - Gb)
    lpi <- as.numeric(Xg %*% g[k, sprintf("phi[%d]", 1:npg)])
    for (i in seq_len(n)) {
      no <- dat$n_obs[i]; ao <- dat$age_mat[i, 1:no]; yo <- dat$y_mat[i, 1:no]
      mu  <- sa_inv_monotone(bA$y, G, pmin(pmax(g0[i] + exp(de[i]) * (ao - dat$a_tilde[i]), Gb), max(G)))
      sij <- sa_interp(vA$y, om, pmin(pmax(mu, vA$lo), vA$hi))
      omN <- sa_interp(vA$y, om, pmin(pmax(BB[k, i], vA$lo), vA$hi))
      ## the same two log-weights the likelihood forms, turned into the
      ## posterior membership probability by one logistic step
      aS <- -log1p(exp(-lpi[i])) + sum(dt((yo - mu) / sij, S$t_df, log = TRUE) - log(sij))
      aN <- -log1p(exp( lpi[i])) + sum(dt((yo - BB[k, i]) / omN, S$t_df, log = TRUE) - log(omN))
      W[i, k] <- 1 / (1 + exp(aN - aS))
    }
  }
  list(alpha = AL, delta = DE, p_susc = W, g = g,
       alpha_med = apply(AL, 1, median), delta_med = apply(DE, 1, median),
       p_med = rowMeans(W))
}

## ---------------------------------------------------------------------
##  FIGURE 2. Every accumulator aligned to its own estimated departure age.
##  Left: the calendar axis a - alpha_i; subjects still spread apart because
##  each travels at e^{delta_i}. Right: the CLOCK axis e^{delta_i}(a - alpha_i),
##  which is what the model claims removes the remaining difference. If the
##  model holds, every accumulator collapses onto one shared curve.
##
##  alpha_i and delta_i are posterior medians here. That is a display choice:
##  the estimator coordinates these two, so the alignment inherits that and
##  must NOT be used to read off a relationship between onset and rate.
## ---------------------------------------------------------------------
fig_aligned <- function(fit, sq = NULL, out = "figures/fig2_aligned.pdf", n_draw = 60) {
  bA <- fit$bA; S <- fit$settings; dat <- fit$dat
  if (is.null(sq)) sq <- subject_quantities(fit, n_draw = n_draw)
  g <- sq$g
  keep <- which(sq$p_med >= 0.5)
  obs <- do.call(rbind, lapply(keep, function(i) { no <- dat$n_obs[i]
    data.frame(id = i, t = dat$age_mat[i, 1:no] - sq$alpha_med[i], y = dat$y_mat[i, 1:no],
               tc = exp(sq$delta_med[i]) * (dat$age_mat[i, 1:no] - sq$alpha_med[i])) }))
  thm_ <- apply(g[, sprintf("theta[%d]", 1:bA$K), drop = FALSE], 2, median)
  rgr <- rate_grid(bA, thm_); G <- clock_grid(bA, rgr, S$thresh)
  mb  <- median(g[, "mu_b"]); Gb <- sa_interp(bA$y, G, mb)

  tc  <- seq(0, quantile(obs$tc, .99), length.out = 200)
  cur <- data.frame(tc = tc, y = sa_inv_monotone(bA$y, G, pmin(Gb + tc, max(G))))
  tt  <- seq(0, quantile(obs$t, .99), length.out = 200)
  spd <- do.call(rbind, lapply(c(-1, 0, 1), function(m) {
    de <- m * median(g[, "sigma_delta"])
    data.frame(t = tt, y = sa_inv_monotone(bA$y, G, pmin(Gb + exp(de) * tt, max(G))),
               s = sprintf("delta = %+.0f sd", m)) }))

  pL <- ggplot(obs, aes(t, to_centiloid(y), group = id)) +
    geom_line(colour = "grey60", alpha = .25, linewidth = .25) +
    geom_line(data = spd, aes(t, to_centiloid(y), group = s, colour = s), linewidth = .6) +
    scale_colour_manual(values = c("#2C6FA8", "grey15", "#B03030")) +
    geom_hline(yintercept = to_centiloid(S$thresh), linetype = "dotted") +
    labs(title = "A. Aligned on the calendar axis",
         subtitle = "a - alpha_i. Subjects still spread apart: each travels at its own speed.",
         x = "years since estimated departure", y = "amyloid (Centiloid)") +
    .thm() + theme(legend.position = c(.8, .2))

  pR <- ggplot(obs, aes(tc, to_centiloid(y), group = id)) +
    geom_line(colour = "grey60", alpha = .25, linewidth = .25) +
    geom_line(data = cur, aes(tc, to_centiloid(y), group = 1), colour = "grey10", linewidth = .7) +
    geom_hline(yintercept = to_centiloid(S$thresh), linetype = "dotted") +
    labs(title = "B. Aligned on the clock axis",
         subtitle = "e^{delta_i}(a - alpha_i). One shared curve, if the model holds.",
         x = "disease time (years at unit speed)", y = "amyloid (Centiloid)") + .thm()

  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  pdf(out, width = 9, height = 4); print(pL | pR); dev.off()
  cat(sprintf("wrote %s   (%d accumulators of %d)\n", out, length(keep), dat$n_subj))
  invisible(sq)
}

## ---------------------------------------------------------------------
##  FIGURE 3. Observed against fitted, on the calendar-age axis, for the
##  accumulators. The third panel is the one that matters: a residual that
##  drifts with age would mean the shared curve is the wrong shape.
## ---------------------------------------------------------------------
fig_fit_check <- function(fit, sq = NULL, out = "figures/fig3_fit_check.pdf", n_draw = 60) {
  bA <- fit$bA; S <- fit$settings; dat <- fit$dat
  if (is.null(sq)) sq <- subject_quantities(fit, n_draw = n_draw)
  g <- sq$g
  thm_ <- apply(g[, sprintf("theta[%d]", 1:bA$K), drop = FALSE], 2, median)
  rgr <- rate_grid(bA, thm_); G <- clock_grid(bA, rgr, S$thresh)
  mb  <- median(g[, "mu_b"]); Gb <- sa_interp(bA$y, G, mb)
  rb  <- max(sa_interp(bA$y, rgr, mb), S$r_min)
  keep <- which(sq$p_med >= 0.5)

  d <- do.call(rbind, lapply(keep, function(i) {
    no <- dat$n_obs[i]; ao <- dat$age_mat[i, 1:no]
    g0 <- Gb + exp(sq$delta_med[i]) * (dat$a_tilde[i] - sq$alpha_med[i])
    gt <- pmin(pmax(g0 + exp(sq$delta_med[i]) * (ao - dat$a_tilde[i]), Gb), max(G))
    data.frame(id = i, age = ao, obs = dat$y_mat[i, 1:no],
               fit = sa_inv_monotone(bA$y, G, gt)) }))
  d$res <- to_centiloid(d$obs) - to_centiloid(d$fit)   # a DIFFERENCE: scales only

  p1 <- ggplot(d, aes(age, to_centiloid(obs), group = id)) +
    geom_line(colour = "grey55", alpha = .3, linewidth = .25) +
    geom_hline(yintercept = to_centiloid(S$thresh), linetype = "dotted") +
    labs(title = "A. Observed", x = "age", y = "amyloid (Centiloid)") + .thm()
  p2 <- ggplot(d, aes(age, to_centiloid(fit), group = id)) +
    geom_line(colour = "#2C6FA8", alpha = .3, linewidth = .25) +
    geom_hline(yintercept = to_centiloid(S$thresh), linetype = "dotted") +
    labs(title = "B. Fitted", x = "age", y = "amyloid (Centiloid)") + .thm()
  p3 <- ggplot(d, aes(age, res)) +
    geom_point(colour = "grey40", alpha = .25, size = .5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "#B03030", linewidth = .6) +
    labs(title = "C. Observed minus fitted",
         subtitle = "A drift with age would mean the shared curve is the wrong shape.",
         x = "age", y = "Centiloid") + .thm()

  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  pdf(out, width = 12, height = 4); print(p1 | p2 | p3); dev.off()
  cat(sprintf("wrote %s   (residual SD %.2f CL, mean %.3f)\n", out, sd(d$res), mean(d$res)))
  invisible(d)
}

## ---- everything at once ---------------------------------------------
make_figures <- function(fit, out_dir = "figures", n_draw = 60) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  fig_survival(fit, file.path(out_dir, "fig1_survival.pdf"))
  sq <- subject_quantities(fit, n_draw = n_draw)
  fig_aligned(fit, sq, file.path(out_dir, "fig2_aligned.pdf"))
  fig_fit_check(fit, sq, file.path(out_dir, "fig3_fit_check.pdf"))
  invisible(sq)
}

## ---------------------------------------------------------------------
##  Mean age at which the positivity threshold is crossed, by APOE4 status.
##
##      crossing age = alpha_i + T * exp(-delta_i),     T = -G_A(mu_b)
##
##  i.e. depart, then spend the transit budget at your own speed. Computed
##  WITHIN each posterior draw and averaged over the accumulators of that
##  draw -- averaging point estimates instead would mix the coordination
##  between alpha and delta into the answer.
## ---------------------------------------------------------------------
crossing_ages <- function(fit, sq = NULL, n_draw = 40) {
  bA <- fit$bA; S <- fit$settings
  if (is.null(sq)) sq <- subject_quantities(fit, n_draw = n_draw)
  g <- sq$g; K <- ncol(sq$alpha)
  Tb <- vapply(seq_len(K), function(k) {
    r <- rate_grid(bA, g[k, sprintf("theta[%d]", 1:bA$K)])
    -sa_interp(bA$y, clock_grid(bA, r, S$thresh), g[k, "mu_b"]) }, 0)
  CR   <- sq$alpha + sweep(exp(-sq$delta), 2, Tb, "*")
  sus  <- sq$p_susc >= 0.5
  apoe <- fit$dat$X[, 1] + fit$dat$x_centre[1]
  grp  <- function(sel) {
    v <- vapply(seq_len(K), function(k) mean(CR[which(sel & sus[, k]), k]), 0)
    c(median = median(v), lo = unname(quantile(v, .025)), hi = unname(quantile(v, .975))) }
  out <- rbind(`APOE4 negative` = grp(apoe == 0), `APOE4 positive` = grp(apoe == 1))
  cat(sprintf("mean threshold-crossing age\n  APOE4-  %.1f [%.1f, %.1f]\n  APOE4+  %.1f [%.1f, %.1f]\n  gap     %.1f years\n  transit budget %.1f years\n",
      out[1,1], out[1,2], out[1,3], out[2,1], out[2,2], out[2,3],
      out[1,1] - out[2,1], median(Tb)))
  invisible(list(ages = out, transit = median(Tb)))
}
