# OUTLINE FOR GENERATING SIMULATION DATA

# Pipeline: truth (posterior medians) + design data (resample real data/distribution)
# -> generate latent quantities per subject (random)
# -> build true trajectories via the clock
# -> add measurement noise to get y_mat
# -> package data into same list format as adni_amyloid.rds
# -> save in rds files

# Run from the package root:
#     Rscript simulated_data/simulate_data.R     # check + plot ONE baseline dataset, then all scenarios
#     source("simulated_data/simulate_data.R")   # just load the functions

# needs SETTINGS, the integrator and the bases (not the NIMBLE model)
if (!exists("SETTINGS")) {
  for (f in file.path("R", c("00_settings.R", "01_integrator.R", "02_bases.R"))) source(f)
}

###############################################
# 1. TRUTH: posterior medians for population parameters
################################################
# Source: data/fitted_posterior.rds (use $g = 252 pooled draws x 33 params)
# Extract: phi[1:3], gamma[1:3], sigma_alpha, psi[1:2], sigma_delta, mu_b, sigma_b, theta[1:10], nu[1:8]
# Apply median to all these values
# Return in a truth list

get_baseline_truth <- function(path="data/fitted_posterior.rds") {
  post <- readRDS(path)
  med <- apply(post$g, 2, median) # one median per parameter
  get <- function(nm) unname(med[nm])
  list(
    phi = get(c("phi[1]", "phi[2]", "phi[3]")),          # susceptibility
    gamma = get(c("gamma[1]", "gamma[2]", "gamma[3]")),  # onset age (Weibull AFT)
    sigma_alpha = get("sigma_alpha"),                     # onset-age spread
    psi = get(c("psi[1]", "psi[2]")),                     # rate covariate effects
    sigma_delta = get("sigma_delta"),                     # between-person speed spread
    mu_b = get("mu_b"),                                   # shared floor
    sigma_b = get("sigma_b"),                         # non-accumulator level spread
    theta = get(sprintf("theta[%d]", 1:10)),        # rate curve coeffs
    nu = get(sprintf("nu[%d]", 1:8)),            # measurement-scale curve omega_A(y)
    N = 1101L                                    # cohort size
  )
}

########################################
# 2. SCENARIOS
########################################
# Each scenarios is baseline truth with one modification
#   baseline: no change
#   sample size: n * 0.5 (try half)
#   onset spread lo/hi: sigma_alpha * 0.5 or * 2
#   speed spread lo/hi: sigma_delta * 0.5 or * 2
#   measure error lo/hi: omega * 0.5 or * 2
# Check every scenario stays inside the prior supports in SETTINGS

make_scenarios <- function(truth) {
  mod <- function(...) modifyList(truth, list(...))
  # the var basis rows sum to 1, so omega * c is the same as nu + log(c):
  # the scaled curve is still inside the model's own parameter space
  sc <- list(
    baseline = truth,
    sample_size_half = mod(N = as.integer(round(truth$N * 0.5))),
    onset_spread_lo = mod(sigma_alpha = truth$sigma_alpha * 0.5),
    onset_spread_hi = mod(sigma_alpha = truth$sigma_alpha * 2),
    speed_spread_lo = mod(sigma_delta = truth$sigma_delta * 0.5),
    speed_spread_hi = mod(sigma_delta = truth$sigma_delta * 2),
    meas_error_lo = mod(nu = truth$nu + log(0.5)),
    meas_error_hi = mod(nu = truth$nu + log(2))
  )
  for (nm in names(sc)) check_truth_support(sc[[nm]], nm)
  sc
}

# stop if a truth sits outside the prior supports the model is fitted with
check_truth_support <- function(truth, name = "truth", S = SETTINGS) {
  ok <- c(
    sigma_alpha = truth$sigma_alpha > 0 && truth$sigma_alpha < 2,
    sigma_delta = truth$sigma_delta > 0 && truth$sigma_delta < S$sigma_delta_max,
    sigma_b     = truth$sigma_b > 0 && truth$sigma_b < 0.5,
    mu_b        = truth$mu_b > S$grid_lo && truth$mu_b < S$thresh,
    nu          = all(truth$nu > S$nu_lo & truth$nu < S$nu_hi),
    theta1      = isTRUE(all.equal(truth$theta[1], S$theta1)),
    N           = truth$N >= 2
  )
  if (!all(ok)) stop(sprintf("scenario '%s' outside prior support: %s",
                             name, paste(names(ok)[!ok], collapse = ", ")))
  invisible(TRUE)
}

#######################################
# 3. DESIGN: Covariates + visited ages
######################################
# Covariates drawn independently each subject:
#   APOE4 ~ Bernoulli(0.407)
#   Female ~ Bernoulli(0.489)
# x centre = (apoe4 = 0.427, female = 0.503)
# X = APOE4 - x centre, female - x centre (add centering)
# Baseline age ~ Normal(72.5, 7.4), truncated [51, 94] - carriers are 1.7 years younger
# num subjects = 1101
# num obs = sample from visit ranges based on probability --> mean of 3 scans
# gap between scans (n_obs - 1 gaps): normal(2.05, 0.15)
# visit ages = baseline age + cumulative sum of gaps
# write visit ages into age_mat, padding with 0s after
# CHECK: simulated follow-up length (last scan age - first) = max ~14, median ~4

# Normal truncated to [lo, hi], by inverse CDF
rtnorm <- function(n, mean, sd, lo = -Inf, hi = Inf) {
  u <- runif(n, pnorm(lo, mean, sd), pnorm(hi, mean, sd))
  qnorm(u, mean, sd)
}

get_design <- function(N = 1101, p_apoe4 = 0.407, p_female = 0.489,
                       x_centre = c(apoe4 = 0.427, female = 0.503),
                       age_mean = 72.5, age_sd = 7.4, age_range = c(51, 94),
                       apoe4_age_shift = -1.7,
                       visits_range = 2:7,
                       visit_prob = c(487, 296, 155, 83, 64, 16) / 1101,
                       gap_mean = 2.05, gap_sd = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ## covariates: 0/1, then centred like the real data
  apoe4  <- rbinom(N, 1, p_apoe4)
  female <- rbinom(N, 1, p_female)
  X <- cbind(apoe4 = apoe4 - x_centre[["apoe4"]], female = female - x_centre[["female"]])

  # baseline age: carriers are younger, shifted so the cohort mean stays age_mean
  mu_age <- age_mean + apoe4_age_shift * (apoe4 - p_apoe4)
  base_age <- rtnorm(N, mu_age, age_sd, age_range[1], age_range[2])
  n_obs <- sample(visits_range, N, replace = TRUE, prob = visit_prob)
  max_n <- max(visits_range)

  # visit ages: baseline + cumulative gaps, with 0-padding
  age_mat <- matrix(0, N, max_n)
  for (i in seq_len(N)) {
    gaps <- rnorm(n_obs[i] - 1L, gap_mean, gap_sd)
    age_mat[i, seq_len(n_obs[i])] <- c(base_age[i], base_age[i] + cumsum(gaps))
  }

  list(n_subj = as.integer(N), n_obs = as.integer(n_obs), age_mat = age_mat,
       X = X, x_centre = x_centre, apoe4 = apoe4, female = female)
}

#######################################
# 4. CURVES implied by the truth
######################################
# Rate curve + clock: bA <- make_rate_basis()
#                     r <- rate_grid(bA, truth$theta)
#                     G <- clock_grid(bA, r, SETTINGS$thresh)
# vA = 4 knots: 0.563, 0.628, 0.794, 1.125 (from fit)
# omega: exp(vA$B %*% truth$nu)

# hardcoded knots based on fit
true_var_knots <- c(0.563, 0.628, 0.794, 1.125)

true_var_basis <- function(kn = true_var_knots, S = SETTINGS) {
  # replicate make_var_basis
  lo <- S$grid_lo; hi <- S$grid_hi
  y  <- seq(lo, hi, length.out = 201L)
  B  <- bSpline(y, knots = kn, degree = 3L, intercept = TRUE,
                Boundary.knots = c(lo, hi))
  list(y = y, n = length(y), step = y[2] - y[1], lo = lo, hi = hi,
       B = B, K = ncol(B), knots = kn)
}

truth_curves <- function(truth) {
  bA <- make_rate_basis()
  r <- rate_grid(bA, truth$theta)
  G <- clock_grid(bA, r, SETTINGS$thresh)
  vA <- true_var_basis()
  omega <- as.numeric(exp(vA$B %*% truth$nu))
  list(bA = bA, r = r, G = G, vA = vA, omega = omega)
}

#########################################
# 5. LATENT QUANTITIES per subject (random)
#########################################
# Draw these for every subject:
#   susceptible: Z_i ~ Bernoulli(pi) -> let Xg be X with intercept column of 1, pi = plogis(Xg %*% phi)
#   onset age: alpha_i = alpha_min + exp(Xg %*% gamma) * E_i^sigma_alpha, E_i ~ Exp(1), truncated (alpha_min, alpha_max)
#   speed: delta_i = X %*% psi + sigma_delta * z_i, z_i ~ N(0, 1),
#   non-acc level b_i ~ N(mu_b, sigma_b) truncated above at thresh
# For accumulator (Z=1)
#   g_ij = G(mu_b) + exp(delta_i) * (age_ij - alpha_i)
#   g_ij = max(g_ij, G(mu_b))            # flat at mu_b before onset
#   A_ij = G^{-1}(g_ij)
# For non accmulator (Z=0): A_ij = b_i at every visit

# Flag any A_ij beyond grid_hi (return -Inf) - make sure not happening frequently

# Accumulator trajectory at each visit; -Inf where the clock runs off the grid
accum_levels <- function(alpha, delta, design, curves, mu_b) {
  bA <- curves$bA
  Gb <- sa_interp(bA$y, curves$G, mu_b)
  A  <- matrix(0, design$n_subj, ncol(design$age_mat))
  for (i in seq_len(design$n_subj)) {
    j <- seq_len(design$n_obs[i])
    g <- pmax(Gb + exp(delta[i]) * (design$age_mat[i, j] - alpha[i]), Gb)
    a <- sa_inv_monotone(bA$y, curves$G, g)
    A[i, j] <- ifelse(is.na(a), -Inf, a)
  }
  A
}

draw_latent <- function(truth, design, curves, S = SETTINGS, max_redraw = 50L) {
  N  <- design$n_subj
  Xg <- cbind(1, design$X)

  pi_i <- as.numeric(plogis(Xg %*% truth$phi))
  Z    <- rbinom(N, 1, pi_i)
  # truncated above at thresh as in the model; below at the grid edge (5 sd away)
  b    <- rtnorm(N, truth$mu_b, truth$sigma_b, S$grid_lo, S$thresh)

  # Weibull AFT onset truncated to (alpha_min, alpha_max] + speed, same draw as rxtilde
  lam <- as.numeric(exp(Xg %*% truth$gamma))
  draw_onset <- function(k) {
    Fmax <- 1 - exp(-((S$alpha_max - S$alpha_min) / lam[k])^(1 / truth$sigma_alpha))
    E <- -log(1 - runif(length(k), 0, Fmax))
    S$alpha_min + lam[k] * E^truth$sigma_alpha
  }
  draw_z <- function(k) rtnorm(length(k), 0, 1, -S$z_bound, S$z_bound)
  speed  <- function(z, k) as.numeric(design$X[k, , drop = FALSE] %*% truth$psi) +
    truth$sigma_delta * z

  all_i <- seq_len(N)
  alpha <- draw_onset(all_i); z <- draw_z(all_i); delta <- speed(z, all_i)
  A_acc <- accum_levels(alpha, delta, design, curves, truth$mu_b)

  # The model gives zero density to a subject whose accumulator trajectory
  # leaves the grid (in either branch), so those are redrawn, and counted.
  off <- which(apply(A_acc, 1, function(a) any(a == -Inf)))
  n_offgrid <- length(off); it <- 0L
  while (length(off) > 0 && it < max_redraw) {
    alpha[off] <- draw_onset(off); z[off] <- draw_z(off); delta[off] <- speed(z[off], off)
    A_acc <- accum_levels(alpha, delta, design, curves, truth$mu_b)
    off <- which(apply(A_acc, 1, function(a) any(a == -Inf))); it <- it + 1L
  }
  if (length(off) > 0) stop(sprintf("%d subjects still off the grid after %d redraws",
                                    length(off), max_redraw))

  # true level at each visit: the trajectory if Z = 1, flat at b_i if Z = 0
  A_mat <- A_acc
  for (i in which(Z == 0)) A_mat[i, seq_len(design$n_obs[i])] <- b[i]

  list(Z = Z, pi = pi_i, alpha = alpha, delta = delta, z = z, b = b,
       A_mat = A_mat, A_acc = A_acc, n_offgrid = n_offgrid)
}

#########################################
# 6. Add measurement noise to get y
#########################################
# y_ij = A_ij + omega(A_ij) * t_ij,   t_ij ~ Student-t(df = t_df = 4)
# omega is a scale
get_y <- function(latent, design, curves, S = SETTINGS) {
  vA <- curves$vA
  y_mat <- matrix(0, design$n_subj, ncol(design$age_mat))
  for (i in seq_len(design$n_subj)) {
    j  <- seq_len(design$n_obs[i])
    A  <- latent$A_mat[i, j]
    om <- sa_interp(vA$y, curves$omega, A)
    y  <- A + om * rt(length(j), df = S$t_df)
    # a t_4 tail can throw a scan off the SUVR grid: redraw those few
    bad <- which(y <= S$grid_lo | y >= S$grid_hi)
    while (length(bad) > 0) {
      y[bad] <- A[bad] + om[bad] * rt(length(bad), df = S$t_df)
      bad <- bad[y[bad] <= S$grid_lo | y[bad] >= S$grid_hi]
    }
    y_mat[i, j] <- y
  }
  y_mat
}

#########################################
# 7. PACKAGE into real data format
#########################################
# Match field for what build_model() reads:
# n_subj, n_obs, y_mat, age_mat, a_tilde, x_init, X, ids, x_centre, thresh
#   a_tilde = mean visit age per subject
#   x_init  = mean simulated y per subject
#unused cells in y_mat / age_mat = 0

get_dat <- function(y_mat, design, S = SETTINGS) {
  N <- design$n_subj
  x_init <- sapply(seq_len(N), function(i) mean(y_mat[i, 1:design$n_obs[i]]))
  a_tilde <- sapply(seq_len(N), function(i) mean(design$age_mat[i, 1:design$n_obs[i]]))
  list(n_subj = N, n_obs = design$n_obs, y_mat = y_mat, age_mat = design$age_mat,
       a_tilde = a_tilde, x_init = x_init, X = design$X, ids = seq_len(N),
       x_centre = design$x_centre, thresh = S$thresh)
}

# The model's own coordinate for the onset: the accumulator level at the visit
# centroid dat$a_tilde, continued linearly below mu_b if onset comes later (see dxtilde)
true_x_tilde <- function(latent, dat, truth, curves, S = SETTINGS) {
  bA <- curves$bA
  Gb <- sa_interp(bA$y, curves$G, truth$mu_b)
  rb <- max(sa_interp(bA$y, curves$r, truth$mu_b), S$r_min)
  clock <- exp(latent$delta) * (dat$a_tilde - latent$alpha)
  ifelse(clock >= 0, sa_inv_monotone(bA$y, curves$G, Gb + pmax(clock, 0)),
         truth$mu_b + rb * clock)
}


#########################################
# 8. CHECK dataset for validity
#########################################
# no NA/non-finite values; all y within the SUVR grid
# share of accumulators ~ pi_bar; share above thresh ~ real data
# onset ages within (alpha_min, alpha_max]
# spaghetti plot vs the real data (check before fitting)

# observed scans of a dat list, in long format
dat_long <- function(dat) {
  i <- rep(seq_len(dat$n_subj), dat$n_obs)
  j <- sequence(dat$n_obs)
  data.frame(id = i, visit = j, age = dat$age_mat[cbind(i, j)], y = dat$y_mat[cbind(i, j)])
}

# HARD checks stop with an error (the dataset is unusable); SOFT checks only
# warn (the dataset is valid but does not look like ADNI -- expected off baseline).
check_dat <- function(dat, latent, real_path = "data/adni_amyloid.rds",
                      S = SETTINGS, verbose = TRUE) {
  res <- data.frame(check = character(), hard = logical(), pass = logical(), detail = character())
  chk <- function(name, pass, detail = "", hard = TRUE)
    res[nrow(res) + 1L, ] <<- list(name, hard, isTRUE(pass), detail)

  N <- dat$n_subj; L <- dat_long(dat)
  obs <- col(dat$y_mat) <= dat$n_obs                      # which cells are real scans
  fu  <- tapply(L$age, L$id, function(a) max(a) - min(a))

  # ---- format
  chk("fields match adni_amyloid.rds",
      identical(names(dat), c("n_subj", "n_obs", "y_mat", "age_mat", "a_tilde",
                              "x_init", "X", "ids", "x_centre", "thresh")))
  chk("dimensions agree", all(dim(dat$y_mat) == dim(dat$age_mat)) && nrow(dat$y_mat) == N &&
        length(dat$n_obs) == N && length(dat$a_tilde) == N && length(dat$x_init) == N &&
        all(dim(dat$X) == c(N, 2)) && length(dat$ids) == N)
  chk("2+ scans per subject, none beyond max_n",
      min(dat$n_obs) >= 2 && max(dat$n_obs) <= ncol(dat$y_mat),
      sprintf("n_obs %d-%d, %d scans", min(dat$n_obs), max(dat$n_obs), sum(dat$n_obs)))
  chk("unused cells are 0", all(dat$y_mat[!obs] == 0) && all(dat$age_mat[!obs] == 0))

  # ---- values
  chk("no NA / non-finite values",
      all(is.finite(dat$y_mat)) && all(is.finite(dat$age_mat)) && all(is.finite(dat$X)) &&
        all(is.finite(dat$a_tilde)) && all(is.finite(dat$x_init)))
  chk("all y inside the SUVR grid", all(L$y > S$grid_lo & L$y < S$grid_hi),
      sprintf("y range %.3f-%.3f", min(L$y), max(L$y)))
  chk("visit ages strictly increasing", all(L$age[-1][diff(L$id) == 0] > L$age[-nrow(L)][diff(L$id) == 0]))
  chk("a_tilde / x_init are subject means",
      max(abs(dat$a_tilde - tapply(L$age, L$id, mean))) < 1e-10 &&
        max(abs(dat$x_init - tapply(L$y, L$id, mean))) < 1e-10)
  chk("covariates are centred 0/1",
      all(round(dat$X[, 1] + dat$x_centre[1], 8) %in% c(0, 1)) &&
        all(round(dat$X[, 2] + dat$x_centre[2], 8) %in% c(0, 1)))

  # ---- latent truth
  A_obs <- latent$A_mat[obs]
  chk("true levels finite and inside the grid",
      all(is.finite(A_obs)) && all(A_obs >= S$grid_lo & A_obs <= S$grid_hi),
      sprintf("A range %.3f-%.3f", min(A_obs), max(A_obs)))
  chk("onset ages within (alpha_min, alpha_max]",
      all(latent$alpha > S$alpha_min & latent$alpha <= S$alpha_max),
      sprintf("alpha range %.1f-%.1f", min(latent$alpha), max(latent$alpha)))
  chk("|z| within z_bound", all(abs(latent$z) < S$z_bound))
  chk("non-accumulator levels below thresh", all(latent$b < S$thresh))
  se <- sqrt(sum(latent$pi * (1 - latent$pi))) / N
  chk("share of accumulators ~ mean pi_i (within 4 SE)",
      abs(mean(latent$Z) - mean(latent$pi)) < 4 * se,
      sprintf("%.3f vs %.3f", mean(latent$Z), mean(latent$pi)))
  chk("off-grid redraws are rare (< 2% of subjects)", latent$n_offgrid / N < 0.02,
      sprintf("%d of %d", latent$n_offgrid, N))

  # ---- soft: does it look like the real cohort?
  chk("follow-up length: median ~4, max ~14", abs(median(fu) - 4) < 1 && max(fu) > 10 && max(fu) < 16,
      sprintf("median %.1f, max %.1f", median(fu), max(fu)), hard = FALSE)
  if (file.exists(real_path)) {
    R <- dat_long(readRDS(real_path))
    ever <- function(D) mean(tapply(D$y > S$thresh, D$id, any))
    chk("share of scans above thresh ~ real (+/- 0.08)",
        abs(mean(L$y > S$thresh) - mean(R$y > S$thresh)) < 0.08,
        sprintf("%.3f vs %.3f", mean(L$y > S$thresh), mean(R$y > S$thresh)), hard = FALSE)
    chk("share of subjects ever above thresh ~ real (+/- 0.08)",
        abs(ever(L) - ever(R)) < 0.08, sprintf("%.3f vs %.3f", ever(L), ever(R)), hard = FALSE)
    chk("median SUVR ~ real (+/- 0.05)", abs(median(L$y) - median(R$y)) < 0.05,
        sprintf("%.3f vs %.3f", median(L$y), median(R$y)), hard = FALSE)
  }

  if (verbose) for (k in seq_len(nrow(res)))
    cat(sprintf("  %s  %s%s\n", if (res$pass[k]) "PASS" else if (res$hard[k]) "FAIL" else "WARN",
                res$check[k], if (nzchar(res$detail[k])) sprintf("  [%s]", res$detail[k]) else ""))
  bad <- res$hard & !res$pass
  if (any(bad)) stop("invalid simulated dataset: ", paste(res$check[bad], collapse = "; "))
  invisible(res)
}

# Simulated vs real, side by side. Look at this BEFORE fitting anything.
plot_sim_check <- function(sim, real_path = "data/adni_amyloid.rds",
                           file = "simulated_data/figures/sim_check_baseline.pdf",
                           S = SETTINGS) {
  suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
  thm <- theme_minimal(base_size = 8.5) +
    theme(plot.title = element_text(face = "bold", size = 9), legend.position = "bottom")
  cols <- c(real = "grey35", simulated = "#B03030")
  real <- readRDS(real_path)
  Ls <- dat_long(sim$dat); Lr <- dat_long(real)
  both <- rbind(cbind(Lr, src = "real"), cbind(Ls, src = "simulated"))
  per_subj <- function(dat, src) {
    L <- dat_long(dat)
    data.frame(src = src, base_age = dat$age_mat[, 1], n_obs = dat$n_obs,
               fu = as.numeric(tapply(L$age, L$id, function(a) max(a) - min(a))),
               slope = as.numeric(by(L, L$id, function(d) coef(lm(y ~ age, d))[2])))
  }
  subj <- rbind(per_subj(real, "real"), per_subj(sim$dat, "simulated"))

  spag <- function(L, col, title)
    ggplot(L, aes(age, y, group = id)) + geom_line(alpha = 0.25, linewidth = 0.25, colour = col) +
    geom_hline(yintercept = S$thresh, linetype = 2, linewidth = 0.3) +
    coord_cartesian(xlim = c(50, 100), ylim = c(0.35, 1.65)) +
    labs(title = title, x = "age", y = "SUVR") + thm
  dens <- function(D, v, xlab, title)
    ggplot(D, aes(.data[[v]], colour = src, fill = src)) + geom_density(alpha = 0.15, linewidth = 0.4) +
    scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
    labs(title = title, x = xlab, y = "density") + thm

  # latent truth: A_ij coloured by branch
  La <- dat_long(within(sim$dat, y_mat <- sim$latent$A_mat))
  La$branch <- ifelse(sim$latent$Z[La$id] == 1, "accumulator", "non-accumulator")
  p_lat <- ggplot(La, aes(age, y, group = id, colour = branch)) +
    geom_line(alpha = 0.3, linewidth = 0.25) +
    geom_hline(yintercept = S$thresh, linetype = 2, linewidth = 0.3) +
    scale_colour_manual(values = c(accumulator = "#B03030", "non-accumulator" = "#2C6FA8"), name = NULL) +
    coord_cartesian(xlim = c(50, 100), ylim = c(0.35, 1.65)) +
    labs(title = "simulated TRUE levels A_ij, by branch", x = "age", y = "SUVR") + thm

  p_nobs <- ggplot(subj, aes(factor(n_obs), fill = src)) + geom_bar(position = "dodge") +
    scale_fill_manual(values = cols, name = NULL) +
    labs(title = "scans per subject", x = "n_obs", y = "subjects") + thm

  cu <- truth_curves(sim$truth)
  p_rate <- ggplot(data.frame(y = cu$bA$y, r = cu$r), aes(y, r)) + geom_line(linewidth = 0.4) +
    geom_vline(xintercept = c(sim$truth$mu_b, S$thresh), linetype = 2, linewidth = 0.3) +
    labs(title = "true rate curve r_A(y)", x = "SUVR (dashed: mu_b, thresh)", y = "SUVR / year") + thm
  p_om <- ggplot(data.frame(y = cu$vA$y, om = cu$omega), aes(y, om)) + geom_line(linewidth = 0.4) +
    labs(title = "true measurement scale omega_A(y)", x = "SUVR", y = "omega") + thm

  fig <- (spag(Lr, cols[["real"]], sprintf("real: %d subjects, %d scans", real$n_subj, nrow(Lr))) |
            spag(Ls, cols[["simulated"]], sprintf("simulated: %d subjects, %d scans", sim$dat$n_subj, nrow(Ls))) |
            p_lat) /
    (dens(both, "y", "SUVR", "all scans") | dens(subj, "base_age", "age", "baseline age") |
       dens(subj, "fu", "years", "follow-up length")) /
    (p_nobs | dens(subj[abs(subj$slope) < 0.1, ], "slope", "SUVR / year", "per-subject OLS slope") |
       p_rate | p_om)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(file, fig, width = 11, height = 9.5)
  invisible(fig)
}

# Simulated TRUTH vs what the fitted model says about the real cohort:
# true SUVR, the shared curves, and the per-subject latent quantities.
#
# The real side has no observed truth, so it comes from the shipped posterior.
# Per-subject quantities are compared WITHIN posterior draws (one thin grey line
# per draw, membership sampled from its posterior probability), never as
# posterior medians: medians shrink the spread and would make a correct
# simulation look too dispersed. Thin red lines are simulated replicates, so the
# two bundles show posterior uncertainty against sampling variability.
plot_truth_check <- function(sim, n_draw = 30, n_rep = 15,
                             file = "simulated_data/figures/truth_check_baseline.pdf",
                             S = SETTINGS, seed = 2) {
  suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
  if (!exists("subject_quantities")) source("run_model.R")
  thm <- theme_minimal(base_size = 8.5) +
    theme(plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 6.8, colour = "grey30"),
          legend.position = "bottom", legend.title = element_blank())
  cols <- c("real (fitted)" = "grey35", simulated = "#B03030")
  bcol <- c(accumulator = "#B03030", "non-accumulator" = "#2C6FA8")
  truth <- sim$truth; cu <- truth_curves(truth); bA <- cu$bA; vA <- cu$vA

  # ---- real side: per-draw latent quantities from the shipped posterior
  fit <- load_fitted_model(); real <- fit$dat; n <- real$n_subj
  sq  <- subject_quantities(fit, n_draw = n_draw, seed = seed)
  pp  <- pool_paired(fit)
  set.seed(seed); d <- sort(sample.int(nrow(pp$g), min(n_draw, nrow(pp$g))))
  stopifnot("draws must line up with subject_quantities" = isTRUE(all.equal(pp$g[d, ], sq$g)))
  K  <- length(d)
  BB <- t(pp$u[d, sprintf("b[%d]", 1:n), drop = FALSE])            # n x K, like sq$alpha
  XT <- t(pp$u[d, sprintf("x_tilde[%d]", 1:n), drop = FALSE])
  ZZ <- matrix(rbinom(n * K, 1, sq$p_susc), n, K)                   # membership, per draw

  # true level at every scan for posterior draw k
  real_levels <- function(k) {
    r <- rate_grid(fit$bA, sq$g[k, sprintf("theta[%d]", 1:fit$bA$K)])
    G <- clock_grid(fit$bA, r, S$thresh); Gb <- sa_interp(fit$bA$y, G, sq$g[k, "mu_b"])
    L <- dat_long(real)
    g <- pmin(pmax(Gb + exp(sq$delta[L$id, k]) * (L$age - sq$alpha[L$id, k]), Gb), max(G))
    L$A <- ifelse(ZZ[L$id, k] == 1, sa_inv_monotone(fit$bA$y, G, g), BB[L$id, k])
    L$branch <- ifelse(ZZ[L$id, k] == 1, "accumulator", "non-accumulator")
    L
  }
  sim_levels <- function(s) {
    L <- dat_long(within(s$dat, y_mat <- s$latent$A_mat)); names(L)[4] <- "A"
    L$branch <- ifelse(s$latent$Z[L$id] == 1, "accumulator", "non-accumulator")
    L
  }

  # ---- simulated side: this dataset plus replicates from the same truth
  sims <- c(list(sim), lapply(seq_len(n_rep - 1L), function(r)
    simulate_one(truth, seed = sim$seed + 7919L * r, verbose = FALSE)))

  # long table of one latent quantity: one row per subject per draw / replicate
  lat <- function(name, real_mat, real_keep, sim_get) rbind(
    do.call(rbind, lapply(seq_len(K), function(k)
      data.frame(q = name, src = "real (fitted)", rep = k, v = real_mat[real_keep[, k], k]))),
    do.call(rbind, lapply(seq_along(sims), function(r)
      data.frame(q = name, src = "simulated", rep = r, v = sim_get(sims[[r]])))))
  acc <- ZZ == 1
  D <- rbind(
    lat("onset age alpha (accumulators)", sq$alpha, acc, function(s) s$latent$alpha[s$latent$Z == 1]),
    lat("speed delta (accumulators)", sq$delta, acc, function(s) s$latent$delta[s$latent$Z == 1]),
    lat("x_tilde: level at centroid (accumulators)", XT, acc, function(s) s$latent$x_tilde[s$latent$Z == 1]),
    lat("level b (non-accumulators)", BB, !acc, function(s) s$latent$b[s$latent$Z == 0]))
  LV <- rbind(
    do.call(rbind, lapply(seq_len(K), function(k) data.frame(src = "real (fitted)", rep = k, v = real_levels(k)$A))),
    do.call(rbind, lapply(seq_along(sims), function(r) data.frame(src = "simulated", rep = r, v = sim_levels(sims[[r]])$A))))

  bundle <- function(DD, title, xlab, sub = NULL)
    ggplot(DD, aes(v, colour = src, group = interaction(src, rep))) +
    geom_line(stat = "density", alpha = 0.35, linewidth = 0.3) +
    scale_colour_manual(values = cols) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 0.8))) +
    labs(title = title, subtitle = sub, x = xlab, y = "density") + thm

  # ---- row 1: true SUVR
  spag <- function(L, title, sub)
    ggplot(L, aes(age, A, group = id, colour = branch)) + geom_line(alpha = 0.3, linewidth = 0.25) +
    geom_hline(yintercept = S$thresh, linetype = 2, linewidth = 0.3) +
    scale_colour_manual(values = bcol) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 0.8))) +
    coord_cartesian(xlim = c(50, 100), ylim = c(0.35, 1.65)) +
    labs(title = title, subtitle = sub, x = "age", y = "true SUVR") + thm
  p_real <- spag(real_levels(1), "A. Real cohort: fitted true SUVR", "One posterior draw; membership sampled.")
  p_sim  <- spag(sim_levels(sim), "B. Simulated: true SUVR A_ij", "This dataset.")
  p_lev  <- bundle(LV, "C. True SUVR at every scan", "true SUVR",
                   sprintf("%d posterior draws vs %d simulated replicates.", K, n_rep))

  # ---- row 2: the shared curves. Band = posterior, line = the simulated truth
  band <- function(M, x) { q <- apply(M, 1, quantile, c(.025, .5, .975), na.rm = TRUE)
    data.frame(x = x, lo = q[1, ], mid = q[2, ], hi = q[3, ]) }
  curve_panel <- function(bd, tr, title, sub, xlab, ylab)
    ggplot(bd, aes(x, mid)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey35", alpha = 0.2) +
    geom_line(aes(colour = "real (fitted)"), linewidth = 0.4) +
    geom_line(data = tr, aes(x, y, colour = "simulated"), linewidth = 0.5) +
    scale_colour_manual(values = cols) +
    labs(title = title, subtitle = sub, x = xlab, y = ylab) + thm
  Rg <- sapply(seq_len(K), function(k) rate_grid(bA, sq$g[k, sprintf("theta[%d]", 1:bA$K)]))
  qd <- quantile(dat_long(real)$y, c(.02, .98))
  p_rate <- curve_panel(band(Rg, bA$y), data.frame(x = bA$y, y = cu$r), "D. Rate curve r_A(y)",
                        "Band: posterior 95%. Grey: outside the real data.", "SUVR", "SUVR / year") +
    annotate("rect", xmin = -Inf, xmax = qd[1], ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.15) +
    annotate("rect", xmin = qd[2], xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.15)
  tt <- seq(0, 60, by = 0.5)
  path <- function(G, mb) sa_inv_monotone(bA$y, G, pmin(sa_interp(bA$y, G, mb) + tt, max(G)))
  Pg <- sapply(seq_len(K), function(k) path(clock_grid(bA, Rg[, k], S$thresh), sq$g[k, "mu_b"]))
  p_path <- curve_panel(band(Pg, tt), data.frame(x = tt, y = path(cu$G, truth$mu_b)),
                        "E. Shared trajectory from the floor", "Level against years since onset, at unit speed.",
                        "years since onset", "true SUVR") +
    geom_hline(yintercept = S$thresh, linetype = 2, linewidth = 0.3)
  Og <- sapply(seq_len(K), function(k) as.numeric(exp(fit$vA$B %*% sq$g[k, sprintf("nu[%d]", 1:fit$vA$K)])))
  p_om <- curve_panel(band(Og, fit$vA$y), data.frame(x = vA$y, y = cu$omega),
                      "F. Measurement scale omega_A(y)", "Band: posterior 95%.", "SUVR", "omega")
  # model-free: the per-subject OLS slope at each level, straight from the scans
  slopes <- function(dat, src) {
    L <- dat_long(dat)
    data.frame(src = src, lev = as.numeric(tapply(L$y, L$id, mean)),
               sl = as.numeric(by(L, L$id, function(d) coef(lm(y ~ age, d))[2])))
  }
  SL <- rbind(slopes(real, "real (fitted)"), do.call(rbind, lapply(sims, function(s) slopes(s$dat, "simulated"))))
  SL$bin <- cut(SL$lev, seq(0.4, 1.6, by = 0.1))
  SB <- do.call(rbind, lapply(split(SL, list(SL$src, SL$bin), drop = TRUE), function(d) if (nrow(d) >= 10)
    data.frame(src = d$src[1], lev = mean(d$lev), mid = median(d$sl),
               lo = quantile(d$sl, .25), hi = quantile(d$sl, .75))))
  SB$src <- ifelse(SB$src == "simulated", "simulated", "real (observed)")
  p_slope <- ggplot(SB, aes(lev, mid, colour = src, fill = src)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.4) + geom_point(size = 0.8) +
    scale_colour_manual(values = c("real (observed)" = "grey35", simulated = "#B03030")) +
    scale_fill_manual(values = c("real (observed)" = "grey35", simulated = "#B03030")) +
    labs(title = "G. Observed slope by level (model-free)", subtitle = "Per-subject OLS slope: median and IQR by mean SUVR.",
         x = "subject mean SUVR", y = "SUVR / year") + thm

  # ---- row 3: latent quantities
  shr <- rbind(data.frame(src = "real (fitted)", v = colMeans(ZZ)),
               data.frame(src = "simulated", v = sapply(sims, function(s) mean(s$latent$Z))))
  p_shr <- ggplot(shr, aes(src, v, colour = src)) + geom_jitter(width = 0.15, height = 0, size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = cols) + guides(colour = "none") +
    labs(title = "H. Share of accumulators", subtitle = "One point per draw / replicate.", x = NULL, y = "share") + thm
  qs <- unique(D$q)
  p_lat <- Map(function(q, lab, xl) bundle(D[D$q == q, ], lab, xl), qs,
               c("I. Onset age alpha", "J. Speed delta", "K. Level at centroid x_tilde", "L. Non-accumulator level b"),
               c("age (accumulators)", "log speed (accumulators)", "SUVR (accumulators)", "SUVR (non-accumulators)"))

  fig <- (p_real | p_sim | p_lev) / (p_rate | p_path | p_om | p_slope) /
    (p_shr | p_lat[[1]] | p_lat[[2]] | p_lat[[3]] | p_lat[[4]]) + plot_layout(heights = c(1.25, 1, 1))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(file, fig, width = 14, height = 11)

  # ---- the same comparison in numbers: mean over draws / replicates of each summary
  summ <- function(v) c(mean = mean(v), sd = sd(v), q05 = unname(quantile(v, .05)), q95 = unname(quantile(v, .95)))
  tab <- do.call(rbind, lapply(split(D, list(D$q, D$src)), function(d)
    data.frame(quantity = d$q[1], src = d$src[1],
               t(rowMeans(sapply(split(d$v, d$rep), summ))))))
  tab <- rbind(tab, data.frame(quantity = "share of accumulators", src = c("real (fitted)", "simulated"),
                               mean = tapply(shr$v, shr$src, mean), sd = tapply(shr$v, shr$src, sd), q05 = NA, q95 = NA))
  tab <- tab[order(tab$quantity, tab$src), ]; rownames(tab) <- NULL
  print(format(tab, digits = 3), row.names = FALSE)
  invisible(list(fig = fig, table = tab))
}

#########################################
# 9. Generate one dataset
#########################################
# Check: plot this against the real data before run simulation
simulate_one <- function(truth, N = NULL, seed, verbose = TRUE) {
  set.seed(seed)
  if (is.null(N)) N <- truth$N
  design <- get_design(N = N)
  curves <- truth_curves(truth)
  latent <- draw_latent(truth, design, curves)
  y_mat  <- get_y(latent, design, curves)
  dat    <- get_dat(y_mat, design)
  latent$x_tilde <- true_x_tilde(latent, dat, truth, curves)
  check_dat(dat, latent, verbose = verbose)
  list(dat = dat, truth = truth, latent = latent, seed = seed)
}

#########################################
# Simulate
#########################################
# Generate different scenarios and replicates
simulate_all <- function(n_rep = 1L, seed0 = 1000L, out_dir = "simulated_data/datasets") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  scenarios <- make_scenarios(get_baseline_truth())
  for (s in seq_along(scenarios)) for (rep in seq_len(n_rep)) {
    nm <- names(scenarios)[s]
    cat(sprintf("\n%s, replicate %d\n", nm, rep))
    sim <- simulate_one(scenarios[[s]], seed = seed0 + 100L * s + rep)
    sim$scenario <- nm; sim$rep <- rep
    saveRDS(sim, file.path(out_dir, sprintf("sim_%s_rep%02d.rds", nm, rep)))
  }
  invisible(names(scenarios))
}

# simulated counterpart of load_amyloid_data(): the dat list build_model() reads
load_simulated_data <- function(scenario = "baseline", rep = 1L, out_dir = "simulated_data/datasets")
  readRDS(file.path(out_dir, sprintf("sim_%s_rep%02d.rds", scenario, rep)))$dat

if (sys.nframe() == 0L) {
  # FIRST: one baseline dataset, checked and plotted against the real data
  cat("baseline check dataset\n")
  sim <- simulate_one(get_baseline_truth(), seed = 1L)
  plot_sim_check(sim)
  cat("-> simulated_data/figures/sim_check_baseline.pdf\n")
  # its TRUTH against the fitted real cohort: true SUVR, curves, latent quantities
  plot_truth_check(sim)
  cat("-> simulated_data/figures/truth_check_baseline.pdf\n")
  # THEN: every scenario
  simulate_all(n_rep = 1L)
}
