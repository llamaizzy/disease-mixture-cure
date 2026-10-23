## =====================================================================
##  The data-generating process: simulate one cohort from the model.
##
##  The DESIGN is real and fixed: the ADNI visit ages, number of scans per
##  subject and covariates are kept exactly as observed, and only the
##  amyloid measurements are simulated. Two reasons. First, the question
##  is whether the model recovers truth under THIS study's follow-up
##  pattern (median four years onto a twenty-year process), not under an
##  idealised one. Second, a fixed design is what lets one compiled model
##  serve every replicate -- which is where all the compute savings come
##  from (see sim_fit.R).
## =====================================================================

## Draw the whole cohort: membership, departure age, speed, non-accumulator
## level, and the scans. Returns the data list the model fits (same shape as
## load_amyloid_data()) plus the subject-level truth for scoring.
simulate_cohort <- function(dat, truth, bA, vA, seed, S = SETTINGS) {
  set.seed(seed)
  n  <- dat$n_subj; Xg <- cbind(1, dat$X); max_n <- ncol(dat$y_mat)
  r  <- rate_grid(bA, truth$theta); G <- clock_grid(bA, r, S$thresh)
  Gb <- sa_interp(bA$y, G, truth$mu_b)
  om <- as.numeric(exp(vA$B %*% truth$nu))
  omega_at <- function(y) sa_interp(vA$y, om, pmin(pmax(y, vA$lo), vA$hi))

  ## WHO: susceptibility
  pi_i <- plogis(as.numeric(Xg %*% truth$phi))
  susc <- runif(n) < pi_i

  ## the non-accumulators' lifelong level: N(mu_b, sigma_b) truncated above at
  ## the threshold (a subject above it has certainly accumulated)
  b <- truth$mu_b + truth$sigma_b *
    qnorm(runif(n) * pnorm((S$thresh - truth$mu_b) / truth$sigma_b))

  ## WHEN: Weibull AFT on the departure age, truncated to (alpha_min, alpha_max]
  lam  <- exp(as.numeric(Xg %*% truth$gamma)); kap <- 1 / truth$sigma_alpha
  Fmax <- 1 - exp(-((S$alpha_max - S$alpha_min) / lam)^kap)
  alpha <- S$alpha_min + lam * (-log(1 - runif(n) * Fmax))^truth$sigma_alpha

  ## HOW FAST: non-centred, z truncated at +-z_bound as in the model
  z     <- qnorm(runif(n, pnorm(-S$z_bound), pnorm(S$z_bound)))
  delta <- as.numeric(dat$X %*% truth$psi) + truth$sigma_delta * z

  ## the trajectories, on the clock; flat at mu_b before departure
  A <- y <- matrix(0, n, max_n)
  for (i in seq_len(n)) {
    no <- dat$n_obs[i]; ages <- dat$age_mat[i, 1:no]
    if (susc[i]) {
      gt <- pmin(pmax(Gb + exp(delta[i]) * (ages - alpha[i]), Gb), max(G))
      mu <- sa_inv_monotone(bA$y, G, gt)
    } else mu <- rep(b[i], no)
    A[i, 1:no] <- mu
    y[i, 1:no] <- mu + omega_at(mu) * rt(no, df = S$t_df)   # scaled t_4, as fitted
  }

  ## state at the visit centroid: C cured, L susceptible-not-departed, P progressing
  state <- ifelse(!susc, "C", ifelse(alpha > dat$a_tilde, "L", "P"))

  dat_sim <- dat
  dat_sim$y_mat  <- y
  dat_sim$x_init <- vapply(seq_len(n), function(i) mean(y[i, 1:dat$n_obs[i]]), 0)

  ## realised-cohort truths for the derived estimands
  apoe  <- dat$X[, 1] + dat$x_centre[1]
  cross <- alpha + (-Gb) * exp(-delta)
  realised <- c(
    cross_apoe0 = mean(cross[susc & apoe == 0]),
    cross_apoe1 = mean(cross[susc & apoe == 1]),
    cure_apoe0  = 1 - plogis(sum(c(1, 0 - dat$x_centre[1], 0) * truth$phi)),
    cure_apoe1  = 1 - plogis(sum(c(1, 1 - dat$x_centre[1], 0) * truth$phi)),
    frac_susc   = mean(susc), frac_P = mean(state == "P"), frac_L = mean(state == "L"))

  list(dat = dat_sim,
       subj = list(susc = susc, state = state, alpha = alpha, delta = delta,
                   z = z, b = b, A = A, cross = cross),
       realised = realised)
}

## The cohort design used by a tier: the full ADNI design, or a fixed random
## subsample of SUBJECTS (never visits). seed = 7 matches run_model.R's
## "smoke" tier, so n_sub = 250 gives the identical subset.
sim_design <- function(n_sub = NA_integer_, seed = 7L) {
  dat <- load_amyloid_data()
  if (!is.na(n_sub) && n_sub < dat$n_subj) {
    set.seed(seed); k <- sort(sample(dat$n_subj, n_sub))
    dat <- subset_subjects(dat, k)
  }
  dat
}

## The variance basis is built from the FULL real cohort's scans for every
## tier and every replicate, so that "the truth nu" means the same function
## regardless of subsample or simulated outcome. (run_model.R rebuilds it from
## the data being fitted; on simulated data that would move the knots with
## each replicate and make nu incomparable across them.)
sim_var_basis <- function() {
  dat <- load_amyloid_data()
  make_var_basis(unlist(lapply(seq_len(dat$n_subj), function(i) dat$y_mat[i, 1:dat$n_obs[i]])))
}
