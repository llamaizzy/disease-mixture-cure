## =====================================================================
##  Truth scenarios for the simulation study.
##
##  Every scenario starts from the posterior medians of the shipped fit and
##  changes ONE thing, so a difference in recovery between scenarios has one
##  cause. Truth is always expressed in the SAME bases the model fits with
##  (theta on the rate basis, nu on the variance basis); nothing is
##  simulated from a curve the model cannot represent, so any bias found
##  here is the estimator's, not a misspecification's.
## =====================================================================

SIM_SCENARIOS <- c("fitted", "rate_only", "low_cure", "null")

sim_truth <- function(scenario = SIM_SCENARIOS, post = NULL) {
  scenario <- match.arg(scenario)
  if (is.null(post)) post <- readRDS("data/fitted_posterior.rds")$g
  med <- apply(post, 2, median)
  tr <- list(
    theta       = unname(med[sprintf("theta[%d]", 1:10)]),
    nu          = unname(med[sprintf("nu[%d]", 1:8)]),
    mu_b        = unname(med["mu_b"]),
    gamma       = unname(med[sprintf("gamma[%d]", 1:3)]),
    sigma_alpha = unname(med["sigma_alpha"]),
    psi         = unname(med[sprintf("psi[%d]", 1:2)]),
    sigma_delta = unname(med["sigma_delta"]),
    phi         = unname(med[sprintf("phi[%d]", 1:3)]),
    sigma_b     = unname(med["sigma_b"]))
  tr$theta[1] <- SETTINGS$theta1          # pinned in the model; pinned in the truth

  switch(scenario,
    ## the model's own answer: APOE4 acts on WHETHER (OR 5.3) and WHEN (TR 0.73),
    ## not on HOW FAST (RR 1.09). Parameter recovery at the real design.
    fitted = tr,

    ## the headline claim inverted: APOE4 carriers depart at the SAME age but
    ## accumulate 50% faster. If the model reports this as "earlier" it cannot
    ## tell latency from rate, and the real-data headline is not trustworthy.
    rate_only = { tr$gamma[2] <- 0; tr$psi[1] <- log(1.5); tr },

    ## almost everyone is susceptible (pi = 0.95, no covariate effect on it).
    ## The question is whether the estimated cure fraction stays near 5% or
    ## gets inflated -- the failure mode the extended clock was built against.
    ## 0.95 rather than 1 because phi[1] has a N(0.95, 1) prior; a truth 5 sd
    ## into that prior's tail would test the prior, not the sampler.
    low_cure = { tr$phi[1] <- qlogis(0.95); tr$phi[2:3] <- 0; tr },

    ## no covariate does anything. Every contrast interval should cover 0
    ## about 95% of the time; this is the false-positive check.
    null = { tr$gamma[2:3] <- 0; tr$psi[] <- 0; tr$phi[2:3] <- 0; tr })
}

## The quantities the study scores, on the scales the report uses. Truth
## values that depend on the REALISED cohort (mean crossing age among the
## subjects who actually are accumulators) are filled in by sim_generate().
sim_estimands <- function(truth, bA, S = SETTINGS, y_eval = seq(0.6, 1.5, by = 0.1)) {
  r  <- rate_grid(bA, truth$theta); G <- clock_grid(bA, r, S$thresh)
  Gb <- sa_interp(bA$y, G, truth$mu_b)
  est <- c(
    OR_apoe4 = exp(truth$phi[2]), OR_female = exp(truth$phi[3]), pi_bar = plogis(truth$phi[1]),
    TR_apoe4 = exp(truth$gamma[2]), TR_female = exp(truth$gamma[3]),
    RR_apoe4 = exp(truth$psi[1]),   RR_female = exp(truth$psi[2]),
    mu_b = truth$mu_b, sigma_b = truth$sigma_b, sigma_alpha = truth$sigma_alpha,
    sigma_delta = truth$sigma_delta, gamma1 = truth$gamma[1],
    transit = -Gb)
  lr <- setNames(sa_interp(bA$y, log(r), y_eval), sprintf("logr_%.1f", y_eval))
  c(est, lr)
}
