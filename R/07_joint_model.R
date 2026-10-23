## =====================================================================
##  The joint amyloid-tau mixture cure ODE model, in full.
##  (docs/context/coupling_and_copula.pdf; densities in R/06_joint_densities.R)
##
##  THE MODEL IN WORDS
##
##  Each subject is in one of three lifetime states, never observed and
##  marginalised in the sampler:
##      00  accumulates neither        weight (1 - pi_i)
##      10  amyloid only               weight pi_i (1 - q_i)
##      11  amyloid then tau           weight pi_i q_i
##  There is no tau-only state (Z_T <= Z_A).
##
##      logit pi_i = w_i' phi_A           amyloid susceptibility
##      logit q_i  = w_i' phi_T           tau susceptibility GIVEN amyloid
##
##  Amyloid is exactly the model in R/04_model.R. Tau has the same
##  structure (own floor, own rate curve, own clock, own AFT departure and
##  speed) with two links to amyloid:
##
##   * a Gaussian copula with correlation rho on the two departure ages,
##     with both Weibull AFT marginals preserved exactly;
##   * the tau rate is multiplied by the subject's CURRENT amyloid burden,
##         dT/da = e^{delta_T} m_i(a) rbar_T(T),
##         m_i(a) = exp{ beta [H_A(A_i(a)) - 1/2] },  H_A = expit((y - c_A)/d),
##     so m_i = 1 at the amyloid threshold, and beta = 0 switches it off.
##
##  Residual speeds z_A, z_T are independent; neither speed regression has
##  an intercept, which is what fixes the scale of each rate curve.
##
##  The amyloid-only model is untouched: this file adds a second builder.
## =====================================================================

build_joint_model <- function(dat, bA, vA, bT, vT,
                              SA = SETTINGS, ST = SETTINGS_TAU, SJ = SETTINGS_JOINT) {

  code <- nimbleCode({

    ## ================= amyloid population (as R/04_model.R) =============
    theta[2] ~ dnorm(theta1_mean, sd = theta1_sd)
    sp_th ~ T(dgamma(1, 1), sp_min, sp_max)
    for (j in 3:K) {
      lam_sp_th[j] ~ T(dexp(1), lam_min, lam_max)
      sd_sp_th[j] <- sqrt(lam_sp_th[j] / sp_th)
      theta[j] ~ dnorm(theta[j-1], sd = sd_sp_th[j])
    }
    r_grid[1:ngrid] <- exp(Bh[1:ngrid, 1:K] %*% theta[1:K])
    G_grid[1:ngrid] <- nfClock(r_grid[1:ngrid], y_step, ngrid, y_lo, thresh, r_min)
    ## THE EXPOSURE CLOCK: how much tau-driving exposure an accumulator
    ## collects per unit of amyloid level (see R/06_joint_densities.R)
    E_grid[1:ngrid] <- nfExpClock(r_grid[1:ngrid], y_step, ngrid, y_lo, r_min,
                                  beta, thresh, d_coup)

    nu[1] ~ T(dnorm(nu1_mean, sd = nu1_sd), nu_lo, nu_hi)
    for (v in 2:Kv) nu[v] ~ T(dnorm(nu[v-1], sd = nu_rw_sd), nu_lo, nu_hi)
    omega_grid[1:nvgrid] <- exp(Bv[1:nvgrid, 1:Kv] %*% nu[1:Kv])

    mu_b ~ T(dnorm(mu_b_mean, sd = mu_b_sd), Y_L, thresh)
    gamma[1] ~ dnorm(g1_mean, sd = g1_sd)
    for (q in 2:pg) gamma[q] ~ dnorm(0, sd = 0.25)
    sigma_alpha ~ T(dgamma(2, 8), 0, 2)
    for (q in 1:p) psi[q] ~ dnorm(0, sd = 0.25)
    sigma_delta ~ T(dnorm(0, sd = 1), 0, sigma_delta_max)
    phi[1] ~ dnorm(phi1_mean, sd = phi1_sd)
    phi[2] ~ dnorm(0.5, sd = 1)
    phi[3] ~ dnorm(0.1, 0.7)
    sigma_b ~ T(dnorm(0, sd = 0.15), 0, 0.5)
    pi_bar <- 1 / (1 + exp(-phi[1]))

    ## ================= tau population ===================================
    ## theta_T[1] is free (see SETTINGS_TAU): the random walk starts at it
    theta_T[1] ~ dnorm(theta1_mean_T, sd = theta1_sd_T)
    sp_th_T ~ T(dgamma(1, 1), sp_min, sp_max)
    for (j in 2:K_T) {
      lam_sp_th_T[j] ~ T(dexp(1), lam_min, lam_max)
      sd_sp_th_T[j] <- sqrt(lam_sp_th_T[j] / sp_th_T)
      theta_T[j] ~ dnorm(theta_T[j-1], sd = sd_sp_th_T[j])
    }
    r_grid_T[1:ngrid_T] <- exp(Bh_T[1:ngrid_T, 1:K_T] %*% theta_T[1:K_T])
    G_grid_T[1:ngrid_T] <- nfClock(r_grid_T[1:ngrid_T], y_step_T, ngrid_T, y_lo_T,
                                   thresh_T, r_min)

    nu_T[1] ~ T(dnorm(nu1_mean_T, sd = nu1_sd), nu_lo, nu_hi)
    for (v in 2:Kv_T) nu_T[v] ~ T(dnorm(nu_T[v-1], sd = nu_rw_sd), nu_lo, nu_hi)
    omega_grid_T[1:nvgrid_T] <- exp(Bv_T[1:nvgrid_T, 1:Kv_T] %*% nu_T[1:Kv_T])

    mu_b_T ~ T(dnorm(mu_b_mean_T, sd = mu_b_sd_T), y_lo_T, thresh_T)
    gamma_T[1] ~ dnorm(g1_mean_T, sd = g1_sd)
    for (q in 2:pg) gamma_T[q] ~ dnorm(0, sd = 0.25)
    sigma_alpha_T ~ T(dgamma(2, 8), 0, 2)
    for (q in 1:p) psi_T[q] ~ dnorm(0, sd = 0.25)
    sigma_delta_T ~ T(dnorm(0, sd = 1), 0, sigma_delta_max)
    phi_T[1] ~ dnorm(phi1_mean_T, sd = phi1_sd_T)
    for (q in 2:pg) phi_T[q] ~ dnorm(0, sd = 0.5)
    sigma_b_T ~ T(dnorm(0, sd = sigma_b_sd_T), 0, sigma_b_max_T)
    q_bar <- 1 / (1 + exp(-phi_T[1]))         # P(tau | amyloid) at covariate means

    ## ================= coupling and dependence ==========================
    beta ~ T(dnorm(0, sd = beta_sd), 0, )     # amyloid -> tau rate multiplier
    u_rho ~ dbeta(rho_a, rho_b)
    rho_alpha <- 2 * u_rho - 1                # departure-age copula correlation

    ## ================= subjects =========================================
    for (i in 1:n_subj) {
      ## amyloid, as before
      z[i] ~ T(dnorm(0, sd = 1), -z_bound, z_bound)
      delta[i] <- inprod(X[i, 1:p], psi[1:p]) + sigma_delta * z[i]
      xg[i]    <- inprod(Xg[i, 1:pg], gamma[1:pg])
      lp_i[i]  <- inprod(Xg[i, 1:pg], phi[1:pg])
      x_tilde[i] ~ dxtilde(delta = delta[i], atil = a_tilde[i], mu_b = mu_b,
                           xg = xg[i], sigma_alpha = sigma_alpha,
                           alpha_min = alpha_min, alpha_max = alpha_max,
                           G = G_grid[1:ngrid], rgrid = r_grid[1:ngrid],
                           lo = y_lo, step = y_step, ngrid = ngrid, r_min = r_min)
      b[i] ~ T(dnorm(mu_b, sd = sigma_b), , thresh)

      ## tau
      z_T[i] ~ T(dnorm(0, sd = 1), -z_bound, z_bound)
      delta_T[i] <- inprod(X[i, 1:p], psi_T[1:p]) + sigma_delta_T * z_T[i]
      xg_T[i]    <- inprod(Xg[i, 1:pg], gamma_T[1:pg])
      lq_i[i]    <- inprod(Xg[i, 1:pg], phi_T[1:pg])
      ## the tau location coordinate: its density is the AFT marginal, the
      ## copula with THIS subject's amyloid departure age, and the exposure
      ## clock Jacobian -- all in one place, added once
      x_tilde_T[i] ~ dxtildeT(
        deltaT = delta_T[i], atilT = a_tilde_T[i], mu_bT = mu_b_T, xgT = xg_T[i],
        sigma_alphaT = sigma_alpha_T, aminT = alpha_min, amaxT = alpha_max,
        GT = G_grid_T[1:ngrid_T], rT = r_grid_T[1:ngrid_T], loT = y_lo_T,
        stepT = y_step_T, ngridT = ngrid_T, r_minT = r_min,
        xA = x_tilde[i], deltaA = delta[i], atilA = a_tilde[i], mu_bA = mu_b,
        xgA = xg[i], sigma_alphaA = sigma_alpha, aminA = alpha_min, amaxA = alpha_max,
        GA = G_grid[1:ngrid], rA = r_grid[1:ngrid], EA = E_grid[1:ngrid],
        loA = y_lo, stepA = y_step, ngridA = ngrid, r_minA = r_min,
        beta = beta, cA = thresh, d_coup = d_coup, rho = rho_alpha)
      b_T[i] ~ T(dnorm(mu_b_T, sd = sigma_b_T), , thresh_T)

      ## the three states, combined by logsumexp inside the density
      y_joint[i, 1:max_nn] ~ dsubjJoint(
        nA = n_obs[i], agesA = age_mat[i, 1:max_n], nT = n_obs_T[i],
        agesT = age_mat_T[i, 1:max_n_T], offT = max_n,
        xA = x_tilde[i], atilA = a_tilde[i], deltaA = delta[i], bA = b[i],
        logit_pi = lp_i[i], mu_bA = mu_b,
        GA = G_grid[1:ngrid], rA = r_grid[1:ngrid], EA = E_grid[1:ngrid],
        loA = y_lo, stepA = y_step, ngridA = ngrid, r_minA = r_min,
        omA = omega_grid[1:nvgrid], vloA = v_lo, vstepA = v_step, nvA = nvgrid,
        xT = x_tilde_T[i], atilT = a_tilde_T[i], deltaT = delta_T[i], bT = b_T[i],
        logit_q = lq_i[i], mu_bT = mu_b_T,
        GT = G_grid_T[1:ngrid_T], rT = r_grid_T[1:ngrid_T], loT = y_lo_T,
        stepT = y_step_T, ngridT = ngrid_T, r_minT = r_min,
        omT = omega_grid_T[1:nvgrid_T], vloT = v_lo_T, vstepT = v_step_T, nvT = nvgrid_T,
        beta = beta, cA = thresh, d_coup = d_coup, df = t_df)
    }
  })

  Xg <- cbind(1, dat$X)
  ## rate-curve inits: rise by a factor of 20 from the floor to the threshold,
  ## then flatten (see R/04_model.R for why)
  init_theta <- function(b, S) {
    lr <- S$theta1 + log(S$init_rise) * pmin((b$y - b$lo) / (S$thresh - b$lo), 1)
    as.numeric(qr.solve(b$B, lr))
  }
  thA <- init_theta(bA, SA); thA[1] <- SA$theta1
  thT <- init_theta(bT, ST)

  ## tau scans need at least one column even if a subset has none
  max_n_T <- max(1L, ncol(dat$y_tau_mat))
  y_tau   <- dat$y_tau_mat[, seq_len(max_n_T), drop = FALSE]
  age_tau <- dat$age_tau_mat[, seq_len(max_n_T), drop = FALSE]

  consts <- list(
    n_subj = dat$n_subj, n_obs = dat$n_obs, max_n = ncol(dat$y_mat),
    age_mat = dat$age_mat, a_tilde = dat$a_tilde,
    n_obs_T = dat$n_obs_tau, max_n_T = max_n_T, max_nn = ncol(dat$y_mat) + max_n_T,
    age_mat_T = age_tau, a_tilde_T = dat$a_tilde_tau,
    X = dat$X, p = ncol(dat$X), Xg = Xg, pg = ncol(Xg),
    ## amyloid grids and priors
    K = bA$K, ngrid = bA$n, Bh = bA$B, y_lo = bA$lo, y_step = bA$step,
    Kv = vA$K, nvgrid = vA$n, Bv = vA$B, v_lo = vA$lo, v_step = vA$step,
    thresh = SA$thresh, Y_L = bA$lo, r_min = SA$r_min, t_df = SA$t_df,
    lam_min = SA$lam_min, lam_max = SA$lam_max, sp_min = SA$sp_min, sp_max = SA$sp_max,
    nu_lo = SA$nu_lo, nu_hi = SA$nu_hi, nu1_mean = SA$nu1_mean, nu1_sd = SA$nu1_sd,
    nu_rw_sd = SA$nu_rw_sd, theta1_mean = SA$theta1, theta1_sd = SA$theta1_sd,
    sigma_delta_max = SA$sigma_delta_max, z_bound = SA$z_bound,
    alpha_min = SA$alpha_min, alpha_max = SA$alpha_max,
    mu_b_mean = SA$mu_b_mean, mu_b_sd = SA$mu_b_sd,
    phi1_mean = SA$phi1_mean, phi1_sd = SA$phi1_sd,
    g1_mean = log(70 - SA$alpha_min), g1_sd = SA$g1_sd,
    ## tau grids and priors
    K_T = bT$K, ngrid_T = bT$n, Bh_T = bT$B, y_lo_T = bT$lo, y_step_T = bT$step,
    Kv_T = vT$K, nvgrid_T = vT$n, Bv_T = vT$B, v_lo_T = vT$lo, v_step_T = vT$step,
    thresh_T = ST$thresh, theta1_mean_T = ST$theta1, theta1_sd_T = ST$theta1_sd,
    nu1_mean_T = ST$nu1_mean, mu_b_mean_T = ST$mu_b_mean, mu_b_sd_T = ST$mu_b_sd,
    phi1_mean_T = ST$phi1_mean, phi1_sd_T = ST$phi1_sd,
    g1_mean_T = log(ST$g1_age - SA$alpha_min),
    sigma_b_sd_T = ST$sigma_b_sd, sigma_b_max_T = ST$sigma_b_max,
    ## coupling
    d_coup = SJ$d_coup, beta_sd = SJ$beta_sd,
    rho_a = SJ$rho_beta[1], rho_b = SJ$rho_beta[2])

  data <- list(y_joint = cbind(dat$y_mat, y_tau),
               theta = c(SA$theta1, rep(NA_real_, bA$K - 1)))

  inits <- list(
    theta = thA, sp_th = 1, lam_sp_th = c(NA, NA, rep(1, bA$K - 2)),
    nu = rep(SA$nu1_mean, vA$K),
    mu_b = 0.60, gamma = c(log(70 - SA$alpha_min), rep(0, ncol(Xg) - 1)),
    sigma_alpha = 0.27, psi = rep(0, ncol(dat$X)), sigma_delta = 0.4,
    phi = c(0.95, rep(0, ncol(Xg) - 1)), sigma_b = 0.10,
    x_tilde = pmin(pmax(dat$x_init, bA$lo + 0.02), bA$hi - 0.02),
    b = pmin(dat$x_init, SA$thresh - 0.01),
    z = rep(0, dat$n_subj),
    theta_T = thT, sp_th_T = 1, lam_sp_th_T = c(NA, rep(1, bT$K - 1)),
    nu_T = rep(ST$nu1_mean, vT$K),
    mu_b_T = min(ST$mu_b_mean, ST$thresh - 0.05),
    gamma_T = c(log(ST$g1_age - SA$alpha_min), rep(0, ncol(Xg) - 1)),
    sigma_alpha_T = 0.27, psi_T = rep(0, ncol(dat$X)), sigma_delta_T = 0.4,
    phi_T = rep(0, ncol(Xg)), sigma_b_T = 0.08,
    x_tilde_T = pmin(pmax(dat$x_init_tau, bT$lo + 0.02), bT$hi - 0.02),
    b_T = pmin(dat$x_init_tau, ST$thresh - 0.01),
    z_T = rep(0, dat$n_subj),
    beta = 0, u_rho = 0.5)                    # start UNCOUPLED and independent

  ## FEASIBLE SPEED INITS, for both biomarkers (see R/04_model.R). With
  ## beta = 0 the tau map is the plain clock map, so the same bound applies
  ## with the tau clock and tau reference age.
  feasible_z <- function(b, S, theta, mu_b, x_tilde, a_tilde, sigma_delta) {
    r0  <- as.numeric(exp(b$B %*% theta))
    G0  <- clock_grid(b, r0, S$thresh, S$r_min)
    Gb0 <- sa_interp(b$y, G0, mu_b)
    rb0 <- max(sa_interp(b$y, r0, mu_b), S$r_min)
    g0v <- ifelse(x_tilde >= mu_b, sa_interp(b$y, G0, x_tilde),
                  Gb0 + (x_tilde - mu_b) / rb0)
    lever <- pmax(g0v - Gb0, 1e-8)
    room  <- pmax(a_tilde - S$alpha_min, 1e-8)
    dini  <- pmax(log(lever / room) + 0.15, 0)
    pmin(pmax(dini / sigma_delta, -S$z_bound + 0.01), S$z_bound - 0.01)
  }
  inits$z   <- feasible_z(bA, SA, inits$theta, inits$mu_b, inits$x_tilde,
                          dat$a_tilde, inits$sigma_delta)
  inits$z_T <- feasible_z(bT, ST, inits$theta_T, inits$mu_b_T, inits$x_tilde_T,
                          dat$a_tilde_tau, inits$sigma_delta_T)

  list(code = code, constants = consts, data = data, inits = inits)
}
