## =====================================================================
##  The mixture cure ODE model for amyloid accumulation, in full.
##
##  Flattened deliberately: the research version builds this through a
##  chain of nested files so each layer can be tested against the one
##  below. That is good for development and bad for reading. Everything
##  is here.
##
##  THE MODEL IN WORDS
##
##  Each subject is either a NON-ACCUMULATOR, who sits at a constant level
##  for life, or an ACCUMULATOR, who is flat at a shared floor mu_b until a
##  departure age alpha_i and then follows a shared trajectory. Which one
##  is never observed, so the indicator is marginalised out.
##
##      logit pi_i = w_i' phi                        who accumulates
##      b_i ~ N(mu_b, sigma_b^2), b_i < threshold    those who do not
##      dA/da = exp(delta_i) r_A(A)                  those who do
##      alpha_i ~ Weibull AFT with scale exp(x_i' gamma)   when they start
##      delta_i = x_i' psi + sigma_delta z_i               how fast they go
##
##  WHY THERE IS A CLOCK
##
##  Because r_A depends only on the LEVEL and not on age, the ODE separates
##  and integrates exactly:
##
##      G_A(y) = int_{thresh}^{y} du / max(r_A(u), r_min)
##      G_A(A_i(a)) = G_A(mu_b) + exp(delta_i) (a - alpha_i)
##
##  So the trajectory is a straight line on the clock scale and needs no
##  numerical integration. -G_A(mu_b) is the TRANSIT BUDGET: the time an
##  accumulator must spend getting from the floor to positivity.
##
##  WHERE THE INITIAL CONDITION COMES FROM
##
##  Most progression models give each subject a free starting LEVEL. This
##  one does not: everyone launches from the same mu_b, and what varies is
##  WHEN. The reason is that both cannot be free -- the trajectory depends
##  on a subject-specific start b_i and a departure age only through
##  G_A(b_i) - exp(delta_i) alpha_i, so trading one against the other
##  leaves every fitted value unchanged. Fixing the launch level removes
##  that flat direction, and "this subject is further along than expected"
##  is then explained by a mechanism rather than a free parameter.
## =====================================================================

build_model <- function(dat, bA, vA, S = SETTINGS) {

  code <- nimbleCode({

    ## ---- rate curve: log r_A(y) = B_A(y)' theta -----------------------
    ## theta[1] is FIXED. Its basis support lies entirely below mu_b, where
    ## the model never evaluates the curve, so it is not identified.
    theta[2] ~ dnorm(theta1_mean, sd = theta1_sd)
    sp_th ~ T(dgamma(1, 1), sp_min, sp_max)           # GLOBAL smoothing level
    for (j in 3:K) {
      ## lam is the LOCAL multiplier: it lets the curve bend sharply in one
      ## place without forcing it to bend everywhere. Truncated away from
      ## zero because the increment sd is sqrt(lam/sp), so lam -> 0 sends
      ## that sd to zero and the normal density at its own mean to infinity.
      lam_sp_th[j] ~ T(dexp(1), lam_min, lam_max)
      sd_sp_th[j] <- sqrt(lam_sp_th[j] / sp_th)
      theta[j] ~ dnorm(theta[j-1], sd = sd_sp_th[j])
    }
    r_grid[1:ngrid] <- exp(Bh[1:ngrid, 1:K] %*% theta[1:K])
    G_grid[1:ngrid] <- nfClock(r_grid[1:ngrid], y_step, ngrid, y_lo, thresh, r_min)

    ## ---- measurement scale: fixed-width random walk on log omega ------
    nu[1] ~ T(dnorm(nu1_mean, sd = nu1_sd), nu_lo, nu_hi)
    for (v in 2:Kv) nu[v] ~ T(dnorm(nu[v-1], sd = nu_rw_sd), nu_lo, nu_hi)
    omega_grid[1:nvgrid] <- exp(Bv[1:nvgrid, 1:Kv] %*% nu[1:Kv])

    ## ---- the shared floor --------------------------------------------
    mu_b ~ T(dnorm(mu_b_mean, sd = mu_b_sd), Y_L, thresh)

    ## ---- latency (when accumulation starts) --------------------------
    gamma[1] ~ dnorm(g1_mean, sd = g1_sd)
    for (q in 2:pg) gamma[q] ~ dnorm(0, sd = 0.25)
    sigma_alpha ~ T(dgamma(2, 8), 0, 2)

    ## ---- rate (how fast, once started) -------------------------------
    for (q in 1:p) psi[q] ~ dnorm(0, sd = 0.25)
    sigma_delta ~ T(dnorm(0, sd = 1), 0, sigma_delta_max)

    ## ---- susceptibility (whether at all) -----------------------------
    phi[1] ~ dnorm(phi1_mean, sd = phi1_sd)
    for (q in 2:pw) phi[q] ~ dnorm(0, sd = 0.5)
    sigma_b ~ T(dnorm(0, sd = 0.15), 0, 0.5)
    pi_bar <- 1 / (1 + exp(-phi[1]))          # susceptibility at covariate means

    for (i in 1:n_subj) {
      z[i] ~ T(dnorm(0, sd = 1), -z_bound, z_bound)
      delta[i] <- inprod(X[i, 1:p], psi[1:p]) + sigma_delta * z[i]
      xg[i]    <- inprod(Xg[i, 1:pg], gamma[1:pg])     # AFT log scale
      lp_i[i]  <- inprod(Xg[i, 1:pw], phi[1:pw])       # logit susceptibility

      ## x_tilde is the subject's level at their own visit centroid. It has
      ## NO free prior: its density is implied by the AFT on the departure
      ## age through a change of variables (see R/05_onset.R).
      x_tilde[i] ~ dxtilde(delta = delta[i], atil = a_tilde[i], mu_b = mu_b,
                           xg = xg[i], sigma_alpha = sigma_alpha,
                           alpha_min = alpha_min, alpha_max = alpha_max,
                           G = G_grid[1:ngrid], rgrid = r_grid[1:ngrid],
                           lo = y_lo, step = y_step, ngrid = ngrid, r_min = r_min)

      ## the non-accumulator level, truncated below threshold: a subject
      ## seen above it has certainly accumulated
      b[i] ~ T(dnorm(mu_b, sd = sigma_b), , thresh)

      ## the two branches, combined by logsumexp inside the density
      y_obs[i, 1:max_n] ~ dsubjMix(n_obs = n_obs[i], ages = age_mat[i, 1:max_n],
                           x_tilde = x_tilde[i], a_tilde = a_tilde[i],
                           delta = delta[i], b = b[i], logit_pi = lp_i[i],
                           mu_b = mu_b, G = G_grid[1:ngrid], rgrid = r_grid[1:ngrid],
                           lo = y_lo, step = y_step, ngrid = ngrid, r_min = r_min,
                           omega = omega_grid[1:nvgrid], vlo = v_lo,
                           vstep = v_step, nvgrid = nvgrid, df = t_df)
    }
  })

  Xg <- cbind(1, dat$X)                    # intercept + covariates, already centred
  ## Init the rate curve by least-squares fitting a target log-rate that rises
  ## from the floor to the threshold by a factor of 20 and then FLATTENS. A
  ## curve that keeps rising across the whole grid compresses the clock at high
  ## levels and sends long-follow-up subjects off the end of it -- which shows
  ## up not as a bad fit but as an -Inf at the initial values.
  lr_target <- S$theta1 + log(20) *
    pmin((bA$y - bA$lo) / (S$thresh - bA$lo), 1)
  th0 <- as.numeric(qr.solve(bA$B, lr_target))

  consts <- list(
    n_subj = dat$n_subj, max_n = ncol(dat$y_mat), n_obs = dat$n_obs,
    age_mat = dat$age_mat, a_tilde = dat$a_tilde,
    X = dat$X, p = ncol(dat$X), Xg = Xg, pg = ncol(Xg), pw = ncol(Xg),
    K = bA$K, ngrid = bA$n, Bh = bA$B, y_lo = bA$lo, y_step = bA$step,
    Kv = vA$K, nvgrid = vA$n, Bv = vA$B, v_lo = vA$lo, v_step = vA$step,
    thresh = S$thresh, Y_L = bA$lo, r_min = S$r_min, t_df = S$t_df,
    lam_min = S$lam_min, lam_max = S$lam_max, sp_min = S$sp_min, sp_max = S$sp_max,
    nu_lo = S$nu_lo, nu_hi = S$nu_hi, nu1_mean = S$nu1_mean, nu1_sd = S$nu1_sd,
    nu_rw_sd = S$nu_rw_sd, theta1_mean = S$theta1, theta1_sd = S$theta1_sd,
    sigma_delta_max = S$sigma_delta_max, z_bound = S$z_bound,
    alpha_min = S$alpha_min, alpha_max = S$alpha_max,
    mu_b_mean = S$mu_b_mean, mu_b_sd = S$mu_b_sd,
    phi1_mean = S$phi1_mean, phi1_sd = S$phi1_sd,
    g1_mean = log(70 - S$alpha_min), g1_sd = S$g1_sd)

  data <- list(y_obs = dat$y_mat, theta = c(S$theta1, rep(NA_real_, bA$K - 1)))

  inits <- list(
    theta = th0, sp_th = 1, lam_sp_th = c(NA, NA, rep(1, bA$K - 2)),
    nu = rep(S$nu1_mean, vA$K),
    mu_b = 0.60, gamma = c(log(70 - S$alpha_min), rep(0, ncol(Xg) - 1)),
    sigma_alpha = 0.27, psi = rep(0, ncol(dat$X)), sigma_delta = 0.4,
    phi = c(0.95, rep(0, ncol(Xg) - 1)), sigma_b = 0.10,
    x_tilde = pmin(pmax(dat$x_init, bA$lo + 0.02), bA$hi - 0.02),
    b = pmin(dat$x_init, S$thresh - 0.01),
    z = rep(0, dat$n_subj))
  inits$theta[1] <- S$theta1

  ## FEASIBLE z INITS.  Since
  ##     alpha_i = atil_i - e^{-delta_i} { G_ext(x_i) - G(mu_b) }
  ## must exceed alpha_min, a subject seen at a high level but not yet very
  ## old requires a MINIMUM speed:
  ##     delta_i > log[ {G_ext(x_i) - G(mu_b)} / (atil_i - alpha_min) ].
  ## At delta = 0 a small fraction of subjects violate this -- the ones whose
  ## lever arm exceeds their own age -- and each contributes -Inf, so the model
  ## cannot be initialised at all. Starting them at the required speed is an
  ## INITIALISATION, not a prior: it changes where the chain begins, not what
  ## it targets.
  r0  <- as.numeric(exp(bA$B %*% inits$theta))
  G0  <- sa_cumtrapz(1 / pmax(r0, S$r_min), bA$step)
  G0  <- G0 - sa_interp(bA$y, G0, S$thresh)
  Gb0 <- sa_interp(bA$y, G0, inits$mu_b)
  rb0 <- max(sa_interp(bA$y, r0, inits$mu_b), S$r_min)
  g0v <- ifelse(inits$x_tilde >= inits$mu_b,
                sa_interp(bA$y, G0, inits$x_tilde),
                Gb0 + (inits$x_tilde - inits$mu_b) / rb0)
  lever <- pmax(g0v - Gb0, 1e-8)
  room  <- pmax(dat$a_tilde - S$alpha_min, 1e-8)
  dini  <- pmax(log(lever / room) + 0.15, 0)
  inits$z <- pmin(pmax(dini / inits$sigma_delta, -S$z_bound + 0.01), S$z_bound - 0.01)

  list(code = code, constants = consts, data = data, inits = inits)
}
