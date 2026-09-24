## =====================================================================
##  Fixed settings. Everything a student might reasonably want to change
##  lives here, in one place, with a comment saying what it does.
## =====================================================================
SETTINGS <- list(
  ## --- the amyloid level grid (SUVR). Trajectories and the clock live here.
  grid_lo    = 0.30,
  grid_hi    = 1.70,
  grid_step  = 0.002,
  thresh     = 0.7464406,   # positivity threshold = 24 Centiloid

  ## --- rate curve r_A(y): a cubic B-spline in LEVEL
  K_rate     = 10L,
  r_min      = 1e-4,        # numerical floor in the clock integrand

  ## --- smoothing on the rate curve (adaptive random walk)
  lam_min    = 0.01, lam_max = 100,   # truncation on the LOCAL multiplier
  sp_min     = 0.01, sp_max  = 100,   # truncation on the GLOBAL level
  theta1     = -6.9,        # theta[1] is fixed: its support is below mu_b
  theta1_sd  = 2.0,         # prior sd for theta[2], which starts the walk

  ## --- measurement scale omega_A(y): a spline in LEVEL, fixed-width walk
  K_var      = 5L,          # -> 4 interior knots, basis of 8 columns
  var_knot_q = c(0.15, 0.85),
  nu1_mean   = -3.4,        # exp(-3.4) = 0.033, the assay repeatability floor
  nu1_sd     = 1.0,
  nu_rw_sd   = 0.5,
  nu_lo      = -9, nu_hi = 1,
  t_df       = 4,           # Student-t degrees of freedom for residuals

  ## --- the onset layer
  alpha_min  = 25,  alpha_max = 125,   # departure age bounds
  mu_b_mean  = 0.45, mu_b_sd  = 0.15,  # prior on the shared floor
  phi1_mean  = 0.95, phi1_sd  = 1.0,   # susceptibility intercept, logit scale
  g1_sd      = 0.75,                   # AFT intercept prior sd
  sigma_delta_max = 2,
  z_bound    = 4
)
