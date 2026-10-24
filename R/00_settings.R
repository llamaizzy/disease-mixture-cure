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

## =====================================================================
##  The joint amyloid-tau model (docs/context/coupling_and_copula.pdf).
##
##  Tau gets its own grid, threshold, bases and prior scales. Everything
##  marked PROVISIONAL is a calibration choice for the mesial temporal
##  meta-ROI SUVR summary in data/adni_joint.rds that has not been settled
##  scientifically; it is here so that it is frozen and visible, not
##  because it is known.
## =====================================================================
## ---------------------------------------------------------------------
##  Tau = the NEO-T composite: volume-weighted [18F]flortaucipir SUVR of
##  bilateral inferior and middle temporal gyri (Ossenkoppele et al. 2022,
##  Nat Med 28:2381), inferior cerebellar grey reference, no partial-volume
##  correction, 6 mm (the UC Berkeley ADNI pipeline). Priors below are set
##  from the literature on THAT measurement, not from this dataset:
##
##   [J20] Jack et al. 2020, JAMA Neurol 77:1-10 (doi 10.1001/jamaneurol.
##         2020.3156): flortaucipir, cerebellar-crus reference, no PVC.
##         Inferior temporal SUVR in cognitively unimpaired with < 8 CL:
##         ADNI 1.184 (SD 0.078), MCSA 1.193 (0.094). Annualised change:
##         0.004-0.005 /yr at 9-67 CL, 0.030 (SD 0.03-0.05) /yr at > 68 CL.
##   [O22] Ossenkoppele et al. 2022: NEO-T positivity = mean + 2 SD of the
##         amyloid-negative participants in each cohort; 65 of 448 (14.5%)
##         A+ CU were NEO-T+. BioFINDER-2's temporal-neocortex cut-off 1.36.
##   [N25] Ossenkoppele et al. 2025, Nat Neurosci (doi 10.1038/s41593-025-
##         02000-6): GMM cut-off for flortaucipir temporal cortex 1.40;
##         T+ among A+ CU 18.5% at 65, 20.9% at 75, 23.5% at 85.
##   [D18] Devous et al. 2018, J Nucl Med 59:937: test-retest SD of %
##         change, temporal 2.0% (PERSI), cerebellar references 3.7-4.2%.
##   [T20] Timmers et al. 2020, JCBFM 40:2464: cerebellar-grey SUVr
##         test-retest, lateral temporal 0.7-2.2%, global 1.6-2.8% (median).
## ---------------------------------------------------------------------
SETTINGS_TAU <- modifyList(SETTINGS, list(
  ## --- the tau level grid. NEO-T observed range in the cohort 0.92 - 4.01;
  ##     the grid runs well past both ends so that auxiliary tau coordinates
  ##     (subjects with no tau scans) are almost never cut off by the table.
  grid_lo    = 0.70,
  grid_hi    = 4.50,
  grid_step  = 0.002,
  ## --- positivity threshold: [O22]'s rule applied to [J20]'s numbers,
  ##     1.18 + 2 x 0.08 = 1.34, which is also BioFINDER-2's 1.36 and below
  ##     the GMM 1.40 of [N25]. The MODEL does not use this as a hard
  ##     classifier: it anchors the clock and truncates the flat levels.
  thresh     = 1.34,

  ## --- rate curve r_T(y). theta_T[1] is NOT pinned: its basis support
  ##     [0.70, 1.08] overlaps the plausible launch level, so unlike the
  ##     amyloid theta[1] it is (weakly) informed. Prior mean is the slow
  ##     regime of [J20], 0.005 SUVR/yr = exp(-5.3), with a wide sd.
  theta1     = -5.3,
  theta1_sd  = 2.0,
  init_rise  = 6,           # initial curve rises x6 from floor to threshold
                            # (0.005 -> 0.03 /yr, the two regimes of [J20])
  ## The adaptive random walk only where there are data to adapt to. Basis
  ## functions 7-10 have support above 2.54 SUVR, where 1.5% of tau scans
  ## sit; there the walk has a FIXED increment sd (as the noise curve does).
  ## With the adaptive prior the tail coefficients made excursions to
  ## log-rates of +0.7 and -12 (increment sd up to 30 allowed by the lambda
  ## truncation), which produced heavy tails and R-hat 1.06 for theta_T[7:8].
  rw_fixed_from = 7L,
  rw_tail_sd    = 0.5,

  ## --- measurement scale omega_T(y): a 3-4% test-retest SD of the
  ##     DIFFERENCE of two scans ([D18] cerebellar, [T20]) at ~1.2 SUVR is
  ##     0.04-0.05, i.e. 0.03 per scan; [J20]'s amyloid-negative scatter
  ##     implies ~0.04 per scan. The t4 scale is SD / 1.41, so 0.03 = exp(-3.5).
  nu1_mean   = -3.5,

  ## --- the onset layer for tau
  ## launch level = non-accumulator mean: [J20]'s amyloid-negative inferior
  ## temporal mean 1.18-1.19, less a few hundredths because the middle
  ## temporal gyrus runs lower; sd 0.05 covers pipeline/reference differences
  mu_b_mean  = 1.17, mu_b_sd = 0.05,   # truncated to (grid_lo, thresh)
  ## between-subject spread of flat levels: [J20]'s 0.08-0.09 includes scan
  ## noise, so the latent spread is ~0.07; half-normal scale 0.10, cap 0.3
  sigma_b_sd = 0.10, sigma_b_max = 0.30,
  ## conditional susceptibility logit q_i = v_i' phi_T. Cross-sectional
  ## NEO-T positivity among A+ CU is 15-24% ([O22], [N25]) and is a LOWER
  ## bound on the lifetime fraction, so the intercept is centred at
  ## q = 0.62 with a 95% prior range of 0.18-0.93
  phi1_mean  = 0.5,  phi1_sd = 1.0,
  ## tau departure: cortical tau follows amyloid by years ([J20]: increases
  ## confined to > 68 CL; [N25]: positivity rising through 65-85), so the
  ## AFT intercept is centred at age 75 rather than amyloid's 70
  g1_age     = 75
))
SETTINGS$init_rise   <- 20     # the amyloid init curve rises x20 (as before)
SETTINGS$g1_age      <- 70
SETTINGS$sigma_b_sd  <- 0.15; SETTINGS$sigma_b_max <- 0.5   # amyloid, as before

## Coupling and dependence (shared, not biomarker-specific).
SETTINGS_JOINT <- list(
  d_coup     = 0.10,        # PROVISIONAL sigmoid width of H_A(y), in amyloid SUVR
  beta_sd    = 1.0,         # beta ~ N(0, 1) truncated to beta >= 0
  rho_beta   = c(2, 2),     # (rho + 1)/2 ~ Beta(2, 2): either sign, mild pull to 0
  subj_block = "four",      # "four": one RW block per subject over
                            # (x_tilde_A, z_A, x_tilde_T, z_T) when the subject
                            # has tau scans; "pairs": two biomarker blocks
  ## sampler changes after the first full fit (run_joint_model.R, changes 1-3)
  amy_block  = TRUE,        # AF_slice block over mu_b, sigma_alpha, phi[1]
  tau_curve_split = 7L,     # theta_T[7:10] (basis support above ~2.5 SUVR, where
                            # <2% of tau scans are) sampled in their own block
  subj_scales = c(x_tilde = 0.05, z = 1, x_tilde_T = 0.15, z_T = 1)
                            # initial proposal sds of the subject RW blocks
)
