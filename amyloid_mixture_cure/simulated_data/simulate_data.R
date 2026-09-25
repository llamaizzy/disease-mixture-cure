# OUTLINE FOR GENERATING SIMULATION DATA

# Pipeline: truth (posterior medians) + design data (resample real data/distribution)
# -> generate latent quantities per subject (random)
# -> build true trajectories via the clock
# -> add measurement noise to get y_mat
# -> package data into same list format as adni_amyloid.rds
# -> save in rds files

###############################################
# 1. TRUTH: posterior medians for population parameters
################################################
# Source: data/fitted_posterior.rds (use $g = 252 pooled draws x 33 params)
# Extract: phi[1:3], gamma[1:3], sigma_alpha, psi[1:2], sigma_delta, mu_b, sigma_b, theta[1:10], nu[1:8]
# Apply median to all these values
# Return in a truth list 

get_baseline_truth <- function(path="data/fitted_posterior.rds") {
  post <- readRDS(path)
  median <- apply(post$g, 2, median) # one median per parameter
  list(
    phi = med[c("phi[1]", "phi[2]", "phi[3]")],          # susceptibility
    gamma = med[c("gamma[1]", "gamma[2]", "gamma[3]")],  # onset age (Weibull AFT)
    sigma_alpha = med["sigma_alpha"],                     # onset-age spread
    psi = med[c("psi[1]", "psi[2]")],                     # rate covariate effects
    sigma_delta = med["sigma_delta"],                     # between-person speed spread
    mu_b = med["mu_b"],                                   # shared floor
    sigma_b = med["sigma_b"],                         # non-accumulator level spread
    theta = med[sprintf("theta[%d]", 1:10)],        # rate curve coeffs
    nu = med[sprintf("nu[%d]", 1:8)]             # measurement-scale curve omega_A(y)
  )
  # return median truth list
}

########################################
# 2. SCENARIOS 
########################################
# Each scenarios is baseline truth with one modification
#   baseline: no change
#   onset spread lo/hi: sigma_alpha * 0.5 or * 2
#   speed spread lo/hi: sigma_delta * 0.5 or * 2
#   measure error lo/hi: omega * 0.5 or * 2
# Check every scenario stays inside the prior supports in SETTINGS

make_scenarios <- function(truth) {
  list(
    baseline = truth,
    onset_spread_lo = ,
    onset_spread_hi = ,
    speed_spread_lo = ,
    speed_spread_hi = ,
    meas_error_lo =  ,
    meas_error_hi = 
  )
  # return name list of modified truth lists
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
# gap between scans (n_obs - 1 gaps): two-part mixture with 0.76 probability of normal(2.05, 0.15), otherwise ~2.5 + exp(mean=1)
#   old gap = pmax(rnorm(Ji - 1, visit_gap_mean, visit_gap_sd), min_gap)
# visit ages = baseline age + cumulative sum of gaps
# write visit ages into age_mat, padding with 0s after
# CHECK: simulated follow-up length (last scan age - first) = max ~14, median ~4

get_design <- function(N = 1101, p_apoe4 = 0.407, p_female = 0.489, 
                       x_centre = c(apoe4 = 0.427, female = 0.503),
                       age_mean = 72.5, age_sd = 7.4, age_range = c(50, 95),
                       visits_range = 2:7,
                       visit_prob = c(487, 296, 155, 83, 64, 16) / 1101,
                       p_ontime = 0.76, seed) {
  ## covariates: 0/1, then centred like the real data
  apoe4  <- rbinom(N, 1, p_apoe4)
  female <- rbinom(N, 1, p_female)
  X <- cbin(apoe = apoe4 - x_centre["apoe"], female = female - x_centre["female"])
  
  # baseline age
  base_age <- rnorm(N, age_mean, age_sd) # add truncation?
  num_obs <- 
  max_n <- max(visits_range)
 
   # visit ages: on-schedule / delayed mixture gaps with 0-padding
  age_mat <- matrix(0, N, max_n)
  gaps <- 
  visit_ages <- c(base_age, base_age + cumsum(gaps))
  
  # write these visit ages into age_mat, leaving zeroes after

  # return list(n_subj, n_obs, age_mat, X, x_centre, apoe4, female)
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
  omega <- exp(vA$B %*% truth$nu)
  list(bA, r, G, vA, omega)
  # return list(bA, r, G, vA, omega)
}

#########################################
# LATENT QUANTITIES per subject (random) 
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
# Flag any A_ij beyond grid_hi (return -Inf)

draw_latent <- function(truth, design, curves) {
  # return data frame(Z, alpha, delta, z, b, A_ij)
}

#########################################
# Add measurement noise to get y 
#########################################
# y_ij = A_ij + omega(A_ij) * t_ij,   t_ij ~ Student-t(df = t_df = 4)
# omega is a scale
get_y <- function(latent, design, curves) {
  # return y_mat
}

#########################################
# PACKAGE into real data format
#########################################
# Match field for what build_model() reads:
# n_subj, n_obs, y_mat, age_mat, a_tilde, x_init, X, ids, x_centre, thresh
#   a_tilde = mean visit age per subject
#   x_init  = mean simulated y per subject
#unused cells in y_mat / age_mat = 0

get_dat <- function(y_mat, design) {
  x_init <- sapply(seq_len(design$N), function(i) mean(y_mat[i, 1:design$n_obs[i]]))
  a_tilde = sapply(seq_len(design$N), function(i) mean(design$age_mat[i, 1:design$n_obs[i]]))
  # return in fixed format
}


#########################################
# CHECK dataset for validity
#########################################
# no NA/non-finite values; all y within the SUVR grid
# share of accumulators ~ pi_bar; share above thresh ~ real data
# onset ages within (alpha_min, alpha_max]
# spaghetti plot vs the real data (check before fitting)
check_dat <- function(dat, latent) {
  # implement
}

#########################################
# Generate one dataset
#########################################
# Check: plot this against the real data before run simulation
simulate_one <- function(truth, N = NULL, seed) {
# set seed
# design <- get_design(N, x_centre, seed)
# curves <- truth_curves(truth)
# latent <- draw_latent(truth, design, curves)
# y_mat  <- get_y(latent, design, curves)
# dat    <- get_dat(y_mat, design)
# check_dat(dat, latent)
# return  list(dat = dat, truth = truth, latent = latent, seed = seed)
}

#########################################
# Simulate
#########################################
# Generate different scenarios and replicates