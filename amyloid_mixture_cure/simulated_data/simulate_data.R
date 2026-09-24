# OUTLINE FOR SIMULATING DATA

# Pipeline: truth (posterior medians) + design data (resample real data/distribution)
# -> generate latent quantities per subject (random)
# -> build true trajectories via the clock
# -> add measurement noise to get y_mat
# -> package data into same list format as adni_amyloid.rds
# -> save in rds files

###############################################
# TRUTH: posterior medians for population parameters
################################################
# Source: data/fitted_posterior.rds (use $g = 252 pooled draws x 33 params)
# Extract: phi[1:3], gamma[1:3], sigma_alpha, psi[1:2], sigma_delta, mu_b, sigma_b, theta[1:10], nu[1:8]
# Apply median to all these values
# Return in a truth list 

get_baseline_truth <- function(path="data/fitted_posterior.rds") {
  post <- readRDS(path)
  median <- apply(post$g, 2, median) # one median per parameter
  truth <- list(
    phi = med[c("phi[1]", "phi[2]", "phi[3]")],          # susceptibility
    gamma = med[c("gamma[1]", "gamma[2]", "gamma[3]")],  # onset age (Weibull AFT)
    sigma_alpha = med["sigma_alpha"],                     # onset-age spread
    psi = med[c("psi[1]", "psi[2]")],                     # rate covariate effects
    sigma_delta = med["sigma_delta"],                     # between-person speed spread
    mu_b = med["mu_b"],                                   # shared floor
    sigma_b = med["sigma_b"],                         # non-accumulator level spread
    theta = med[sprintf("theta[%d]", 1:10)],        # rate curve r_A(y)
    nu = med[sprintf("nu[%d]", 1:8)]             # measurement-scale curve omega_A(y)
  )
  # return median truth list
  }

#######################################
# DESIGN: Covariates + visited ages 
######################################
# Covariates drawn independently each subject:
#   APOE4 ~ Bernoulli(0.407) 
#   Female ~ Bernoulli(0.489)
# x centre = (apoe4 = 0.427, female = 0.503)
# X = APOE4 - x centre, female - x centre (add centering)
# Baseline age ~ Normal(72.5, 7.4), truncated [51, 94] - carriers are 1.7 years younger
# num subjects = 1101
# num obs = 2 + K, K ~ Geom(0.48), truncated at 7 --> give mean of 3 scans (can control follow-up intensity by lowering p) 
# gap between scans ~ normal(2.05, 0.15) - the median is 2.0 years [0.8-11.3]
# visit ages = baseline age + cumulative sum of gaps for each subject
# simulated follow-up length (last scan age - first) = max ~14, median ~4

get_design <- function(N = 1101, x_centre, seed) {
  apoe <- rbinom(N, 1, 0.407)
  female <- rbinom(N, 1, 0.489)
  X <- cbin(apoe = apoe - x_centre["apoe"], female = female - x_centre["female"]) # center X
  base_age <- rnorm(N, 72.5, 7.4) # need truncation
  num_obs <- pmin(2 + rgeom(N, prob=0.48), 7)
  gaps <- 
  visit_ages <- c(base_age, base_age + cumsum(gaps))
  # write these visit ages into age_mat, leaving zeroes after
  # return list(n_subj, n_obs, age_mat, X, x_centre, ids)
}

#######################################
# CURVES implied by the truth 
######################################
# Rate curve + clock: bA <- make_rate_basis()
#                     r <- rate_grid(bA, truth$theta)
#                     G <- clock_grid(bA, r, SETTINGS$thresh)
# vA built from real data: make_var_basis(y_obs_real)
# omega: exp(vA$B %*% truth$nu)

truth_curves <- function(truth) {
  bA <- make_rate_basis()
  r <- rate_grid(bA, truth$theta)
  G <- clock_grid(bA, r, SETTINGS$thresh)
  vA <- # fill in
  omega <- exp(vA$B %*% truth$nu)
  
  c(bA, r, G, vA, omega)
  # return list(bA, r, G, vA, omega)
}

#########################################
# LATENT QUANTITIES per subject (random) 
#########################################
# Draw these for every subject:
#   susceptible: Z_i ~ Bernoulli
#   onset age: same formula
#   speed: delta_i
#   non-acc level b_i ~ N(mu_b, sigma_b) truncated above at thresh
# For accumulator:
#   g_ij = G(mu_b) + exp(delta_i) * (age_ij - alpha_i)
#   g_ij = max(g_ij, G(mu_b))            # flat at mu_b before onset
#   A_ij = G^{-1}(g_ij)
# For non accmulator: A_ij = b_i at every visit
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
  x_init <- sapply(seq_len(design$N), function(i) mean(y_mat[i:design$n_obs[i]]))
  a_tilde = sapply(seq_len(design$N), function(i) mean(design$age_mat[i:design$n_obs[i]]))
  
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
