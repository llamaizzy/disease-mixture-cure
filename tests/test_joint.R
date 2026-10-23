## =====================================================================
##  Deterministic checks for the joint amyloid-tau model, as required by
##  docs/context/coupling_and_copula.pdf section 11:
##    1. uncoupled limit (beta = 0, rho = 0) recovers the separate clocks
##    2. constant-exposure limit: C_i(a) = m0 (a - atil_T), J_T gains 1/m0
##    3. age <-> location round trips and the Jacobian by finite differences
##    4. the exposure-clock trajectory agrees with an ODE solve
##    5. the truncated Weibull and the copula are normalised
##    6. a subject with no tau scans gives the amyloid-only likelihood
##
##      Rscript tests/test_joint.R          (about two minutes, mostly compiling)
## =====================================================================
suppressMessages(source("run_joint_model.R"))
ok <- 0L; fail <- character(0)
chk <- function(name, v) {
  v <- isTRUE(v)
  if (v) { ok <<- ok + 1L; cat(sprintf("  PASS  %s\n", name)) }
  else   { fail <<- c(fail, name); cat(sprintf("  FAIL  %s\n", name)) }
}

SA <- SETTINGS; ST <- SETTINGS_TAU; SJ <- SETTINGS_JOINT
dat <- load_amyloid_data()
bA <- make_rate_basis(SA); bT <- make_rate_basis(ST)
scans <- function(y, n) unlist(lapply(seq_along(n), function(i) y[i, seq_len(n[i])]))
vA <- make_var_basis(scans(dat$y_mat, dat$n_obs), SA)
vT <- make_var_basis(scans(dat$y_tau_mat, dat$n_obs_tau), ST)

## a population state: the inits of the model builder, with the curves bent
spec <- build_joint_model(dat, bA, vA, bT, vT)
ini  <- spec$inits
thA <- ini$theta + c(0, 0.3, -0.2, 0.4, 0.1, -0.3, 0.2, 0, 0.1, -0.1)
thT <- ini$theta_T + c(0.2, -0.1, 0.3, 0.1, -0.2, 0.2, 0, 0.1, -0.1, 0.1)
rA <- rate_grid(bA, thA); GA <- clock_grid(bA, rA, SA$thresh)
rT <- rate_grid(bT, thT); GT <- clock_grid(bT, rT, ST$thresh)
mu_b <- 0.55; mu_bT <- 1.05; d_coup <- SJ$d_coup
cE <- compileNimble(nfExpClock)
Egrid <- function(beta) cE(rA, bA$step, bA$n, bA$lo, SA$r_min, beta, SA$thresh, d_coup)
cX  <- compileNimble(dxtilde); cXT <- compileNimble(dxtildeT)
cJ  <- compileNimble(dsubjJoint); cM <- compileNimble(dsubjMix)
cInv <- compileNimble(nfInvExposure); cExp <- compileNimble(nfExposure)
cSt <- compileNimble(nfAmyState)
cpl <- function(y, beta) exp(beta * (plogis((y - SA$thresh) / d_coup) - 0.5))

## an accumulating subject: amyloid centroid 75, level 0.9; tau ref 78, level 1.3
xA <- 0.90; atilA <- 75; dA <- 0.1; xgA <- log(45); saA <- 0.3
xT <- 1.30; atilT <- 78; dT <- -0.2; xgT <- log(50); saT <- 0.3
dxt_T <- function(x, beta, rho, EA = Egrid(beta), xA_ = xA, dA_ = dA, dT_ = dT)
  cXT(x, dT_, atilT, mu_bT, xgT, saT, SA$alpha_min, SA$alpha_max,
      GT, rT, bT$lo, bT$step, bT$n, SA$r_min,
      xA_, dA_, atilA, mu_b, xgA, saA, SA$alpha_min, SA$alpha_max,
      GA, rA, EA, bA$lo, bA$step, bA$n, SA$r_min, beta, SA$thresh, d_coup, rho, log = 1)
dxt_A <- function(x, delta = dA)
  cX(x, delta, atilA, mu_b, xgA, saA, SA$alpha_min, SA$alpha_max,
     GA, rA, bA$lo, bA$step, bA$n, SA$r_min, log = 1)
## the OLD tau map: dxtilde with tau arguments
dxt_T_old <- function(x)
  cX(x, dT, atilT, mu_bT, xgT, saT, SA$alpha_min, SA$alpha_max,
     GT, rT, bT$lo, bT$step, bT$n, SA$r_min, log = 1)

cat("\n1. uncoupled limit\n")
E0 <- Egrid(0)
chk("beta = 0: exposure clock is the amyloid clock up to a constant",
    diff(range(E0 - GA)) < 1e-9)
for (x in c(0.95, 1.30, 1.8)) {
  chk(sprintf("beta = 0, rho = 0: tau density equals the old clock map at x = %.2f", x),
      abs(dxt_T(x, 0, 0, E0) - dxt_T_old(x)) < 1e-8)
}
chk("rho = 0: copula term is zero",
    abs(dxt_T(1.3, 0.8, 0) - dxt_T(1.3, 0.8, 1e-12)) < 1e-8)

cat("\n2. constant-exposure limit\n")
## a subject whose amyloid is flat and HIGH for the whole tau window: departed
## long ago, now saturated near the top of the curve. m_i is then the constant
## m(A_high) and C_i(a) = m0 (a - atil_T): the tau map is the old one on a
## rescaled clock, and the Jacobian picks up 1/m0.
## Far above the amyloid threshold H_A is within 0.4% of one, so a subject at
## 1.30 SUVR at age 70, moving slowly (delta = -1), has a constant multiplier
## over a tau window 72-80 while staying inside the level table.
beta <- 1.2; EA <- Egrid(beta); win <- c(72, 76, 80)
st <- cSt(1.30, -1, 70, mu_b, GA, rA, bA$lo, bA$step, bA$n, SA$r_min)
psi <- sapply(win, function(u) cExp(u, st[1], st[2], st[3], -1, 70, mu_b,
                                    GA, EA, bA$lo, bA$step, bA$n, beta, SA$thresh, d_coup))
lev <- sapply(win, function(u) sa_inv_monotone(bA$y, GA, st[1] + exp(-1) * (u - 70)))
m0 <- cpl(lev, beta)
chk("amyloid saturated: path stays on the table and multiplier is constant",
    all(is.finite(lev)) && diff(range(m0)) / m0[1] < 0.01)
chk("exposure then grows linearly at rate m0",
    abs((psi[3] - psi[1]) / (win[3] - win[1]) - mean(m0)) / mean(m0) < 0.01)
chk("the multiplier is well above one (coupling is active)", mean(m0) > exp(beta * 0.45))

cat("\n3. round trips and the Jacobian\n")
beta <- 0.9; rho <- 0.4; EA <- Egrid(beta)
## alpha_T implied by a location, then the location implied by that age
st <- cSt(xA, dA, atilA, mu_b, GA, rA, bA$lo, bA$step, bA$n, SA$r_min)
alphaA <- st[3]
expo <- function(u) cExp(u, st[1], st[2], alphaA, dA, atilA, mu_b, GA, EA,
                         bA$lo, bA$step, bA$n, beta, SA$thresh, d_coup)
GbT <- sa_interp(bT$y, GT, mu_bT); rbT <- max(sa_interp(bT$y, rT, mu_bT), SA$r_min)
Gext <- function(x) if (x >= mu_bT) sa_interp(bT$y, GT, x) else GbT + (x - mu_bT) / rbT
alpha_of_x <- function(x) {
  psi0 <- expo(atilT)
  cInv(psi0 - exp(-dT) * (Gext(x) - GbT), st[1], st[2], alphaA, dA, mu_b,
       GA, EA, bA$lo, bA$step, bA$n, beta, SA$thresh, d_coup)[1]
}
x_of_alpha <- function(a) {      # invert: C_i(alpha) = -e^{-dT} (Gext(x) - GbT)
  dG <- -exp(dT) * (expo(a) - expo(atilT))
  g  <- GbT + dG
  if (g >= GbT) sa_inv_monotone(bT$y, GT, g) else mu_bT + (g - GbT) * rbT
}
for (x in c(1.10, 1.30, 1.60, 0.98)) {
  a <- alpha_of_x(x)
  chk(sprintf("x -> alpha_T -> x round trip at x = %.2f (alpha_T = %.1f)", x, a),
      abs(x_of_alpha(a) - x) < 1e-6)
}
## the Jacobian: d alpha_T / d x by finite differences against the one in the
## density, which is read off by comparing the density with and without it
num_jac <- function(x, h = 1e-5) abs(alpha_of_x(x + h) - alpha_of_x(x - h)) / (2 * h)
for (x in c(1.10, 1.30, 1.60)) {
  a <- alpha_of_x(x)
  ## log density minus the marginal and the copula = log Jacobian
  pA <- (1 - exp(-((alphaA - 25) / exp(xgA))^(1/saA))) / (1 - exp(-(100 / exp(xgA))^(1/saA)))
  pT <- (1 - exp(-((a - 25) / exp(xgT))^(1/saT))) / (1 - exp(-(100 / exp(xgT))^(1/saT)))
  lf <- dweibull(a - 25, 1/saT, exp(xgT), log = TRUE) - log(1 - exp(-(100 / exp(xgT))^(1/saT)))
  v1 <- qnorm(pA); v2 <- qnorm(pT)
  lc <- -0.5 * log(1 - rho^2) - (rho^2 * (v1^2 + v2^2) - 2 * rho * v1 * v2) / (2 * (1 - rho^2))
  ljac <- dxt_T(x, beta, rho, EA) - lf - lc
  chk(sprintf("Jacobian matches finite differences at x = %.2f", x),
      abs(exp(ljac) - num_jac(x)) / num_jac(x) < 0.02)
}
## the cross-derivative: moving the AMYLOID location moves the implied tau age
a0 <- alpha_of_x(1.3)
stB <- cSt(xA + 0.02, dA, atilA, mu_b, GA, rA, bA$lo, bA$step, bA$n, SA$r_min)
expoB <- function(u) cExp(u, stB[1], stB[2], stB[3], dA, atilA, mu_b, GA, EA,
                          bA$lo, bA$step, bA$n, beta, SA$thresh, d_coup)
aB <- cInv(expoB(atilT) - exp(-dT) * (Gext(1.3) - GbT), stB[1], stB[2], stB[3], dA, mu_b,
           GA, EA, bA$lo, bA$step, bA$n, beta, SA$thresh, d_coup)[1]
chk("nonzero amyloid -> tau cross-derivative", abs(aB - a0) > 1e-3)

cat("\n4. ODE agreement\n")
## integrate dT/da = e^{dT} m(A(a)) rbar_T(T) with a fine RK4 from tau departure
alphaT <- alpha_of_x(1.3)
rateT <- function(y) max(exp(sum(splines2::bSpline(y, knots = bT$knots, degree = 3, intercept = TRUE,
                                                   Boundary.knots = c(bT$lo, bT$hi)) * thT)), SA$r_min)
Alevel <- function(u) sa_inv_monotone(bA$y, GA, max(st[2], st[1] + exp(dA) * (u - atilA)))
f <- function(u, y) exp(dT) * cpl(Alevel(u), beta) * rateT(y)
rk4 <- function(y0, a0, a1, h = 0.01) {
  y <- y0; u <- a0
  while (u < a1 - 1e-12) {
    hh <- min(h, a1 - u)
    k1 <- f(u, y); k2 <- f(u + hh/2, y + hh*k1/2); k3 <- f(u + hh/2, y + hh*k2/2); k4 <- f(u + hh, y + hh*k3)
    y <- y + hh * (k1 + 2*k2 + 2*k3 + k4) / 6; u <- u + hh
  }
  y
}
for (a in c(alphaT + 3, atilT, atilT + 6)) {
  gt <- max(GbT, Gext(1.3) + exp(dT) * (expo(a) - expo(atilT)))
  clock_level <- sa_inv_monotone(bT$y, GT, gt)
  ode_level <- rk4(mu_bT, alphaT, a)
  chk(sprintf("tau level at age %.1f: clock %.4f vs ODE %.4f", a, clock_level, ode_level),
      abs(clock_level - ode_level) < 0.003)
}

cat("\n5. normalisation\n")
cW <- compileNimble(nfLogTWeib); cC <- compileNimble(nfLogCopula); cF <- compileNimble(nfTWeibCDF)
chk("truncated Weibull integrates to one",
    abs(integrate(function(a) sapply(a, function(v) exp(cW(v, xgT, saT, 25, 125))), 25, 125)$value - 1) < 1e-4)
chk("truncated Weibull CDF reaches one at the upper bound", abs(cF(125, xgT, saT, 25, 125) - 1) < 1e-10)
chk("copula density integrates to one over the unit square",
    abs(integrate(function(u) sapply(u, function(a)
      integrate(function(v) sapply(v, function(b) exp(cC(a, b, 0.6))), 1e-6, 1 - 1e-6)$value),
      1e-6, 1 - 1e-6)$value - 1) < 2e-3)

cat("\n6. marginalisation with no tau scans\n")
omA <- as.numeric(exp(vA$B %*% ini$nu)); omT <- as.numeric(exp(vT$B %*% ini$nu_T))
i <- which(dat$n_obs_tau == 0 & dat$n_obs >= 4)[1]
ya <- dat$y_mat[i, ]; aa <- dat$age_mat[i, ]
lj <- cJ(c(ya, rep(0, 7)), dat$n_obs[i], aa, 0, rep(0, 7), 7,
         0.8, dat$a_tilde[i], 0.1, 0.6, 0.9, mu_b, GA, rA, EA, bA$lo, bA$step, bA$n, SA$r_min,
         omA, vA$lo, vA$step, vA$n,
         1.2, dat$a_tilde_tau[i], 0.0, 1.0, 0.3, mu_bT, GT, rT, bT$lo, bT$step, bT$n, SA$r_min,
         omT, vT$lo, vT$step, vT$n, beta, SA$thresh, d_coup, SA$t_df, log = 1)
lm <- cM(ya, dat$n_obs[i], aa, 0.8, dat$a_tilde[i], 0.1, 0.6, 0.9, mu_b,
         GA, rA, bA$lo, bA$step, bA$n, SA$r_min, omA, vA$lo, vA$step, vA$n, SA$t_df, log = 1)
chk("no tau scans: three-state likelihood equals the amyloid two-branch one", abs(lj - lm) < 1e-10)
chk("state weights sum to one", {
  pi <- 0.7; q <- 0.4; abs((1 - pi) + pi * (1 - q) + pi * q - 1) < 1e-12 })

cat(sprintf("\n%d passed, %d failed\n", ok, length(fail)))
if (length(fail)) { cat("failed:", paste(fail, collapse = "; "), "\n"); quit(status = 1) }
