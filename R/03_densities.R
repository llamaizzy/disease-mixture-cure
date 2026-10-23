## =====================================================================
##  The custom NIMBLE densities.
##
##  Two things are custom here, and both are custom for the same reason:
##  a subject's contribution is not a product of independent per-visit
##  terms that BUGS syntax can express one at a time. It is a single
##  function of a trajectory, and for the mixture it is a LOG-SUM-EXP
##  over two trajectories.
##
##  A rule that runs through all of it: an infeasible state returns
##  -Inf, never a large finite value. A large finite value is a plateau,
##  and a slice sampler will happily walk out onto a plateau and stay
##  there.
## =====================================================================
suppressPackageStartupMessages(library(nimble))

## ---- numerical primitives -------------------------------------------

## Cumulative trapezoid of f on a uniform grid, starting at 0.
nfCumTrapz <- nimbleFunction(
  run = function(f = double(1), step = double(0), n = double(0)) {
    returnType(double(1))
    out <- nimNumeric(n); out[1] <- 0
    for (i in 2:n) out[i] <- out[i-1] + step * (f[i] + f[i-1]) / 2
    return(out)
  })

## Linear interpolation on a uniform grid. This CLAMPS out of range, so
## every caller bounds-checks first and returns -Inf itself.
nfInterp <- nimbleFunction(
  run = function(xq = double(0), lo = double(0), step = double(0),
                 n = double(0), vals = double(1)) {
    returnType(double(0))
    pos <- (xq - lo) / step + 1
    if (pos <= 1) return(vals[1])
    if (pos >= n) return(vals[n])
    i0 <- floor(pos); w <- pos - i0
    return(vals[i0] * (1 - w) + vals[i0 + 1] * w)
  })

## Inverse of a strictly increasing tabulated G, by bisection. This is how
## a clock time is turned back into a level.
nfInvMono <- nimbleFunction(
  run = function(g = double(0), lo = double(0), step = double(0),
                 n = double(0), G = double(1)) {
    returnType(double(0))
    if (g <= G[1]) return(lo)
    if (g >= G[n]) return(lo + (n - 1) * step)
    a <- 1; b <- n
    while (b - a > 1) { m <- floor((a + b) / 2); if (G[m] <= g) a <- m else b <- m }
    dG <- G[b] - G[a]; w <- 0
    if (dG > 0) w <- (g - G[a]) / dG
    return(lo + (a - 1 + w) * step)
  })

## log density of a scaled Student-t. omega is a SCALE, not an SD: the
## implied SD is omega * sqrt(df / (df - 2)).
nfLogDT <- nimbleFunction(
  run = function(y = double(0), mu = double(0), omega = double(0), df = double(0)) {
    returnType(double(0))
    z <- (y - mu) / omega
    return(lgamma((df + 1) / 2) - lgamma(df / 2) - 0.5 * log(df * pi) -
             log(omega) - (df + 1) / 2 * log(1 + z * z / df))
  })

## THE CLOCK:  G(y) = int_{thresh}^{y} du / max(r(u), r_min),  tabulated.
## Anchored at the positivity threshold, so G = 0 there and disease time is
## measured relative to crossing.
nfClock <- nimbleFunction(
  run = function(r = double(1), step = double(0), n = double(0),
                 lo = double(0), thresh = double(0), r_min = double(0)) {
    returnType(double(1))
    f <- nimNumeric(n)
    for (i in 1:n) f[i] <- 1 / max(r[i], r_min)
    ## The local must be TYPED before taking a nested nimbleFunction's return
    ## value, or NIMBLE infers a matrix map and the by-reference vector
    ## argument of nfInterp will not bind.
    G <- nimNumeric(n)
    G[1:n] <- nfCumTrapz(f, step, n)
    g0 <- nfInterp(thresh, lo, step, n, G)
    out <- nimNumeric(n)
    for (i in 1:n) out[i] <- G[i] - g0
    return(out)
  })

## ---- the implied density of x_tilde ---------------------------------
##
##  x_tilde_i is the accumulator's level at their OWN visit centroid. It has
##  no free prior. The departure age is a deterministic function of it,
##
##      alpha_i = atil_i - exp(-delta_i) { G_ext(x_i) - G_ext(mu_b) },
##
##  so the AFT prior on alpha plus a change of variables FIXES the density
##  of x_tilde. Putting a separate prior on x_tilde as well would either
##  double-count or displace the AFT, leaving gamma unidentified.
##
##  WHY THE CLOCK IS EXTENDED BELOW mu_b. Restricting x >= mu_b would put a
##  point mass at exactly mu_b, of mass S_alpha(atil_i) -- the subjects who
##  are susceptible but have not departed yet. A slice sampler cannot reach
##  an atom, so the model would silently rule that state out and push every
##  flat subject into the cured branch, inflating the cure fraction BY
##  CONSTRUCTION. Continuing the clock linearly below the floor,
##
##      G_ext(y) = G(mu_b) + (y - mu_b) / rbar(mu_b)   for y < mu_b,
##
##  spreads that atom into a density. Below mu_b, x is ONLY a coordinate
##  meaning "departs after the centroid"; no level below the floor is implied.
dxtilde <- nimbleFunction(
  run = function(x = double(0), delta = double(0), atil = double(0),
                 mu_b = double(0), xg = double(0), sigma_alpha = double(0),
                 alpha_min = double(0), alpha_max = double(0),
                 G = double(1), rgrid = double(1), lo = double(0),
                 step = double(0), ngrid = double(0), r_min = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    hi <- lo + (ngrid - 1) * step
    Gb <- nfInterp(mu_b, lo, step, ngrid, G)
    rb <- max(nfInterp(mu_b, lo, step, ngrid, rgrid), r_min)
    clock <- 0; ljac <- 0
    if (x >= mu_b) {
      if (x > hi) { if (log) return(-Inf) else return(0) }
      clock <- nfInterp(x, lo, step, ngrid, G) - Gb
      ljac  <- -log(max(nfInterp(x, lo, step, ngrid, rgrid), r_min))
    } else {
      clock <- (x - mu_b) / rb
      ljac  <- -log(rb)
    }
    alpha <- atil - exp(-delta) * clock
    if (alpha <= alpha_min | alpha > alpha_max) { if (log) return(-Inf) else return(0) }
    ## truncated Weibull, AFT parameterisation: kappa = 1/sigma_alpha shape,
    ## lambda = exp(x'gamma) scale, so gamma reads as a log TIME RATIO.
    kappa <- 1 / sigma_alpha
    lam   <- exp(xg)
    t     <- alpha - alpha_min
    lf   <- log(kappa) - log(lam) + (kappa - 1) * log(t / lam) - (t / lam)^kappa
    lnrm <- log(1 - exp(-((alpha_max - alpha_min) / lam)^kappa))
    ld <- lf - lnrm - delta + ljac      # -delta and ljac are the Jacobian
    if (is.nan(ld)) { if (log) return(-Inf) else return(0) }
    if (log) return(ld) else return(exp(ld))
  })

rxtilde <- nimbleFunction(
  run = function(n = integer(0), delta = double(0), atil = double(0),
                 mu_b = double(0), xg = double(0), sigma_alpha = double(0),
                 alpha_min = double(0), alpha_max = double(0),
                 G = double(1), rgrid = double(1), lo = double(0),
                 step = double(0), ngrid = double(0), r_min = double(0)) {
    returnType(double(0))
    lam  <- exp(xg)
    Fmax <- 1 - exp(-((alpha_max - alpha_min) / lam)^(1 / sigma_alpha))
    u     <- runif(1, 0, Fmax)
    alpha <- alpha_min + lam * (-log(1 - u))^sigma_alpha
    Gb <- nfInterp(mu_b, lo, step, ngrid, G)
    rb <- max(nfInterp(mu_b, lo, step, ngrid, rgrid), r_min)
    clock <- exp(delta) * (atil - alpha)
    if (clock >= 0) {
      if (Gb + clock >= G[ngrid]) return(lo + (ngrid - 1) * step + 1)   # off grid
      return(nfInvMono(Gb + clock, lo, step, ngrid, G))
    }
    return(mu_b + rb * clock)
  })

## ---- the two-branch subject likelihood -------------------------------
##
##  Branch S (accumulator): substituting alpha into the clamped trajectory
##  makes alpha CANCEL,
##
##     G(A_i(a)) = G(mu_b) + e^delta (a - alpha_i)
##               = G_ext(x_tilde_i) + e^delta (a - atil_i),
##
##  so "flat before departure" is one max() on the clock, not a new
##  integrator.  Branch N (non-accumulator): flat at b_i for life.
##
##  The two are combined by LOG-SUM-EXP. The membership indicator is never
##  sampled: marginalising it is both cheaper and better mixing than a
##  discrete latent variable, and it is what makes the per-subject
##  posterior probability of susceptibility available as a by-product.
##
##  Both branches use the SAME measurement-scale function omega. Letting
##  them differ would let the mixture explain a flat subject by a smaller
##  error rather than by a different mechanism.
dsubjMix <- nimbleFunction(
  run = function(x = double(1), n_obs = double(0), ages = double(1),
                 x_tilde = double(0), a_tilde = double(0), delta = double(0),
                 b = double(0), logit_pi = double(0), mu_b = double(0),
                 G = double(1), rgrid = double(1), lo = double(0),
                 step = double(0), ngrid = double(0), r_min = double(0),
                 omega = double(1), vlo = double(0), vstep = double(0),
                 nvgrid = double(0), df = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    hi <- lo + (ngrid - 1) * step
    if (x_tilde >= hi) { if (log) return(-Inf) else return(0) }
    Gb <- nfInterp(mu_b, lo, step, ngrid, G)
    rb <- max(nfInterp(mu_b, lo, step, ngrid, rgrid), r_min)
    g0 <- 0
    if (x_tilde >= mu_b) g0 <- nfInterp(x_tilde, lo, step, ngrid, G)
    else                 g0 <- Gb + (x_tilde - mu_b) / rb
    ed <- exp(delta)
    lS <- 0; lN <- 0
    omN <- nfInterp(b, vlo, vstep, nvgrid, omega)
    for (j in 1:n_obs) {
      gt <- g0 + ed * (ages[j] - a_tilde)
      if (gt < Gb) gt <- Gb                     # flat at mu_b before departure
      if (gt > G[ngrid]) { if (log) return(-Inf) else return(0) }
      mu <- nfInvMono(gt, lo, step, ngrid, G)
      om <- nfInterp(mu, vlo, vstep, nvgrid, omega)
      lS <- lS + nfLogDT(x[j], mu, om, df)
      lN <- lN + nfLogDT(x[j], b,  omN, df)
    }
    lp <- logit_pi
    aS <- -log(1 + exp(-lp)) + lS               # log(pi)     + loglik_S
    aN <- -log(1 + exp(lp))  + lN               # log(1 - pi) + loglik_N
    mx <- max(aS, aN)
    ll <- mx + log(1 + exp(-abs(aS - aN)))      # logsumexp, stable form
    if (is.nan(ll)) { if (log) return(-Inf) else return(0) }
    if (ll == Inf)  { if (log) return(-Inf) else return(0) }
    if (log) return(ll) else return(exp(ll))
  })

rsubjMix <- nimbleFunction(
  run = function(n = integer(0), n_obs = double(0), ages = double(1),
                 x_tilde = double(0), a_tilde = double(0), delta = double(0),
                 b = double(0), logit_pi = double(0), mu_b = double(0),
                 G = double(1), rgrid = double(1), lo = double(0),
                 step = double(0), ngrid = double(0), r_min = double(0),
                 omega = double(1), vlo = double(0), vstep = double(0),
                 nvgrid = double(0), df = double(0)) {
    returnType(double(1)); return(nimNumeric(length(ages)))
  })

registerDistributions(list(
  dxtilde = list(
    BUGSdist = "dxtilde(delta, atil, mu_b, xg, sigma_alpha, alpha_min, alpha_max, G, rgrid, lo, step, ngrid, r_min)",
    types = c("G = double(1)", "rgrid = double(1)"), pqAvail = FALSE),
  dsubjMix = list(
    BUGSdist = "dsubjMix(n_obs, ages, x_tilde, a_tilde, delta, b, logit_pi, mu_b, G, rgrid, lo, step, ngrid, r_min, omega, vlo, vstep, nvgrid, df)",
    types = c("value = double(1)", "ages = double(1)", "G = double(1)",
              "rgrid = double(1)", "omega = double(1)"), pqAvail = FALSE)))
