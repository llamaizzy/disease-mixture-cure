## =====================================================================
##  Custom NIMBLE densities for the JOINT amyloid-tau model
##  (docs/context/coupling_and_copula.pdf). Builds on the primitives in
##  R/03_densities.R: nfInterp, nfInvMono, nfLogDT, nfClock, dxtilde.
##
##  WHAT IS NEW
##
##   1. Tau has its own clock G_T, launch level mu_bT and threshold c_T,
##      and its progression is MULTIPLIED by the subject's current amyloid
##      burden:   dT/da = e^{delta_T} m_i(a) rbar_T(T),
##      m_i(a) = exp{ beta [H_A(A_i(a)) - 1/2] },  H_A = expit((y - c_A)/d).
##
##   2. THE EXPOSURE CLOCK. The note integrates m_i over age per subject.
##      That is unnecessary: m_i depends on age only through the amyloid
##      LEVEL, and on the accumulating piece du = e^{-delta_A} dy / rbar_A(y),
##      so
##          int_{alpha_A}^{a} m_i(u) du = e^{-delta_A} [ E(A_i(a)) - E(mu_bA) ],
##          E(y) = int_{lo}^{y} m(v) / rbar_A(v) dv,
##      a POPULATION table like G_A (it depends on theta_A and beta only).
##      Before amyloid departure m_i is the constant m(mu_bA). A subject's
##      exposure Psi_i(u) (exposure accumulated since amyloid departure) is
##      therefore two lookups, and its inverse is a monotone table inverse
##      followed by the amyloid clock -- no quadrature, no root finding.
##      At beta = 0, E = G_A + const and Psi_i(u) = u - alpha_A exactly, so
##      the uncoupled model is recovered to rounding.
##
##   3. The tau location coordinate x_tilde_T maps to a departure age
##      through the exposure clock, not calendar time:
##          C_i(alpha_T) = -e^{-delta_T} [G_T,ext(x_T) - G_T(mu_bT)],
##          C_i(a) = Psi_i(a) - Psi_i(atil_T),
##      and its Jacobian gains the factor 1 / m_i(alpha_T)  (note, eq. 20).
##
##   4. One Gaussian copula on the two departure ages, added ONCE, in the
##      density of x_tilde_T (which is the node that sees both ages).
##
##   5. The subject likelihood is a log-sum-exp over THREE states:
##      00 (neither), 10 (amyloid only), 11 (both), with Z_T <= Z_A.
##
##  CONVENTIONS KEPT FROM THE AMYLOID MODEL
##   * infeasible -> -Inf, never a plateau.
##   * off the end of a level table -> -Inf (no tail continuation yet; see
##     the note's "TEST: numerical and tail assumptions"). The tau grid is
##     deliberately wide so this is rare for auxiliary coordinates.
##   * the Jacobian uses 1/rbar at the point, as dxtilde does.
## =====================================================================

## exp(x) - 1 without cancellation near zero (expm1 is not in the DSL)
nfExpm1 <- nimbleFunction(
  run = function(x = double(0)) {
    returnType(double(0))
    if (abs(x) < 1e-5) return(x + x * x / 2)
    return(exp(x) - 1)
  })

## The coupling multiplier m(y) = exp{ beta [expit((y - cA)/d) - 1/2] }.
nfCouple <- nimbleFunction(
  run = function(y = double(0), beta = double(0), cA = double(0), d = double(0)) {
    returnType(double(0))
    return(exp(beta * (1 / (1 + exp(-(y - cA) / d)) - 0.5)))
  })

## THE EXPOSURE CLOCK  E(y) = int_{lo}^{y} m(v) / max(r(v), r_min) dv.
nfExpClock <- nimbleFunction(
  run = function(r = double(1), step = double(0), n = double(0), lo = double(0),
                 r_min = double(0), beta = double(0), cA = double(0), d = double(0)) {
    returnType(double(1))
    f <- nimNumeric(n)
    for (i in 1:n) f[i] <- nfCouple(lo + (i - 1) * step, beta, cA, d) / max(r[i], r_min)
    E <- nimNumeric(n)
    E[1:n] <- nfCumTrapz(f, step, n)
    return(E)
  })

## ---- truncated Weibull (AFT) pieces, shared by both biomarkers ----------
## shape kappa = 1/sigma_alpha, scale lam = exp(x'gamma), support (amin, amax].

## log density
nfLogTWeib <- nimbleFunction(
  run = function(alpha = double(0), xg = double(0), sigma_alpha = double(0),
                 amin = double(0), amax = double(0)) {
    returnType(double(0))
    kappa <- 1 / sigma_alpha; lam <- exp(xg); t <- alpha - amin
    lf   <- log(kappa) - log(lam) + (kappa - 1) * log(t / lam) - (t / lam)^kappa
    lnrm <- log(-nfExpm1(-((amax - amin) / lam)^kappa))
    return(lf - lnrm)
  })

## CDF, in (0, 1) on the interior of the support
nfTWeibCDF <- nimbleFunction(
  run = function(alpha = double(0), xg = double(0), sigma_alpha = double(0),
                 amin = double(0), amax = double(0)) {
    returnType(double(0))
    kappa <- 1 / sigma_alpha; lam <- exp(xg)
    num <- -nfExpm1(-((alpha - amin) / lam)^kappa)
    den <- -nfExpm1(-((amax - amin) / lam)^kappa)
    return(num / den)
  })

## log Gaussian copula density at percentiles (p1, p2) with correlation rho
nfLogCopula <- nimbleFunction(
  run = function(p1 = double(0), p2 = double(0), rho = double(0)) {
    returnType(double(0))
    v1 <- qnorm(p1, 0, 1); v2 <- qnorm(p2, 0, 1); r2 <- rho * rho
    return(-0.5 * log(1 - r2) -
             (r2 * (v1 * v1 + v2 * v2) - 2 * rho * v1 * v2) / (2 * (1 - r2)))
  })

## ---- the subject's amyloid path, in clock coordinates -------------------
##
##  Returns (g0, Gb, alpha_A): the clock value at the amyloid centroid, the
##  clock at the floor, and the implied amyloid departure age. The amyloid
##  level at age u is then  G^{-1}( max(Gb, g0 + e^delta (u - atil)) ).
##  Returns alpha_A = -Inf if x_tilde_A is off the top of the table.
nfAmyState <- nimbleFunction(
  run = function(xA = double(0), deltaA = double(0), atilA = double(0),
                 mu_bA = double(0), G = double(1), rgrid = double(1),
                 lo = double(0), step = double(0), ngrid = double(0),
                 r_min = double(0)) {
    returnType(double(1))
    out <- nimNumeric(3)
    hi <- lo + (ngrid - 1) * step
    Gb <- nfInterp(mu_bA, lo, step, ngrid, G)
    rb <- max(nfInterp(mu_bA, lo, step, ngrid, rgrid), r_min)
    g0 <- 0
    if (xA >= mu_bA) {
      if (xA > hi) { out[3] <- -Inf; return(out) }
      g0 <- nfInterp(xA, lo, step, ngrid, G)
    } else g0 <- Gb + (xA - mu_bA) / rb
    out[1] <- g0; out[2] <- Gb
    out[3] <- atilA - exp(-deltaA) * (g0 - Gb)
    return(out)
  })

## Psi_i(u): exposure accumulated since amyloid departure, at age u.
##   u <= alpha_A :  m(mu_bA) (u - alpha_A)                  (<= 0)
##   u >  alpha_A :  e^{-delta_A} [ E(A_i(u)) - E(mu_bA) ]   (>  0)
## Returns NaN if the amyloid path at u is off the table.
nfExposure <- nimbleFunction(
  run = function(u = double(0), g0 = double(0), Gb = double(0), alphaA = double(0),
                 deltaA = double(0), atilA = double(0), mu_bA = double(0),
                 G = double(1), E = double(1), lo = double(0), step = double(0),
                 ngrid = double(0), beta = double(0), cA = double(0), d = double(0)) {
    returnType(double(0))
    if (u <= alphaA) return(nfCouple(mu_bA, beta, cA, d) * (u - alphaA))
    gt <- g0 + exp(deltaA) * (u - atilA)
    if (gt > G[ngrid]) return(NaN)
    y  <- nfInvMono(gt, lo, step, ngrid, G)
    Eb <- nfInterp(mu_bA, lo, step, ngrid, E)
    return(exp(-deltaA) * (nfInterp(y, lo, step, ngrid, E) - Eb))
  })

## Inverse of Psi_i: the age u at which exposure since amyloid departure
## equals psi. Returns (u, A_i(u)) so the caller gets m_i(u) for free.
## u = NaN if the required amyloid level is off the table.
nfInvExposure <- nimbleFunction(
  run = function(psi = double(0), g0 = double(0), Gb = double(0), alphaA = double(0),
                 deltaA = double(0), mu_bA = double(0),
                 G = double(1), E = double(1), lo = double(0), step = double(0),
                 ngrid = double(0), beta = double(0), cA = double(0), d = double(0)) {
    returnType(double(1))
    out <- nimNumeric(2)
    if (psi <= 0) {
      out[1] <- alphaA + psi / nfCouple(mu_bA, beta, cA, d); out[2] <- mu_bA
      return(out)
    }
    Eb <- nfInterp(mu_bA, lo, step, ngrid, E)
    et <- Eb + exp(deltaA) * psi
    if (et > E[ngrid]) { out[1] <- NaN; out[2] <- NaN; return(out) }
    y  <- nfInvMono(et, lo, step, ngrid, E)
    out[1] <- alphaA + exp(-deltaA) * (nfInterp(y, lo, step, ngrid, G) - Gb)
    out[2] <- y
    return(out)
  })

## ---- the implied density of x_tilde_T ------------------------------------
##
##  x_tilde_T is the tau location at the tau reference age atil_T. Its
##  departure age solves  C_i(alpha_T) = -e^{-delta_T} DeltaG_T,  where
##  C_i(a) = Psi_i(a) - Psi_i(atil_T) is exposure since the reference age.
##  The density is the truncated-Weibull marginal of alpha_T, times the
##  copula with the amyloid departure age, times the Jacobian
##      J_T = e^{-delta_T} / ( rbar_T(.) m_i(alpha_T) ).
##  The amyloid marginal and its Jacobian live in dxtilde, so this is the
##  only place the copula is added. Nothing else is a prior on x_tilde_T.
dxtildeT <- nimbleFunction(
  run = function(x = double(0),
                 ## tau
                 deltaT = double(0), atilT = double(0), mu_bT = double(0),
                 xgT = double(0), sigma_alphaT = double(0),
                 aminT = double(0), amaxT = double(0),
                 GT = double(1), rT = double(1), loT = double(0),
                 stepT = double(0), ngridT = double(0), r_minT = double(0),
                 ## amyloid
                 xA = double(0), deltaA = double(0), atilA = double(0),
                 mu_bA = double(0), xgA = double(0), sigma_alphaA = double(0),
                 aminA = double(0), amaxA = double(0),
                 GA = double(1), rA = double(1), EA = double(1), loA = double(0),
                 stepA = double(0), ngridA = double(0), r_minA = double(0),
                 ## coupling and dependence
                 beta = double(0), cA = double(0), d_coup = double(0),
                 rho = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    hiT <- loT + (ngridT - 1) * stepT
    GbT <- nfInterp(mu_bT, loT, stepT, ngridT, GT)
    rbT <- max(nfInterp(mu_bT, loT, stepT, ngridT, rT), r_minT)
    dG <- 0; ljac <- 0
    if (x >= mu_bT) {
      if (x > hiT) { if (log) return(-Inf) else return(0) }
      dG   <- nfInterp(x, loT, stepT, ngridT, GT) - GbT
      ljac <- -log(max(nfInterp(x, loT, stepT, ngridT, rT), r_minT))
    } else {
      dG   <- (x - mu_bT) / rbT
      ljac <- -log(rbT)
    }
    ## the amyloid path this subject is currently on
    st <- nimNumeric(3)
    st[1:3] <- nfAmyState(xA, deltaA, atilA, mu_bA, GA, rA, loA, stepA, ngridA, r_minA)
    alphaA <- st[3]
    if (alphaA == -Inf) { if (log) return(-Inf) else return(0) }
    if (alphaA <= aminA | alphaA > amaxA) { if (log) return(-Inf) else return(0) }
    ## exposure at the reference age, and the departure age that solves
    ## C_i(alpha_T) = -e^{-delta_T} dG
    psi0 <- nfExposure(atilT, st[1], st[2], alphaA, deltaA, atilA, mu_bA,
                       GA, EA, loA, stepA, ngridA, beta, cA, d_coup)
    if (is.nan(psi0)) { if (log) return(-Inf) else return(0) }
    inv <- nimNumeric(2)
    inv[1:2] <- nfInvExposure(psi0 - exp(-deltaT) * dG, st[1], st[2], alphaA, deltaA,
                              mu_bA, GA, EA, loA, stepA, ngridA, beta, cA, d_coup)
    alphaT <- inv[1]
    if (is.nan(alphaT)) { if (log) return(-Inf) else return(0) }
    if (alphaT <= aminT | alphaT > amaxT) { if (log) return(-Inf) else return(0) }
    m_at <- nfCouple(inv[2], beta, cA, d_coup)
    ## marginal, copula, Jacobian
    ld <- nfLogTWeib(alphaT, xgT, sigma_alphaT, aminT, amaxT)
    pA <- nfTWeibCDF(alphaA, xgA, sigma_alphaA, aminA, amaxA)
    pT <- nfTWeibCDF(alphaT, xgT, sigma_alphaT, aminT, amaxT)
    if (pA <= 0 | pA >= 1 | pT <= 0 | pT >= 1) { if (log) return(-Inf) else return(0) }
    ld <- ld + nfLogCopula(pA, pT, rho) - deltaT + ljac - log(m_at)
    if (is.nan(ld)) { if (log) return(-Inf) else return(0) }
    if (log) return(ld) else return(exp(ld))
  })

rxtildeT <- nimbleFunction(
  run = function(n = integer(0),
                 deltaT = double(0), atilT = double(0), mu_bT = double(0),
                 xgT = double(0), sigma_alphaT = double(0),
                 aminT = double(0), amaxT = double(0),
                 GT = double(1), rT = double(1), loT = double(0),
                 stepT = double(0), ngridT = double(0), r_minT = double(0),
                 xA = double(0), deltaA = double(0), atilA = double(0),
                 mu_bA = double(0), xgA = double(0), sigma_alphaA = double(0),
                 aminA = double(0), amaxA = double(0),
                 GA = double(1), rA = double(1), EA = double(1), loA = double(0),
                 stepA = double(0), ngridA = double(0), r_minA = double(0),
                 beta = double(0), cA = double(0), d_coup = double(0),
                 rho = double(0)) {
    returnType(double(0)); return(mu_bT)      # never used: no simulation here
  })

## ---- the three-state subject likelihood ----------------------------------
##
##  x = c(amyloid scans [1:maxA], tau scans [maxA + 1:maxT]), zero-padded.
##
##     l00 = log(1-pi)          + lA0 + lT0
##     l10 = log pi + log(1-q)  + lA1 + lT0
##     l11 = log pi + log q     + lA1 + lT1
##
##  lA1 is the accumulator amyloid path; lT1 is the tau path driven by that
##  SAME amyloid path through the exposure clock. A subject with no tau scans
##  has lT0 = lT1 = 0 and q drops out: the amyloid-only mixture is recovered.
dsubjJoint <- nimbleFunction(
  run = function(x = double(1),
                 nA = double(0), agesA = double(1), nT = double(0), agesT = double(1),
                 offT = double(0),
                 ## amyloid
                 xA = double(0), atilA = double(0), deltaA = double(0), bA = double(0),
                 logit_pi = double(0), mu_bA = double(0),
                 GA = double(1), rA = double(1), EA = double(1), loA = double(0),
                 stepA = double(0), ngridA = double(0), r_minA = double(0),
                 omA = double(1), vloA = double(0), vstepA = double(0), nvA = double(0),
                 ## tau
                 xT = double(0), atilT = double(0), deltaT = double(0), bT = double(0),
                 logit_q = double(0), mu_bT = double(0),
                 GT = double(1), rT = double(1), loT = double(0),
                 stepT = double(0), ngridT = double(0), r_minT = double(0),
                 omT = double(1), vloT = double(0), vstepT = double(0), nvT = double(0),
                 ## coupling
                 beta = double(0), cA = double(0), d_coup = double(0),
                 df = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    hiA <- loA + (ngridA - 1) * stepA
    hiT <- loT + (ngridT - 1) * stepT
    if (xA >= hiA | xT >= hiT) { if (log) return(-Inf) else return(0) }

    ## ---- amyloid: the two branches, as in dsubjMix -------------------
    st <- nimNumeric(3)
    st[1:3] <- nfAmyState(xA, deltaA, atilA, mu_bA, GA, rA, loA, stepA, ngridA, r_minA)
    g0 <- st[1]; Gb <- st[2]; alphaA <- st[3]
    edA <- exp(deltaA)
    lA1 <- 0; lA0 <- 0
    omA0 <- nfInterp(bA, vloA, vstepA, nvA, omA)
    for (j in 1:nA) {
      gt <- g0 + edA * (agesA[j] - atilA)
      if (gt < Gb) gt <- Gb
      if (gt > GA[ngridA]) { if (log) return(-Inf) else return(0) }
      mu <- nfInvMono(gt, loA, stepA, ngridA, GA)
      om <- nfInterp(mu, vloA, vstepA, nvA, omA)
      lA1 <- lA1 + nfLogDT(x[j], mu, om, df)
      lA0 <- lA0 + nfLogDT(x[j], bA, omA0, df)
    }

    ## ---- tau: flat branch, and the coupled accumulator branch ---------
    lT0 <- 0; lT1 <- 0
    if (nT > 0) {
      GbT <- nfInterp(mu_bT, loT, stepT, ngridT, GT)
      rbT <- max(nfInterp(mu_bT, loT, stepT, ngridT, rT), r_minT)
      g0T <- 0
      if (xT >= mu_bT) g0T <- nfInterp(xT, loT, stepT, ngridT, GT)
      else             g0T <- GbT + (xT - mu_bT) / rbT
      edT <- exp(deltaT)
      omT0 <- nfInterp(bT, vloT, vstepT, nvT, omT)
      psi0 <- nfExposure(atilT, g0, Gb, alphaA, deltaA, atilA, mu_bA,
                         GA, EA, loA, stepA, ngridA, beta, cA, d_coup)
      if (is.nan(psi0)) { if (log) return(-Inf) else return(0) }
      for (j in 1:nT) {
        psi <- nfExposure(agesT[j], g0, Gb, alphaA, deltaA, atilA, mu_bA,
                          GA, EA, loA, stepA, ngridA, beta, cA, d_coup)
        if (is.nan(psi)) { if (log) return(-Inf) else return(0) }
        gt <- g0T + edT * (psi - psi0)
        if (gt < GbT) gt <- GbT                 # flat at mu_bT before tau departure
        if (gt > GT[ngridT]) { if (log) return(-Inf) else return(0) }
        mu <- nfInvMono(gt, loT, stepT, ngridT, GT)
        om <- nfInterp(mu, vloT, vstepT, nvT, omT)
        lT1 <- lT1 + nfLogDT(x[offT + j], mu, om, df)
        lT0 <- lT0 + nfLogDT(x[offT + j], bT, omT0, df)
      }
    }

    ## ---- the three states ----------------------------------------------
    lpi  <- -log(1 + exp(-logit_pi));  l1pi <- -log(1 + exp(logit_pi))
    lq   <- -log(1 + exp(-logit_q));   l1q  <- -log(1 + exp(logit_q))
    l00 <- l1pi + lA0 + lT0
    l10 <- lpi + l1q + lA1 + lT0
    l11 <- lpi + lq  + lA1 + lT1
    mx <- max(l00, max(l10, l11))
    ll <- mx + log(exp(l00 - mx) + exp(l10 - mx) + exp(l11 - mx))
    if (is.nan(ll)) { if (log) return(-Inf) else return(0) }
    if (ll == Inf)  { if (log) return(-Inf) else return(0) }
    if (log) return(ll) else return(exp(ll))
  })

rsubjJoint <- nimbleFunction(
  run = function(n = integer(0),
                 nA = double(0), agesA = double(1), nT = double(0), agesT = double(1),
                 offT = double(0),
                 xA = double(0), atilA = double(0), deltaA = double(0), bA = double(0),
                 logit_pi = double(0), mu_bA = double(0),
                 GA = double(1), rA = double(1), EA = double(1), loA = double(0),
                 stepA = double(0), ngridA = double(0), r_minA = double(0),
                 omA = double(1), vloA = double(0), vstepA = double(0), nvA = double(0),
                 xT = double(0), atilT = double(0), deltaT = double(0), bT = double(0),
                 logit_q = double(0), mu_bT = double(0),
                 GT = double(1), rT = double(1), loT = double(0),
                 stepT = double(0), ngridT = double(0), r_minT = double(0),
                 omT = double(1), vloT = double(0), vstepT = double(0), nvT = double(0),
                 beta = double(0), cA = double(0), d_coup = double(0),
                 df = double(0)) {
    returnType(double(1)); return(nimNumeric(length(agesA) + length(agesT)))
  })

registerDistributions(list(
  dxtildeT = list(
    BUGSdist = paste0(
      "dxtildeT(deltaT, atilT, mu_bT, xgT, sigma_alphaT, aminT, amaxT, ",
      "GT, rT, loT, stepT, ngridT, r_minT, ",
      "xA, deltaA, atilA, mu_bA, xgA, sigma_alphaA, aminA, amaxA, ",
      "GA, rA, EA, loA, stepA, ngridA, r_minA, beta, cA, d_coup, rho)"),
    types = c("GT = double(1)", "rT = double(1)", "GA = double(1)",
              "rA = double(1)", "EA = double(1)"), pqAvail = FALSE),
  dsubjJoint = list(
    BUGSdist = paste0(
      "dsubjJoint(nA, agesA, nT, agesT, offT, ",
      "xA, atilA, deltaA, bA, logit_pi, mu_bA, GA, rA, EA, loA, stepA, ngridA, r_minA, ",
      "omA, vloA, vstepA, nvA, ",
      "xT, atilT, deltaT, bT, logit_q, mu_bT, GT, rT, loT, stepT, ngridT, r_minT, ",
      "omT, vloT, vstepT, nvT, beta, cA, d_coup, df)"),
    types = c("value = double(1)", "agesA = double(1)", "agesT = double(1)",
              "GA = double(1)", "rA = double(1)", "EA = double(1)", "omA = double(1)",
              "GT = double(1)", "rT = double(1)", "omT = double(1)"), pqAvail = FALSE)))
