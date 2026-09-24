## =====================================================================
##  03_integrator.R -- THE integrator. One cumulative-trapezoid rule on
##  the value grid, one monotone inverse lookup, one substep quadrature
##  on the age axis. The likelihood, the coupling integral Psi, the
##  hazard integral and every diagnostic call these and nothing else.
##
##  A different method for diagnostics would introduce integrator
##  disagreement into the residuals, indistinguishable from model misfit.
## =====================================================================

## Cumulative trapezoid of f over a uniform grid, starting at 0.
sa_cumtrapz <- function(f, step) {
  n <- length(f)
  c(0, cumsum((f[-1] + f[-n]) / 2) * step)
}

## Linear interpolation on a uniform grid. Out of range returns NA_real_,
## which every caller must convert to -Inf rather than clamping: clamping
## would make the likelihood flat outside the grid.
sa_interp <- function(xg, yg, xout) {
  step <- xg[2] - xg[1]; n <- length(xg)
  pos  <- (xout - xg[1]) / step + 1
  out  <- rep(NA_real_, length(xout))
  ok   <- is.finite(pos) & pos >= 1 & pos <= n
  if (any(ok)) {
    p  <- pos[ok]; i0 <- pmin(floor(p), n - 1); w <- p - i0
    out[ok] <- yg[i0] * (1 - w) + yg[i0 + 1] * w
  }
  out
}

## Inverse of a STRICTLY INCREASING tabulated function. Same convention:
## NA outside the representable range, never a clamp.
sa_inv_monotone <- function(xg, yg, yout) {
  n <- length(yg); out <- rep(NA_real_, length(yout))
  ok <- is.finite(yout) & yout >= yg[1] & yout <= yg[n]
  if (any(ok)) {
    j  <- findInterval(yout[ok], yg, rightmost.closed = TRUE)
    j  <- pmin(pmax(j, 1L), n - 1L)
    dy <- yg[j + 1] - yg[j]
    w  <- ifelse(dy > 0, (yout[ok] - yg[j]) / dy, 0)
    out[ok] <- xg[j] * (1 - w) + xg[j + 1] * w
  }
  out
}

## Uniform substep grid over [a0, a1] with width at most SA$substep_yr.
## Always at least two points so the trapezoid rule is defined.
sa_substeps <- function(a0, a1) {
  m <- max(1L, as.integer(ceiling(abs(a1 - a0) / SA$substep_yr)))
  seq(a0, a1, length.out = m + 1L)
}

## ---- the clock, in plain R -------------------------------------------
## Used by every post-processing routine. It must be the SAME rule as the
## one inside the likelihood (nfClock in R/03_densities.R), or the
## difference between the two shows up in the residuals and is
## indistinguishable from model misfit.

## r_A(y) on the level grid, from spline coefficients.
rate_grid <- function(basis, theta) as.numeric(exp(basis$B %*% theta))

## G(y) = int_{thresh}^{y} du / max(r(u), r_min), tabulated on the grid.
clock_grid <- function(basis, r, thresh, r_min = SETTINGS$r_min) {
  G <- sa_cumtrapz(1 / pmax(r, r_min), basis$step)
  G - sa_interp(basis$y, G, thresh)
}
