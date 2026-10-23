## =====================================================================
##  Spline bases on the amyloid LEVEL grid.
##
##  Both curves are functions of level, not of age or disease stage. That
##  is what makes the ODE autonomous, and autonomy is what gives the
##  trajectory a closed form (see R/04_model.R).
##
##  These match the fitted model exactly -- do not substitute another
##  basis unless you intend to refit, since the saved posterior is written
##  in these coefficients.
## =====================================================================
suppressPackageStartupMessages(library(splines2))

## Rate curve r_A(y). Evenly spaced knots. The basis rows sum to one, which
## matters: it means adding a constant to theta scales the whole curve by a
## constant factor, so the level of theta is a single interpretable direction.
make_rate_basis <- function(S = SETTINGS) {
  lo <- S$grid_lo; hi <- S$grid_hi
  y  <- seq(lo, hi, by = S$grid_step)
  nk <- S$K_rate - 3L - 1L                       # cubic
  kn <- seq(lo, hi, length.out = nk + 2L)[-c(1L, nk + 2L)]
  B  <- bSpline(y, knots = kn, degree = 3L, intercept = TRUE,
                Boundary.knots = c(lo, hi))
  stopifnot("rate basis rows must sum to 1" = max(abs(rowSums(B) - 1)) < 1e-8)
  list(y = y, n = length(y), step = S$grid_step, lo = lo, hi = hi,
       B = B, K = ncol(B), knots = kn)
}

## Measurement scale omega_A(y). Knots at INTERIOR QUANTILES of the observed
## levels, so each basis function is informed by a comparable number of scans
## and none is spent where there is almost no data on either side.
make_var_basis <- function(y_obs, S = SETTINGS) {
  lo <- S$grid_lo; hi <- S$grid_hi
  y  <- seq(lo, hi, length.out = 201L)
  kn <- as.numeric(quantile(y_obs, seq(S$var_knot_q[1], S$var_knot_q[2],
                                       length.out = S$K_var - 1L)))
  kn <- unique(kn[kn > lo & kn < hi])
  B  <- bSpline(y, knots = kn, degree = 3L, intercept = TRUE,
                Boundary.knots = c(lo, hi))
  list(y = y, n = length(y), step = y[2] - y[1], lo = lo, hi = hi,
       B = B, K = ncol(B), knots = kn)
}
