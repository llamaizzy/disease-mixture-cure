## =====================================================================
##  Custom NIMBLE samplers for the joint model.
##
##  THE FLOOR-SHIFT (interweaving) SAMPLER
##
##  Diagnosis (posterior 75k NEO modified/rhat_diagnosis.md): the floors
##  mu_b and mu_b_T have full conditionals, given the 1,101 subjects'
##  location coordinates, that are far narrower than their marginals, so
##  they can only move as fast as the subject layer re-sorts. That is the
##  centred-hierarchical geometry, and its standard cure is to ALTERNATE
##  parameterisations (Yu & Meng 2011, "ASIS"): in addition to the usual
##  update of mu_b with the subjects fixed, propose a TRANSLATION of the
##  floor together with every subject coordinate measured from it,
##
##      mu_b -> mu_b + d,   x_tilde[i] -> x_tilde[i] + d,   b[i] -> b[i] + d,
##
##  so each subject's position relative to the floor is unchanged. The
##  map is a translation (Jacobian 1) and the proposal is symmetric, so
##  the Metropolis ratio is the posterior ratio alone. Coordinates below
##  the floor keep their clock value exactly ((x - mu_b)/rbar(mu_b) is
##  invariant); above it the clock difference G(x+d) - G(mu_b+d) changes
##  slightly, which is what makes the move informative rather than a
##  pure relabelling. A translation that pushes a b[i] above its threshold
##  or an x_tilde off the table gets -Inf and is rejected.
##
##  One proposal costs a full-model evaluation (every subject depends on
##  the floor), the same as any global parameter here. The scale adapts
##  toward a 0.44 acceptance rate.
## =====================================================================
suppressPackageStartupMessages(library(nimble))

sampler_floorShift <- nimbleFunction(
  name = "sampler_floorShift",
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    floorNode  <- target                                 # e.g. "mu_b"
    shiftVars  <- control$shiftVars                      # e.g. c("x_tilde", "b")
    stopifnot(length(shiftVars) == 2L)
    nodes1 <- model$expandNodeNames(shiftVars[1])        # e.g. x_tilde[1:n]
    nodes2 <- model$expandNodeNames(shiftVars[2])        # e.g. b[1:n]
    calcNodes  <- model$getDependencies(c(floorNode, nodes1, nodes2))
    scale0     <- if (is.null(control$scale)) 0.005 else control$scale
    scale      <- scale0
    adaptive   <- if (is.null(control$adaptive)) TRUE else control$adaptive
    adaptInterval <- 200
    timesRan <- 0; timesAccepted <- 0; timesAdapted <- 0
    optimalAR <- 0.44; gamma1 <- 0
  },
  run = function() {
    d <- rnorm(1, 0, scale)
    lp0 <- model$getLogProb(calcNodes)
    values(model, floorNode) <<- values(model, floorNode) + d
    values(model, nodes1) <<- values(model, nodes1) + d
    values(model, nodes2) <<- values(model, nodes2) + d
    lp1 <- model$calculate(calcNodes)
    accept <- FALSE
    if (!is.nan(lp1) & lp1 > -Inf) accept <- decide(lp1 - lp0)
    if (accept) {
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
    if (adaptive) adaptiveProcedure(accept)
  },
  methods = list(
    adaptiveProcedure = function(accepted = logical()) {
      timesRan <<- timesRan + 1
      if (accepted) timesAccepted <<- timesAccepted + 1
      if (timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        gamma1 <<- 1 / ((timesAdapted + 3)^0.8)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        timesRan <<- 0; timesAccepted <<- 0
      }
    },
    reset = function() {
      scale <<- scale0
      timesRan <<- 0; timesAccepted <<- 0; timesAdapted <<- 0; gamma1 <<- 0
    }
  )
)
