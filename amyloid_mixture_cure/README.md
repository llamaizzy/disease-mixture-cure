# A mixture cure ODE model for amyloid accumulation

This directory is a complete, self-contained implementation of the model: the
data, the code, the fitted posterior, and the figures. Nothing outside it is
needed.

```r
setwd("amyloid_mixture_cure")
source("run_model.R")

fit <- load_fitted_model()      # the posterior behind the reported results
summarise_fit(fit)              # the headline numbers
make_figures(fit)               # the three report figures -> figures/

fit <- fit_amyloid_model("smoke")   # or fit it yourself; see "Running it"
```

Requires R with `nimble`, `splines2`, `ggplot2`, `patchwork`:

```r
install.packages(c("nimble", "splines2", "ggplot2", "patchwork"))
```

To check your installation, run `Rscript tests/test_package.R` from this
directory. It takes about three minutes and should report 22 passed, 0 failed.

---

## 1. The scientific question

Amyloid-β accumulates in the brain for one to two decades before any cognitive
symptom appears. We observe it through PET scans — here 3,393 scans on 1,101
ADNI subjects: a median of three scans each over a median of four years of
follow-up (mean 4.9, maximum 14.1). Baseline ages run 51 to 94, median 72.
40.7% are APOE4 carriers and 48.9% female.

Four years is a short window onto a twenty-year process. Every subject is seen
through a different slice of it, and we do not know which slice. The modelling
problem is to recover one long curve from many short arcs, and to say how
APOE4 carriage and sex change it.

Three different things could change, and the point of this model is that it
separates them:

| | question | parameter |
|---|---|---|
| **susceptibility** | does this person accumulate at all? | φ |
| **latency** | if so, when do they start? | γ |
| **rate** | once started, how fast do they go? | ψ |

A model without the first of these must attribute a lifelong-flat subject to
"has not started yet", which forces the latency distribution to grow a long
right tail and drags the whole timeline with it. That is the reason for the
cure component.

## 2. The model

### 2.1 Three states

Each subject is in one of three states at any age `a`:

- **C — cured / not susceptible.** Never accumulates. Sits at a constant level
  `b_i` for life. Weight `1 − π_i`.
- **L — latent.** Susceptible, but has not departed yet. Sits at the shared
  floor `μ_b`. Weight `π_i (1 − F_α(a))`.
- **P — progressing.** Departed at age `α_i`, now accumulating. Weight
  `π_i F_α(a)`.

Which state a subject is in is never observed. The membership indicator is
**marginalised out**, not sampled: the likelihood is a log-sum-exp over the two
branches (C, versus L-then-P, which is one continuous trajectory). Marginalising
mixes better than a discrete latent variable, and it gives you the posterior
probability of susceptibility for each subject for free.

### 2.2 The accumulation process

For a susceptible subject, the level follows an autonomous ODE:

```
dA/da = exp(δ_i) · r_A(A(a)),        A(α_i) = μ_b
```

Read it as: `r_A(y)` is a single shared curve saying how fast amyloid moves
when it is at level `y`, and `exp(δ_i)` is subject `i`'s personal speed
multiplier on that curve.

**The rate depends on the level, not on age.** This is the central modelling
choice and it does all the work. Because the ODE is autonomous it separates and
integrates exactly. Define the **clock**

```
G_A(y) = ∫ from y* to y  du / max(r_A(u), r_min)
```

(`y*` is the positivity threshold, so `G_A(y*) = 0`). Then the trajectory is a
straight line on the clock scale:

```
G_A(A_i(a)) = G_A(μ_b) + exp(δ_i) · (a − α_i)
```

So the model needs **no numerical ODE solver**. Everything is one tabulated
integral and one monotone inverse lookup. `−G_A(μ_b)` is the **transit budget**:
the number of disease-years an accumulator must spend getting from the floor up
to positivity, at unit speed. In the fitted model it is 10.7 years.

### 2.3 Where the initial condition comes from

Most disease-progression models give each subject a free starting *level*. This
one does not. Everyone launches from the same `μ_b`, and what varies between
subjects is **when** they launch.

The reason is identifiability, and it is worth being precise about. The
trajectory depends on a subject-specific launch level `b_i` and a departure age
`α_i` only through the combination `G_A(b_i) − exp(δ_i) α_i`. Trading one
against the other along that combination leaves every fitted value unchanged —
a perfectly flat direction in the posterior. Fixing the launch level removes it.

The payoff is interpretive as well as numerical: "this subject is further along
than expected" is then explained by a *mechanism* (they departed earlier, or
they are moving faster) rather than by a free parameter.

### 2.4 Departure age

```
α_i = α_min + exp(x_i'γ) · E_i^{σ_α},      E_i ~ Exp(1)
```

a Weibull accelerated failure time model, truncated to `(α_min, α_max] = (25, 125]`.
The AFT parameterisation is chosen so that `exp(γ_q)` reads directly as a
**time ratio**: a value of 0.73 means carriers depart at 73% of the age
(measured from `α_min`) that non-carriers do.

### 2.5 Susceptibility and the non-accumulators

```
logit π_i = w_i'φ
b_i ~ N(μ_b, σ_b²) truncated above at the threshold
```

The truncation on `b_i` is a modelling statement: a subject observed above the
positivity threshold has certainly accumulated, so the non-accumulator branch
cannot be used to explain them.

### 2.6 Measurement

```
y_ij ~ scaled-t_4( A_i(a_ij), ω_A(A_i(a_ij)) )
```

Two things to notice. The scale `ω_A` is a **function of the level**, because
PET noise grows with signal. And the degrees of freedom are fixed at 4 rather
than using a normal: a handful of subjects genuinely accelerate faster than the
shared curve allows, and under a normal likelihood those subjects bend the
shared curve to accommodate themselves.

Both mixture branches use the **same** `ω_A`. If they did not, the mixture
could explain a flat subject by giving them smaller measurement error rather
than by assigning them to a different mechanism, and the cure fraction would
stop meaning what it says.

### 2.7 Smoothing priors

Both curves are splines in level, and both coefficient vectors get an **adaptive
random walk**:

```
θ_j ~ N(θ_{j−1}, ς_j²),     ς_j = sqrt(λ_j / ς_0)
λ_j ~ Exp(1) truncated to (0.01, 100)      local, per increment
ς_0 ~ Gamma(1, 1) truncated to (0.01, 100) global smoothing level
```

The local multiplier `λ_j` lets the curve bend sharply in one place without
having to bend everywhere — a global smoothing parameter alone would over-smooth
the steep part or under-smooth the flat part.

The truncation away from zero is **not cosmetic**. The increment SD is
`sqrt(λ_j/ς_0)`, so `λ_j → 0` drives it to zero and the normal density at its
own mean to `+∞`: the log-density grows like `−½ log λ`. An `Exp(1)` prior is
finite at zero and does not stop this, so the posterior has an integrable-looking
but sampler-catching singularity at `λ_j = 0`. Truncating at 0.01 removes it.

`θ_1` is **fixed, not estimated**. Its basis function has support entirely below
`μ_b`, where the model never evaluates the rate curve, so it carries no
information. Leaving it free adds a parameter the data cannot see.

## 3. The files

```
run_model.R                 the only file you call
R/00_settings.R             every constant, in one place
R/01_integrator.R           trapezoid, interpolation, monotone inverse, the clock
R/02_bases.R                the two spline bases
R/03_densities.R            the custom NIMBLE densities
R/04_model.R                the model, written out in full
R/05_figures.R              the report figures
data/adni_amyloid.rds       1,101 subjects, 3,393 scans
data/fitted_posterior.rds   the posterior behind the reported results
tests/test_package.R        the self-test
```

Functions you will actually call, all from `run_model.R`:

| | |
|---|---|
| `load_amyloid_data()` | the data, as a list |
| `load_fitted_model()` | the shipped posterior, ready for figures |
| `fit_amyloid_model(tier)` | fit it yourself |
| `summarise_fit(fit)` | the headline table |
| `convergence(fit)` | R-hat, on a fit you ran |
| `make_figures(fit)` | all three figures |
| `crossing_ages(fit)` | mean threshold-crossing age by APOE4 |
| `subject_quantities(fit)` | per-subject α, δ, and P(susceptible), within-draw |

`R/04_model.R` is deliberately **flat**. The research version builds this model
through a chain of nine nested files, so that each layer can be tested against
the one below it by pinning a parameter and checking the log-densities agree to
machine precision. That is the right way to develop a model of this kind and the
wrong way to read one. Everything is in one place here.

### Two custom densities, and why they are custom

**`dsubjMix`** — a subject's contribution is not a product of per-visit terms
that BUGS syntax can write one at a time. It is a single function of a whole
trajectory, and for the mixture it is a log-sum-exp over two of them.

**`dxtilde`** — the more interesting one. The model samples `x̃_i`, the
subject's level at their own visit centroid, rather than `α_i` directly. Given
`x̃_i` the departure age is deterministic, so the AFT prior on `α` plus a change
of variables **fixes** the density of `x̃_i`. There is no free prior to choose.
Putting one there anyway would either double-count the AFT or displace it, and
`γ` would stop being identified.

### The atom at the floor

Restricting `x̃ ≥ μ_b` would put a point mass at exactly `μ_b`, of mass
`S_α(ã_i)` — precisely the subjects in state L, susceptible but not yet
departed. A slice sampler cannot reach an atom. The model would therefore rule
state L out silently and push every flat subject into the cured branch,
**inflating the cure fraction by construction**.

The fix is to continue the clock linearly below the floor at the rate there:

```
G_ext(y) = G_A(μ_b) + (y − μ_b)/r̄_A(μ_b)     for y < μ_b
```

which spreads the atom into a density. Below `μ_b`, `x̃` is *only* a coordinate
meaning "departs after the centroid". No biological level below the floor is
implied. The self-test checks that the resulting density is continuous across
`μ_b`, which is the claim that the gap vanishes.

## 4. Running it

```r
fit <- fit_amyloid_model(tier, n_chains = 3)
```

| tier | subjects | iterations per chain | time (3 chains, M-series Mac) | what it is for |
|---|---|---|---|---|
| `"smoke"`  | 250 subsampled | 300 | ~1 min | does everything wire up |
| `"struct"` | all 1,101 | 800 | ~6 min | does it run on the real cohort |
| `"infer"`  | all 1,101 | 10,000 | ~1 h | results |

Sampling runs at **0.126 s per iteration** on the full cohort, essentially
independent of how many iterations you ask for, plus about 1.5 minutes to compile.
The `"infer"` tier is 3 × 10,000 = 30,000 iterations, hence about an hour. Scale
linearly for anything else: the 40,000–50,000 iterations per chain discussed below
would be 120,000–150,000 iterations in total, or about **four to five hours**.

Fit to `results/<tier>_fit.rds`, and checkpointed every chunk, so a run you
interrupt is not lost.

Check convergence on any fit you run with `convergence(fit)`. At the `"struct"`
length (3 × 800) the worst R-hat is 2.02 and 23 of 31 parameters exceed 1.01 —
800 iterations is a structural check, not an inference run, and its output should
not be interpreted.

**Convergence, stated honestly.** At the reported length (3 × 10,000, about an
hour) the worst R-hat is 1.073 and 12 of 32 parameters exceed 1.01. R-hat falls monotonically
with chain length — 1.560, 1.141, 1.097, 1.073 at 5k, 10k, 15k, 20k draws — so
this is a sampler that needs longer, not one that is stuck. Clearing R-hat < 1.01
would take roughly 40,000–50,000 iterations per chain.

### Why the sampler is configured by hand

NIMBLE's defaults will run this model, but badly, because the parameters that
fight each other are not the ones that share a name. Blocking here is **by
role**, measured on the posterior:

- onset location: `μ_b ~ γ_1` at +0.57, `γ_1 ~ σ_α` at −0.36
- rate: `ψ_1 ~ σ_δ` at +0.32
- contrasts (`γ_2:`, `φ_2:`): max 0.23, all mixing well

and `γ_1` is the only one tied to the curve — `corr(γ_1, θ_3) = 0.535`, against
0.163 for `μ_b`. So `γ_1` joins the `θ` block and `μ_b` does not. Four
configurations were compared at 3 × 2,000 from **stationary** starts:

| | block | worst R-hat |
|---|---|---|
| A | {γ₁, μ_b, θ} | 1.573 |
| **B** | **{γ₁, θ}, μ_b alone** | **1.190** ← used |
| C | {μ_b, γ₁}, θ alone | 2.216 |
| D | B + {ψ, σ_δ} | 1.478 |

Comparing from the *initial values* instead ranks how fast chains escape
burn-in, not how well they mix. That mistake was made twice here before it was
caught; if you rerun this comparison, start the chains from a stationary point.

Three further rules, each learned from a failure:

- **`σ_δ` gets a slice sampler, never a conjugate one.** The conjugate
  inverse-gamma draw is valid only when `δ` itself is the stochastic node with
  prior `N(0, σ_δ)`. Here `δ` is non-centred and the stochastic node is
  `z ~ N(0,1)`, whose density contains no `σ_δ` at all. Using the conjugate
  sampler anyway draws from a conditional the model does not have, and
  Metropolis-Hastings cannot repair it — `σ_δ` ranged over [0.000, 4.26] against
  a (0, 2) prior, R-hat 8.1.
- **`μ_b` gets a plain slice sampler, not `AF_slice`.** It sits behind an `−∞`
  wall (the feasibility constraint) that adaptive-factor adaptation cannot see.
- **Each subject's placement and speed move together** in one `RW_block`. They
  trade off directly — start late and go fast, or start early and go slow — and
  moving them one at a time walks the ridge one step at a time.

### Infeasible states return −∞

Every grid excursion and every infeasible departure age returns `-Inf`, never a
large finite value. A large finite value is a *plateau*, and a slice sampler will
walk out onto a plateau and stay there. This rule runs through all of
`R/03_densities.R`.

The same constraint bites at initialisation. Since

```
α_i = ã_i − exp(−δ_i){ G_ext(x̃_i) − G_A(μ_b) }  >  α_min
```

a subject seen at a high level but not yet very old requires a **minimum speed**.
At `δ = 0` a small fraction of subjects violate this, each contributing `-Inf`,
and the model cannot be initialised at all. `build_model()` starts those subjects
at the required speed. That is an *initialisation*, not a prior: it changes where
the chain begins, not what it targets.

## 5. Results

From `summarise_fit(load_fitted_model())`. Posterior median [95% credible
interval], covariates centred so intercepts read at the cohort average.

**Susceptibility** — does this subject accumulate at all? (odds ratio)

| | |
|---|---|
| APOE4 | **5.32** [3.72, 7.82] |
| female | 1.49 [1.04, 2.14] |
| π at cohort average | 0.748 [0.705, 0.788] |

**Latency** — if so, when do they start? (time ratio; below 1 is earlier)

| | |
|---|---|
| APOE4 | **0.730** [0.681, 0.771] |
| female | 0.952 [0.879, 1.018] |

**Rate** — once started, how fast? (rate ratio)

| | |
|---|---|
| APOE4 | 1.092 [0.971, 1.227] — covers 1 |
| female | 0.903 [0.811, 0.996] |

**Shared structure**

| | |
|---|---|
| μ_b (the floor) | 0.585 [0.576, 0.589] SUVR |
| σ_b | 0.059 | 
| σ_α | 0.286 |
| σ_δ | 0.460 |

Implied (`crossing_ages(fit)`): cure fractions of **40.8%** (APOE4−) and
**11.5%** (APOE4+); mean threshold-crossing ages of **77.5** [76.0, 78.8] and
**65.1** [64.3, 66.0] years, a gap of **12.4 years**; a transit budget of
10.4 years.

The crossing age is `α_i + T·exp(−δ_i)` — depart, then spend the transit budget
`T = −G_A(μ_b)` at your own speed — computed within each draw and averaged over
that draw's accumulators.

**The headline.** APOE4 acts overwhelmingly on *whether* and *when*, not on *how
fast*. The rate ratio covers 1. This is the separation the three-part structure
was built to make, and it is not visible to a model that folds all three into a
single "progression" effect.

## 6. The figures

`make_figures(fit)` writes three PDFs to `figures/`.

**`fig1_survival.pdf`** — the rate curve, and the survival curves with and
without the cure component. Panel A shows `r_A(y)`: accumulation is fastest
around 35 Centiloid, just past the positivity threshold, and has fallen to 13%
of its peak by 137 Centiloid. Panels B and C are the same fit with and without
the cure term, and the contrast is the whole argument for the model: B falls to
zero, because in a standard survival model everyone departs eventually; C
plateaus at `1 − π`.

**`fig2_aligned.pdf`** — every accumulator aligned to its own estimated
departure age. Panel A is the calendar axis `a − α_i`; subjects still spread
apart, because each travels at `exp(δ_i)`. Panel B is the clock axis
`exp(δ_i)(a − α_i)`, which is what the model claims removes the remaining
difference. That the arcs collapse onto one curve in B is the model's central
claim, drawn.

**`fig3_fit_check.pdf`** — observed, fitted, and the residual against calendar
age. Panel C is the one that matters: a residual drifting with age would say the
shared curve is the wrong shape.

## 7. Three things that will mislead you

These are not hypothetical. Each one produced a wrong answer here first.

**Do not correlate posterior medians.** The departure age `α_i` and the speed
`δ_i` have posterior medians that correlate at −0.355 on this data. They also
correlate at −0.335 on data simulated with `α ⊥ δ` **by construction**. Computed
*within each posterior draw*, the real-data correlation is +0.002 [−0.070, 0.063]
— no association at all. The estimator coordinates the two, and taking medians
first bakes that coordination in as if it were a finding. Figure 2 uses medians
to *align* the curves, which is a display choice; it must not be used to read
off a relationship between onset and rate.

**Check predictive summaries against a refit, not just a replicate.** Posterior
replicates cannot see artefacts that are caused by *estimating* the random
effects, because the replicate conditions on the estimates. Simulate from the
fit and then **refit**. A worrying misfit in the top of the amyloid range turned
out this way to be shrinkage plus conditioning on an estimated disease age, and
not a problem with the curve at all.

**A Centiloid transform is affine.** `CL = 169.27 · SUVR − 102.35`. A *level*
transforms with both terms; a *difference* of levels scales only. Applying the
full transform to a residual once emptied a whole diagnostic panel here.

## 8. Things to try

- Drop the cure component (set `φ_1` very large so `π ≈ 1`) and refit. Watch
  what the latency distribution has to do to absorb the flat subjects. This is
  the cleanest way to see what the mixture is buying.
- Drop `b_i` (fix `σ_b = 0`). In the research version this pushed `π` to 0.95
  and cost the APOE4 susceptibility effect its credibility — the subject-level
  intercept is what lets the model tell "flat at a low level" apart from "flat
  at a high level".
- Add an APOE4 × sex interaction to each of `φ`, `γ`, `ψ`. All three came back
  null here; confirming a null is a legitimate exercise.
- Put covariates on `σ_δ`, so that speed *variability* — not just mean speed —
  can depend on genotype.
- Replace the scaled-`t` with a normal and look at what happens to the shared
  curve. The fast accumulators are the test case.
