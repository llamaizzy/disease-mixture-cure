# Simulation study

Does the mixture cure ODE model recover what it claims to recover, under the
follow-up pattern it was actually fitted to, and can it tell *earlier* from
*faster*? This directory answers that with simulated cohorts, at a compute cost
that fits on a laptop.

```r
source("run_simulation.R")
run_simulation("pilot")     # ~2 min:  does everything wire up
run_simulation("bench")     # ~2 min:  measure this machine's throughput
run_simulation("small")     # ~4 h:    400 subjects, all four scenarios
run_simulation("full")      # ~4 h:    the real 1,101-subject design, primary scenario
summarise_sim("small")      # tables + figures from whatever has finished
```

Every replicate is written to `results/sim/<tier>/` the moment it is scored.
Kill the run whenever you like; calling the tier again resumes it.

---

## 1. Protocol (ADEMP)

**Aims.**
1. Parameter recovery at the real design: bias, RMSE and 95% interval coverage
   for every quantity the report states.
2. Separation of the three effects. The headline finding is that APOE4 acts on
   *whether* and *when*, not on *how fast*. A model that cannot distinguish an
   earlier start from a faster rate would produce that headline whether or not
   it were true.
3. The cure fraction. The extended clock (README section 3) exists to stop the
   sampler inflating the cure fraction by construction. Does the estimated cure
   fraction stay honest when the true one is small?
4. Subject-level recovery: membership, departure age, speed.

**Data-generating mechanism.** The *design* is the real ADNI cohort, fixed:
each subject's visit ages, number of scans and covariates are exactly as
observed. Only the amyloid measurements are simulated, from the model in
`R/04_model.R` with the truth below. Two reasons for fixing the design. The
question is whether the model works under *this* study's four-year windows
onto a twenty-year process, not under an idealised design. And a fixed design
is what lets one compiled model serve every replicate, which is where the
compute savings come from (section 3).

Four truths, each the shipped posterior median with one change:

| scenario | what changes | what it tests |
|---|---|---|
| `fitted` | nothing: OR 5.3, TR 0.73, RR 1.09 | recovery at the real design |
| `rate_only` | APOE4 time ratio set to 1, rate ratio to 1.5 | can it tell *faster* from *earlier*? |
| `low_cure` | π = 0.95 for everyone | is a 5% cure fraction reported as 5%? |
| `null` | every covariate effect 0 | false positives; interval coverage at the null |

The variance basis is built from the full real cohort for every tier, so "the
truth ν" means the same function whatever the subsample or replicate.

**Estimands.** On the scales the report uses: the three APOE4 and three sex
effects (odds, time and rate ratios), π at the cohort average, μ_b, σ_b, σ_α,
σ_δ, the transit budget, log r(y) at ten levels from 0.6 to 1.5 SUVR, the
cure fraction by APOE4, and the mean threshold-crossing age by APOE4 among
accumulators. Per subject: P(susceptible) against true membership; α_i and
δ_i for true accumulators, split by whether they had departed by their visit
centroid (state P) or not (state L).

**Methods.** The model fitted exactly as `fit_amyloid_model()` fits it: the
same samplers, the same data-driven initial values (never the truth), the same
jitter for chains 2+, the same burn-in rule. Nothing is tuned to the simulation.

**Performance measures.** Bias, relative bias, empirical SE, RMSE, coverage of
the 95% credible interval and mean interval width, each with its Monte Carlo
standard error. Per-replicate worst R-hat is recorded so that a coverage
failure can be told apart from a convergence failure.

## 2. What the output looks like

`summarise_sim(tier)` writes to `figures/sim/<tier>/`:

- `performance.csv` — one row per scenario × estimand
- `subject_level.csv` — membership AUC, α and δ recovery, the α–δ correlation
  computed two ways, R-hat, minutes per replicate
- `fig_sim1_intervals.pdf` — every replicate's median and interval against truth
- `fig_sim2_coverage_bias.pdf` — coverage and standardised bias with MC error
- `fig_sim3_rate_curve.pdf` — the rate curve: truth, spread of estimates, pointwise coverage
- `fig_sim4_subjects.pdf` — P(susceptible) by true state, AUC, the α–δ pitfall,
  subject-level coverage

and prints a compact table to the console.

## 3. The compute budget, and how it was kept

A fit of the full cohort costs about 70 s to compile and 0.117 s per MCMC
iteration on an 11-core M3 Pro, one process. The reported analysis is
3 × 10,000 iterations, about an hour. A naive simulation, compiling and running
that per replicate, is 100 hours for 100 replicates. Three things bring it down.

**Compile once, fit many.** The design is the only thing baked into a compiled
NIMBLE model. Each worker compiles once, then for every replicate writes the
simulated scans into the compiled model's data node, resets the initial values
and reruns the same compiled MCMC (`simulation/sim_fit.R`). Compile cost drops
from once per replicate to once per worker.

**Parallel workers with checkpoints.** Replicates are independent, so they run
on a PSOCK cluster, default 8 workers (`detectCores() - 3`, capped at 8). The
M3 Pro's six efficiency cores are slower than its five performance cores, so
eight busy workers deliver about six single-process equivalents; tier `bench`
measures the real factor on your machine and the planner uses it.

**Tiers sized to a laptop.** Iteration cost is linear in subjects (0.040 s at
400), so a 400-subject tier is three times cheaper per iteration than the full
design. It is a different design, and its results say what the model does at
n = 400; the `full` tier is what speaks to the paper's n = 1,101.

| tier | subjects | chains × iterations | replicates | one process | 8 workers (measured) |
|---|---|---|---|---|---|
| `pilot` | 250 | 2 × 400 | 5 | 3 min | 2 min |
| `bench` | 400 | 2 × 200 | 1 per worker | — | 1.5 min |
| `small` | 400 | 3 × 3,000 | 100 + 40 + 40 + 40 | 31 h | ~4.1 h |
| `full` | 1,101 | 3 × 6,000 | 40 (`fitted`) | 31 h | ~3.9 h |

Measured on the M3 Pro: with eight workers busy, each runs at 1.32× the
single-process time per iteration, so eight workers deliver about six
single-process equivalents. Eight compiled 400-subject models used 3.6 GB in
total; the full design will need roughly three times that, well inside 18 GB.
The `small` and `full` estimates above come from that bench run.

Two practical points for a laptop. Run a long tier from a terminal with
`caffeinate -i Rscript -e 'source("run_simulation.R"); run_simulation("small")'`
so the machine does not sleep. And if you need the machine back, kill the run
and rerun the same call later: it resumes from the saved replicates, losing at
most one in-flight replicate per worker.

`run_simulation(tier, dry_run = TRUE)` prints the plan and the estimate without
running. Every argument can be overridden:

```r
run_simulation("full", scenarios = "rate_only", n_rep = 20)   # add a scenario to the full tier
run_simulation("small", n_iter = 6000, n_workers = 6)         # longer chains, fewer workers
```

**Adding to the full tier.** The separation test (`rate_only`) runs at 400
subjects in the `small` tier. To run it at the real design, or to match the
reported 3 × 10,000 exactly:

```r
run_simulation("full", scenarios = "rate_only", n_rep = 20)      # ~2 h more
run_simulation("full", n_iter = 10000, n_rep = 24)               # ~4 h; same cost as 40 x 6,000
```

**On chain length.** The `small` and `full` tiers use fewer iterations than the
reported analysis (3 × 3,000 and 3 × 6,000 against 3 × 10,000). The README
reports worst R-hat 1.56 at 5,000 draws and 1.07 at 10,000, so some replicates
will not be fully converged, and the worst offenders are the spline
coefficients rather than the covariate contrasts. Every replicate records its
R-hat; `subject_level.csv` reports the mean worst R-hat per scenario, and
`fig_sim1` lets you see whether the interval misses cluster in poorly-mixed
replicates. If they do, rerun with `n_iter` doubled and half the replicates:
the cost is the same, and the comparison tells you whether the miss was
coverage or convergence.

## 4. Files

```
run_simulation.R             the only file you call
simulation/sim_truth.R       the four scenarios and the estimands
simulation/sim_generate.R    the data-generating process
simulation/sim_fit.R         compile once, fit a replicate on the compiled model
simulation/sim_summarise.R   score a replicate; aggregate a tier; figures
simulation/sim_run.R         tiers, worker pool, checkpoints, the planner
results/sim/<tier>/          one .rds per scored replicate; log.txt
figures/sim/<tier>/          tables and figures
```
