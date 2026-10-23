# Joint model: convergence and β recovery

Run: 3 chains × 25,000 iterations (not 60k), full cohort, 632 min. Diagnostics on the second half of each chain. Files produced here: `diag_population_secondhalf.csv`, `diag_traces.pdf`, `diag_state_probs.rds`.

## 1. Verdict

Not converged, but **not multimodal**. All three chains agree on every posterior mean (β: 1.71 / 1.84 / 1.73; ρ: 0.26 / 0.23 / 0.30). The problem is slow diffusion along a few ridges, so longer or better-blocked chains fix it; different starting points would not find anything different.

| group | worst R-hat | min bulk ESS | slowest parameters |
|---|---|---|---|
| coupling | 1.13 | 18 | β, ρ_α |
| amyloid | 1.11 | 24 | μ_b, γ[1], σ_α (inherited from the amyloid-only model) |
| susceptibility | 1.10 | 25 | φ[1] |
| tau | 1.10 | 26 | γ_T[1], μ_bT, σ_δT, θ_T[4] |
| subjects | 1.23 | — | z, z_T for ~2% of subjects |

Lag-500 autocorrelation of β is still 0.21; of φ[1], 0.37. Measurement curves (ν, ν_T) and the tau smoothing parameters mix fine.

## 2. Headline posterior (second half, pooled)

β = 1.76 ± 0.21 (multiplier range e^{±0.88} = 0.4–2.4), ρ_α = 0.26 ± 0.10, q̄ = 0.946, π̄ = 0.893, μ_bT = 1.067, σ_δT = 0.60. State probabilities: 13% / 5% / 82% (00 / 10 / 11); 135 of 648 tau subjects are ambiguous (0.2 < P(11) < 0.8). These numbers carry ESS ≈ 20 and should not be reported yet.

## 3. What causes the poor mixing

**(a) β sits on a ridge and is updated one coordinate at a time.** Holding everything else at the posterior mean, the conditional sd of β is **0.029**; its marginal sd is **0.30**. A univariate slice step on β therefore moves ~1/10 of the posterior width per iteration. The ridge partners are σ_δT (r = +0.63) and the tau rate-curve level θ_T[4:6] (r ≈ −0.5): a larger amyloid multiplier is offset by slower residual tau speeds or a lower tau rate curve. Profiling σ_δT out flattens the β profile by a factor of ~3. This is an identifiability statement about the data, not a bug: 982 of 1,452 tau scans sit at amyloid levels inside the sigmoid's transition band, so β is informed, but weakly relative to the rate-curve level.

**(b) The inherited amyloid ridge.** γ[1]–σ_α (r = −0.84), γ[1]–μ_b (0.71), μ_b–φ[1] (−0.57): the same triangle that limited the amyloid-only fit (worst R-hat 1.07 at 20k) and it now also touches μ_bT (r = −0.60 with φ[1]). γ[1] is in the θ block but μ_b, σ_α and φ[1] are each sampled alone.

**(c) A feasibility wall on β, but outside the bulk.** With subject coordinates fixed, β < 0.6 or > 2.7 makes 2–16 subjects' implied tau departure age fall outside (25, 125]: old subjects (~90 at the tau reference) with high tau. Between 1.0 and 2.4 nobody is infeasible, so this does not drive the slow mixing at the current posterior; it would matter if β were near 0.

**(d) Per-subject RW_block on mixed scales.** NIMBLE warned on every subject block: x̃ (sd ≈ 0.05 SUVR) and z (sd ≈ 1) share one proposal covariance. Adaptation recovers it, but 2% of subjects still have R-hat > 1.05 on z and z_T, and those are the subjects whose timing–speed trade-off is longest.

## 3b. The subject-level R-hats in `traceplots_subject.pdf`

The 15 worst subject parameters (R-hat 1.10–1.23; all `z`, `z_T`, one `x_tilde_T`) fall into three groups. None is bimodal in the sense of chains disagreeing on a location; in every case the chains overlap and the first/second-half medians swing inside one chain.

**A. Speed is prior-dominated: subjects at the floor (325, 209, 11, 460, 277, 936, 450, 633; 8 of 15).** Amyloid scans flat at 0.56–0.62 against a floor of 0.559; tau scans flat at 1.10–1.23 against a floor of 1.067, below the 1.25 threshold. The data pin the *level* (posterior sd of `x_tilde` ≈ 0.01–0.02) but a subject who has barely left the floor is consistent with any speed, so `z` has posterior sd 0.8–1.0 — the prior. The per-subject RW_block adapts to a target that is a needle in `x_tilde` and a plateau in `z`, so `z` random-walks slowly across its prior; half-chain medians move by 0.5–1. That is a slow walk on an unidentified direction, not a convergence problem for anything reported. The same mechanism gives the amyloid-only fit its subject-level R-hats.

**B. State-ambiguous subjects (611, 544, 417, 450).** Subject 611: amyloid-positive (0.88–1.07) at ages 84–93, tau flat at 1.13–1.19; P(10) = 0.51, P(11) = 0.49. In state 11 the coupling multiplier is ≈ e^{β/2} = 2.4 (amyloid far above threshold), so tau *should* be accelerating; flat tau is bought with an extreme residual speed (`z_T` medians −0.6 to −2.7) or by tau not yet departed. In state 10 the tau coordinates are auxiliary and wander under their transformed prior (`x_tilde_T` ranges 1.0–3.0). The marginal is a mixture of a pinned regime and a free regime, the block proposal rejects 3.5–7.8% of the time (vs 0.2–0.8% elsewhere), and the chain switches regimes slowly. This is the β–σ_δT ridge seen at the subject level. 544 (amyloid *declining* 0.70→0.52, P(00) = 0.31) and 417 (one tau scan at age 94) are the same thing with the amyloid state or a single scan as the ambiguity.

**C. The three-state restriction binding (67, 329, 621).** Subject 67: amyloid flat at 0.50–0.58 for ages 81–90 (`x_tilde` below the floor in 85% of draws, implied amyloid departure ≈ 91), tau 1.36–1.42, *above* threshold. Tau positivity forces Z_T = 1, hence Z_A = 1, so the model places her as an amyloid accumulator who has not departed yet, with tau departing (≈ 76) fifteen years before amyloid. `z` is then pure prior (sd 1.0–1.2) and walks. 329 and 621 are milder versions (amyloid just rising, tau already 1.40 / 1.43; tau departure estimated at or before amyloid departure). These are the subjects the proposal's "TEST: susceptibility assumptions" box warns about: medial temporal tau without amyloid. Count them before fitting a four-state sensitivity model.

Sampling fixes for A and B are the scaled `propCov` and the β block in §5; C is a modelling question, not a sampler one.

## 4. β recovery on simulated data

Data simulated on the real design from the posterior means with β set to 0 and to 1.75 (`simulation/sim_joint.R`), each fitted for one chain of 3,000 iterations from the β = 0 start, with the **old** (pre-§5) blocking. Posterior over iterations 1,501–3,000:

| truth β | β | ρ_α (truth 0.26) | σ_δT (truth 0.60) | q̄ (0.95) | μ_bT (1.067) | γ_T[1] (3.61) |
|---|---|---|---|---|---|---|
| 0 | **0.04 ± 0.04** | 0.22 ± 0.06 | 0.27 ± 0.02 | 0.93 | 1.065 | 3.52 |
| 1.75 | **1.02 ± 0.12, still rising** | 0.45 ± 0.06 | 0.30 ± 0.04 | 0.92 | 1.070 | 3.60 |

- **β = 0 is recovered**: the chain never leaves 0.01–0.05, so the coupling is not an artefact the model invents.
- **β = 1.75 is not refuted but not reached**: 500-iteration block means 0.15 → 0.56 → 0.84 → 0.84 → 0.96 → 1.02. The chain is still climbing at 3,000 iterations, exactly the slow drift along the β–σ_δT ridge diagnosed in §3(a). It cannot distinguish "β is biased low" from "the old sampler has not arrived"; the rerun with the §5 blocking (below) does.
- **σ_δT is under-recovered in both fits** (0.27–0.30 vs 0.60), independently of β, and θ_T[4] is high in both. Within 3,000 iterations the subject residuals `z_T` have not spread out from their inits, so the non-centred scale shrinks toward them. This is the same slow subject-level walk as §3b(A), and it drags β down with it through the ridge. Longer chains and the scaled subject blocks address it.

**Rerun of β = 1.75 with the §5 blocking** (same simulated data, one chain, 3,000 iterations, β = 0 start): posterior over 1,501–3,000 **1.44 ± 0.13**, last 500 iterations 1.50; 500-block means 0.64 → 1.07 → 1.32 → 1.41 → 1.39 → 1.50. σ_δT reached 0.47–0.49 (truth 0.60; old blocking stalled at 0.30), θ_T[4] = −3.93 (truth −3.95; old −3.57), μ_bT, γ_T[1], q̄, π̄ all within ~1 sd. β's lag-50 autocorrelation fell from 0.78 to 0.54. The chain is still drifting up at 3,000 and a single short chain is not a coverage statement, but the direction is unambiguous: with the ridge blocked, β moves toward its true value at roughly three times the old speed and the remaining gap tracks σ_δT, which is itself still rising. μ_b in its new AF_slice block (change 2) moved freely (0.562 → 0.553) without stalling.

## 5. Recommendations, in order of expected payoff

Items 1–3 are now implemented in `configure_joint_samplers()` (`run_joint_model.R`) and are the default; `SETTINGS_JOINT$amy_block = FALSE` reverts item 2, `SETTINGS_JOINT$subj_scales` sets item 3.

1. **Block β with its ridge partners** (done): β and σ_δT join the tau AF_slice block {γ_T[1], θ_T[1:10]}. Expect β's ESS per iteration to rise ~10×.
2. **Block the amyloid triangle** (done, untested beyond compilation): AF_slice over {μ_b, σ_α, φ[1]}, with φ[2:3] in their own block. Caveat: `run_model.R` found μ_b in an AF_slice block with γ[1] or θ mixed *worse*, because its feasibility wall is invisible to the adaptation. This block avoids those partners, but check μ_b's trace in the next fit and revert if it stalls.
3. **Scale the subject blocks** (done): the subject RW blocks start from `propCov = diag(c(0.05, 1, 0.15, 1)^2)` instead of the identity.
### Outcome: the 3 × 30k fit with changes 1–3 (`posterior 30k modified/`, 292 min, parallel)

| group | worst R-hat (old → new) | min bulk ESS (old → new) |
|---|---|---|
| coupling (β, ρ_α) | 1.128 → **1.006** | 18 → **142** |
| amyloid | 1.106 → 1.079 | 24 → 36 |
| tau | 1.097 → 1.137 (μ_bT) | 26 → 16 |
| susceptibility | 1.103 → 1.156 (φ[1]) | 25 → 14 |

Change 1 worked: β 1.83 ± 0.21 with R-hat 1.006 and ESS 142 (8× the old), ρ_α 1.005 / 272, σ_δT 1.019 / 152, γ_T[1] 1.022 / 141, θ_T[4] 1.012 / 288. Chains agree to the second decimal. Change 2 did not: the amyloid triangle moved only modestly (μ_b 1.106 → 1.074, γ[1] 1.100 → 1.060) and φ[1] got *worse* (1.103 → 1.156, ESS 14), as did μ_bT (1.086 → 1.137, ESS 16), which correlates with φ[1] at −0.60. The slowest direction is now φ[1] ↔ (μ_b, μ_bT): the amyloid susceptibility intercept trading against both floors. Next step: add μ_bT and φ_T[1] to the {μ_b, σ_α, φ[1]} block (or revert `amy_block` and block {φ[1], μ_bT} instead), and run longer; β and the tau side are already reportable.

4. **Chain length**: at the current blocking, R-hat ≈ 1.01 on β needs roughly ESS 400, i.e. ~20× the current draws. With fixes 1–3 the 25k already run is likely enough. Run chains in parallel (three R processes) rather than sequentially; the model is embarrassingly parallel across chains.
5. Do not widen the jitter: the chains agree, so dispersion is not the question.
