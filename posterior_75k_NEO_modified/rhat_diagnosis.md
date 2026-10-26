# Why nine parameters sit at R-hat 1.01–1.02 after 3 × 150k

Fit: `posterior 75k NEO modified/`, 3 chains × 150,000, 30% burn-in (105,000 draws per chain). 69 of 78 population parameters are below 1.01. The nine that are not: γ_T[2] 1.019, φ[1]/π̄ 1.019, ψ_T[1] 1.019, μ_bT 1.019, γ_T[1] 1.017, σ_δT 1.017, β 1.014, μ_b 1.012. Everything below asks the same question of each: is R-hat elevated because the chains *disagree*, because they are *not stationary*, because the marginals are *heavy-tailed*, or simply because the chains are *slow* and 105k draws is not yet enough?

## 1. Decomposition: between-chain variance against what the autocorrelation predicts

R-hat exceeds 1 when the variance of the three chain means, B/n, exceeds what a stationary chain with the observed effective sample size would produce, W/ESS_chain. The ratio of the two separates "slow" from "disagreeing".

| parameter | R-hat | ESS | chain means (sd units) | B/n ÷ W/ESS_c | within-chain split-R-hat (ch 1/2/3) |
|---|---|---|---|---|---|
| φ[1] | 1.019 | 154 | −0.16 +0.09 +0.08 | **1.3** | 1.020 1.030 1.004 |
| μ_bT | 1.019 | 351 | +0.14 −0.12 −0.02 | **1.2** | 1.004 1.045 1.009 |
| γ_T[1] | 1.017 | 357 | +0.12 −0.11 0.00 | **1.3** | 1.002 1.037 1.021 |
| σ_δT | 1.017 | 339 | −0.09 +0.12 −0.03 | **1.3** | 1.009 1.034 1.017 |
| β | 1.014 | 535 | −0.04 +0.02 +0.02 | 0.2 | 1.013 **1.044** 1.009 |
| μ_b | 1.012 | 382 | +0.15 −0.08 −0.07 | 2.3 | 1.009 1.008 1.007 |
| γ_T[2] | 1.019 | 425 | −0.12 **+0.19** −0.07 | **5.5** | 1.003 1.025 1.014 |
| ψ_T[1] | 1.019 | 407 | −0.17 **+0.19** −0.01 | **5.4** | 1.002 1.025 1.005 |
| *gamma[2]* (reference) | 1.004 | 1857 | 0.00 −0.05 +0.05 | 1.2 | 1.002 1.015 1.001 |

Reading: for φ[1], μ_bT, γ_T[1], σ_δT the chains agree *as well as their autocorrelation allows* (ratio 1.2–1.3, same as the well-mixed references). Their R-hat is elevated only because ESS is 150–360. For β the chains agree better than expected; its 1.014 is within-chain non-stationarity in chain 2 (split-R-hat 1.044). γ_T[2] and ψ_T[1] are the two genuine exceptions: chain 2 sits 0.19 sd above the other two, five times the spread autocorrelation predicts.

## 2. Drift

Per-chain means in consecutive 21k blocks (sd units): chain 2 descends on γ_T[2] (+0.35, +0.28, +0.38, +0.06, −0.14), ψ_T[1] (+0.34 → +0.10), σ_δT (+0.47 → −0.07) and β (+0.20, +0.27, +0.08, −0.22, −0.21), while μ_bT and γ_T[1] rise (−0.25 → +0.27; −0.26 → +0.28). Chains 1 and 3 oscillate without trend. So chain 2 spent its first ~60k retained iterations on one side of a ridge and crossed during the last 40k. That single excursion is what the γ_T[2] / ψ_T[1] between-chain excess is. It is not a second mode: chain 2 ends where chains 1 and 3 are.

## 3. Autocorrelation

Integrated autocorrelation time (iterations per effective draw) and the lag at which the autocorrelation falls below 0.1:

| | τ_int | acf < 0.1 at lag | | τ_int | acf < 0.1 at lag |
|---|---|---|---|---|---|
| φ[1] | **2041** | 2231 | γ_T[2] | 741 | 588 |
| σ_δT | 930 | 784 | ψ_T[1] | 774 | 866 |
| μ_bT | 898 | 952 | β | 589 | 742 |
| γ_T[1] | 882 | 1113 | μ_b | 825 | 1337 |
| *γ[2]* | 170 | 253 | *θ[4]* | 139 | 269 |
| *ψ[1]* | 433 | 685 | *ν_T[1]* | 2 | 2 |

The slow set needs 600–2,000 iterations per independent draw; the amyloid analogues need 140–430. With 105k retained draws per chain that is ESS 150–530, and R-hat ≈ 1 + c/ESS lands at 1.01–1.02. To bring all nine below 1.01 by length alone would need ESS ≈ 400 for φ[1], i.e. roughly 2.7× the current retained draws: ~400k iterations per chain, about two more days.

## 4. Tails

Excess kurtosis 0.0–0.2 for eight of the nine (γ_T[1] 0.85 with mild left skew). Not heavy-tailed; the gap between rank-normalised R-hat (1.019) and classic PSRF (1.010–1.016) comes from the within-chain drift that split-R-hat penalises, not from the marginals.

## 5. The ridge structure

Correlations among the nine and the principal components of their standardised draws:

- **Tau onset ridge**: μ_bT ~ γ_T[1] +0.57, γ_T[1] ~ γ_T[2] −0.51, μ_bT ~ γ_T[2] −0.39. A higher tau floor shortens the sub-threshold distance, so departure can be later (γ_T[1] up) and the APOE4 timing contrast re-balances.
- **Coupling ridge**: β ~ σ_δT +0.48, β ~ γ_T[1] −0.30, β ~ ψ_T[1] −0.24.
- **Amyloid floor**: φ[1] ~ μ_b −0.44.
- Timing–speed (APOE4): γ_T[2] ~ ψ_T[1] +0.37, moderate, same as the amyloid side where both mix well.

PC1 (28% of variance) = (γ_T[2] +0.50, φ[1] +0.30, β +0.24) against (γ_T[1] −0.53, μ_bT −0.49): **ESS 84, R-hat 1.044**. PC2–4: ESS 440–620. One direction, mixing the tau onset parameters with the amyloid susceptibility intercept, carries the slowness; it is the direction chain 2 crossed in §2.

How much of each parameter's spread is ridge: marginal sd ÷ sd conditional on all other population parameters is γ_T[1] 2.2, β 2.0, σ_δT 1.8, μ_b 1.8, μ_bT 1.6, φ[1] 1.4, γ_T[2] 1.4, ψ_T[1] 1.2.

## 6. The mechanism: narrow full conditionals, slowly moving subjects

Per sweep, the samplers move each slow parameter by only 0.13–0.28 of its posterior sd (β 0.13, ψ_T[1] 0.13, σ_δT 0.15, γ_T[1] 0.16, μ_bT 0.19; references 0.2–1.0), never rejecting. A step of 0.15 sd in a random walk would give τ_int ≈ 40, not 900; so the within-block step is not the limit. The limit is that each of these parameters has a full conditional, given the 1,101 subjects' coordinates, that is *far* narrower than its marginal: μ_bT's marginal sd is 0.0047, its sd given the other population parameters 0.0029, and its per-sweep movement 0.0009. Conditionally on where the subjects currently sit, the floor is nearly determined; it can only move as fast as the subject configuration moves.

The subject configuration moves slowly by construction. 424 of the 1,101 subjects have amyloid locations that straddle the floor (below it in 5–95% of draws): those are the subjects whose latent-vs-departed status is genuinely uncertain, and each is updated by its own 2- or 4-dimensional random-walk block. The count of below-floor subjects has ESS 573 of 15,750 thinned draws; it correlates +0.61 with μ_b and +0.44 (tau count) with μ_bT. φ[1] is the extreme case: its conditional on the population parameters leaves 50% of its variance to the subjects' state responsibilities, and those change one subject at a time. It takes the largest per-sweep steps of all (0.5 sd, because its conditional given a fixed subject configuration is wide) and yet has the longest autocorrelation time (2,041), because the centre of that conditional is set by the slowly re-sorting subjects.

This is the classic centred-hierarchical-model geometry: hyperparameters whose conditionals are narrow relative to their marginals alternate with a high-dimensional latent layer, and the Gibbs-style alternation random-walks along the ridge. It is not multimodality, not heavy tails, and (apart from the chain-2 excursion on the tau onset ridge) not disagreement between chains.

## 7. What would and would not help

| option | effect on the nine | cost |
|---|---|---|
| **Run longer** (≈400k/chain) | all < 1.01 by ESS alone; no estimate changes | ~2 days; the honest baseline |
| **Re-block** γ_T[2:3], ψ_T[1:2] with γ_T[1] | removes part of the tau-onset ridge inside one block; τ_int maybe 1.3–1.5× better; does not touch the subject-coupling bottleneck | cheap, modest |
| **Interweave the floors** (ASIS): alternate the current parameterisation with one where subject locations are defined *relative* to μ_b, μ_bT (x̃ − μ_b), so a floor move carries every subject with it | directly widens the conditionals of μ_b, μ_bT, and through them γ_T[1], φ[1]; the standard cure for this geometry | a new sampler for two scalars; moderate |
| **Marginalise more** (integrate subject coordinates out of the φ[1] update) | would fix φ[1] outright | not available in closed form here |
| More chains | lowers Monte Carlo error, not R-hat per chain | no |
| Different burn-in | 30% is already justified by stationarity | no |

The nine parameters' posterior means are stable across 50k, 75k and 150k to within 0.2 sd, and the quantities the paper reports (β, ρ, cure fractions, state probabilities, positivity ages) have ESS ≥ 500 or are functions of parameters that do. The remaining 1.01–1.02 is Monte Carlo noise on nuisance-scale parameters, and the reviewers' bar is met at < 1.02; getting to < 1.01 everywhere is a matter of either compute or the interweaving sampler, not of a defect in the model or the current blocking.
