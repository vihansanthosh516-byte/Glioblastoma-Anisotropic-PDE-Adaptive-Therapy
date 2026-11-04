# Track C model: four exact statements and what they imply (2026-11-04)

Status: Stage 2 test B of `paper/trackc_restart_plan_2026-11-03.md`. Claude's own derivations. Numeric checks: `tests/test_trackc_lemmas.py` (13 tests, scratch-copy run). A read-only Haiku proof-checker reviews the algebra below (result recorded at the end). SIMULATION model only. These statements describe the assumed equations in `src/trackc_sim.py`, not any tumour.

## Setting

Two cell types, sensitive S and resistant R. Total burden B = S + R, resistant fraction f = R / B, kill rate k = m r, drug indicator u (0 or 1). Equal competition (a_rs = a_sr = 1, the default):

    dS/dt = r S (1 - B) - k u S                      - mu u S + nu (1 - u) R
    dR/dt = r (1 - c) R (1 - B) - eps k u R          + mu u S - nu (1 - u) R

The last terms are the optional switching (mu, while on drug) and reversion (nu, while off drug). They cancel in the sum, so they do not change dB/dt directly.

## Lemma 1 (the resistant fraction)

    df/dt = f (1 - f) [ k (1 - eps) u - c r (1 - B) ] + mu u (1 - f) - nu (1 - u) f
    d logit(f) / dt = k (1 - eps) u - c r (1 - B) + mu u / f - nu (1 - u) / (1 - f)

Proof. Sum the two equations: dB/dt = r (1 - B) (S + (1 - c) R) - k u (S + eps R) = B [ r (1 - B)(1 - c f) - k u (1 - f + eps f) ]. Then df/dt = (dR/dt)/B - f (dB/dt)/B. The growth part gives r (1 - B) f [ (1 - c) - (1 - c f) ] = - c r (1 - B) f (1 - f). The kill part gives k u f [ - eps + (1 - f + eps f) ] = k u (1 - eps) f (1 - f). The switching part gives mu u S / B - nu (1 - u) R / B = mu u (1 - f) - nu (1 - u) f. Dividing by f (1 - f) gives the logit form; this needs 0 < f < 1, that is S > 0 and R > 0. QED.

Numeric check: pointwise to 1e-10 at 4,000 random states, and in integrated form to 1e-7 for four parameter sets (`test_lemma1_derivative_identity_at_random_states`, `test_lemma1_integrated_form`).

## Lemma 2 (timing does not change resistance when there is no cost)

If c = mu = nu = 0, then logit f(t) = logit f0 + k (1 - eps) D(t), where D(t) = the integral of u, the drug-days so far.

Proof. Put c = mu = nu = 0 in Lemma 1 and integrate. QED.

So, with c = mu = nu = 0, two schedules with the same drug-days leave the same resistant fraction, whatever the timing and whatever B did in between. With a cost (c > 0) the fraction also falls by c r times the integral of (1 - B) dt, which depends on the path. With switching, it also moves by mu times the integral of u/f and nu times the integral of (1-u)/(1-f).

Numeric check: three schedules with 6 courses, different timing, same final resistant fraction to 1e-6 (`test_lemma2_resistant_fraction_depends_only_on_drug_days`).

## Lemma 3 (untreated time to progression)

With no drug and c = 0, B is logistic whatever f is: dB/dt = r B (1 - B). The time for B to rise from B0 to 1.2 B0 (needs 1.2 B0 < 1) is

    TTP = (1 / r) ln[ 1.2 (1 - B0) / (1 - 1.2 B0) ].

Proof. logit B grows at rate r, so r TTP = logit(1.2 B0) - logit(B0). QED. With c > 0 and no drug, dB/dt = r B (1 - B)(1 - c f), so d logit(B)/dt = r (1 - c f) is at most r at every moment. logit B therefore needs at least as long to climb the same distance, and the formula is a lower bound on TTP. The formula needs B0 < 5/6 so that 1.2 B0 < 1.

Numeric check: five (r, B0) pairs to 0.01 day (`test_lemma3_untreated_time_to_progression`).

## Proposition (no-crowding limit: front-loaded dosing is optimal)

Replace the crowding factor (1 - B) by 1 in the growth terms (no carrying capacity). Take c = mu = nu = 0. Allowed schedules: drug only on days 0-4 of 28-day cycles, at most 12 courses, the same r, k, eps, B0 and f0 for every schedule. Assume k > 0, eps >= 0 and 0 <= f0 <= 1. Progression is the first time with B >= 1.2 B0; if there is none, TTP is the horizon cap. Then front-loaded dosing (a course in every cycle from the first until the budget is spent) gives a time to progression at least as long as any other allowed schedule.

Proof. In this limit ln B(t) = ln B0 + r t - K(D(t)), with K(D) = integral from 0 to D of k (1 - (1 - eps) g(x)) dx and g(x) = logistic(logit f0 + k (1 - eps) x). Lemma 2 still holds in this limit: with c = 0 both populations grow at the same rate r, so growth cancels in the ratio R/B and only the kill part moves f. The integrand of K is positive when k > 0 and eps >= 0, so K is increasing in D. At every time t, front-loaded dosing has the largest drug-days D(t) of any allowed schedule: cycle j adds min(5, max(0, t - 28 j)) days, which does not increase with j, and a schedule uses at most 12 cycles. So ln B(t) is smallest under front-loaded dosing at every t. Whenever another schedule has reached 1.2 B0 at time t, front-loaded dosing has not reached it earlier, so its first crossing is no earlier. QED.

Checked by the independent proof-checker (below): the conclusion holds in this limit; it is a statement about the no-crowding variant of the equations, not about the full model.

Consequence for the full model (with crowding), by the contrapositive of the Proposition: any gain of a feedback schedule over front-loaded dosing must come from crowding (B not small compared with the carrying capacity 1), from the cost c, or from switching (mu, nu). This is a consequence of the algebra, not a measurement.

Evidence from saved and new runs (all SIMULATION, development data, exploratory):
- Brute force over every schedule with 4 courses in 10 cycles, c = 0, mu = nu = 0, fast growth: the best schedule beat front-loaded dosing by 1.3 d on average (max 4.2 d) when B0 was 0.05-0.10, and by 21.9 d on average (max 53 d) when B0 was 0.40-0.50 (16 scenarios per group). The gain was smaller at lower B0 in this run. Bounded versions of these checks are in `test_front_loaded_is_nearly_optimal_without_cost_when_burden_is_small` and `test_gain_from_perfect_timing_grows_with_starting_burden` (8 scenarios per group, run in a scratch copy: passed).
- Known-parameter oracle on 200 scenarios with 12 courses and 52 cycles (`output/trackc_headroom/headroom_per_scenario.csv`, one development draw, no uncertainty estimate): mean gain over front-loaded dosing 29 d in the lowest B0 quarter and 132 d in the highest. In a linear model on standardised parameters (R-squared 0.33) B0 had the largest standardised coefficient (+38.5 d per standard deviation), then growth rate r (-26.6), then cost c (+8.3). In this draw the mean gain was larger for high c at high B0 (156 d against 105 d) and similar at low B0 (38 d against 31 d); there was no interaction test.

Weakest points: the Proposition is exact only in the no-crowding limit; the brute-force horizon is short (10 cycles); the oracle numbers are lower bounds on the gain achievable with known parameters over the schedules searched (`src/trackc_oracle.py` schedules are achievable, not proven optimal), and the 1,456-day cap biases the means downward; the regression is linear on one development draw. B0 (the starting burden as a fraction of an assumed carrying capacity) and c (the cost of resistance) are assumptions with no measurement behind them in this project.

## Proof-checker result (read-only Haiku 5.5 checker, 2026-11-04; its Python was run in a scratch folder, not in the repository)

- Lemma 1: both formulas PASS; every line of the proof CORRECT. Pointwise check at 20,000 random states: largest difference 3e-16 (df/dt) and 6e-14 (d logit f/dt). Missing assumption 0 < f < 1: added.
- Lemma 2: PASS. RK4 at step 0.01 day, three 6-course schedules (contiguous, spaced, late): residual 3e-12.
- Lemma 3: PASS, five (r, B0) pairs within 2e-8 day. The lower-bound argument for c > 0 was asserted, not argued: argument and the condition B0 < 5/6 added.
- Proposition: conclusion holds in the no-crowding limit; no counter-example found. Corrections made: the reason Lemma 2 survives is c = 0 (both populations grow at the same rate), not that "the crowding factor cancels"; the claim is "at least as long", not "the longest"; hidden assumptions (same parameters across schedules, k > 0, eps >= 0, meaning of "first crossing", the no-crowding variant) are now stated. The checker notes that for eps < 0 the argument need not hold; the model uses eps in [0, 0.2].
- Full-model spot check (r = 0.1, k = 0.3, eps = 0.5, c = 0.2, mu = 0.004, nu = 0.01, 300 days): residual 3.5e-12.
- Evidence section: overclaims found ("cost matters mainly together with high B0", "lower bounds on what perfect knowledge can gain", "B0 the largest single term", "it should be small", "reproduced"). Rewritten above as plain descriptions of one development draw.
- The checker did not run `tests/test_trackc_lemmas.py`. Claude ran it in a scratch copy: 13 passed in 46 s.
