# Amendment 2026-11-27: does ASASSN-14ko's flare clock wobble as the orbit precesses?

Frozen 2026-11-27 (about 18:40 UTC) by AI agents working in MV2's name (Munim, Inc.), **before the
periodogram, false-alarm probability or any upper limit below was computed**. Licence: Apache-2.0
(code), CC BY 4.0 (text).

**Honest status of the blinding:** the 29 flare times used here are already public in this repository
(results/v0.1.0/results.json) and their residuals from the steady-decay fit were plotted in our own
draft. So this is not a blind test of new data. What is frozen is the analysis: the model, grid,
statistics, thresholds and outcome names, written down before running them once.

## The gap this closes

Linial & Quataert (2024, MNRAS 527, 4317; arXiv:2309.15849) explain ASASSN-14ko's shrinking period by
hydrodynamic drag as the star crosses the black hole's accretion disc twice per orbit. In their model
the drag, and therefore the period derivative, is modulated as the orbit precesses, which leaves an
oscillatory term in the flare times with period T_prec/2 (their Fig. 1, about 16 orbits for their
example orbit), with T_prec/P0 = a(1-e^2)/(3 R_g) = (r_p/R_g)(1+e)/3, about 13 (r_p/2r_t) orbits
(their eq. 14 and sec. 4.3). They "encourage further analysis of the flare timing data, and its
deviations from constant Pdot". Their own check of Pdot used flares to 2021 (Payne et al. 2022/2023).

Nobody has published a search for this oscillation. Our v0.1.0 fit (29 flares, 2014-2024) gives
Pdot = -0.00197 +- 0.00014 and a cubic term consistent with steady decay, but does not test a
periodic term.

## Data

`timings_29.csv`: the 29 g-band peak times used in t2_full_baseline of v0.1.0 (21 published peaks to
2021, 8 of ours from 2022-2024), copied from results/v0.1.0/results.json; the job cross-checks every
value against that file. The TESS peak of cycle 38 stays out (its pre-set TESS gate failed in v0.1.0).

## Method (precession.py, run once)

- m = n - 17. Cubic-only model t = a + b m + c m^2 + d m^3 (equal weights; the scatter is dominated by
  the 2.5-day extra jitter, not by the quoted errors). Pdot = 2c/b.
- Wobble model: cubic + A cos(2 pi m/N) + B sin(2 pi m/N), amplitude K = sqrt(A^2 + B^2).
- Grid: frequency 1/N from 1/34 to 1/4 per orbit, step 1/340 (oversampling 10 over the 34-orbit span).
- Statistic: F = ((RSS0 - RSS1)/2) / (RSS1/(n-6)); best = max over the grid.
- Global false-alarm probability: fraction of 2,000 white-noise data sets on the same 29 epochs whose
  max F reaches the observed one ((hits+1)/(2,001)); seed 20261127.
- 95% upper limit K95(N) for N in {6, 8, 10, 13, 16, 20, 25, 30} orbits: the true amplitude at which
  95% of 500 simulated data sets (random phase, Gaussian noise with the cubic fit's rms) return a
  fitted amplitude at least the observed one (bisection, common random numbers).
- Modulation fraction: if Pdot(m) = Pdot0 (1 + f cos(2 pi m/N)), then K = f |Pdot0| P (N/2 pi)^2, so
  f95 = K95 (2 pi/N)^2 / (|Pdot0| P), with Pdot0 and P from the cubic-only fit.
- Mapping to the pericentre (Linial & Quataert's scaling, not ours): N = T_prec/(2 P0), so
  r_p/2r_t = N/6.5.

## Outcomes

1. **WOBBLE_DETECTED** if global FAP < 0.01: report N, K, phase; a claim that the drag (or something
   else) is modulated on that period, to be checked against the next flares.
2. **WOBBLE_HINT** if 0.01 <= FAP < 0.05.
3. **NO_WOBBLE** otherwise. Then report K95 and f95 on the grid and **N\***, the shortest grid period
   such that f95 < 1 at it and at every longer grid period: order-unity modulation of the drag on
   periods >= N\* orbits is excluded, which in Linial & Quataert's scaling means r_p/2r_t < N\*/6.5 if
   their drag is modulated by order unity.

Secondary, not graded: the same search without flare 23 (the largest residual, +8.9 d).

## Known limits, stated in advance

- 29 epochs with gaps; a wobble with N near the 34-orbit span is partly absorbed by the cubic term,
  so limits at N = 25-30 are weaker than white-noise scaling suggests.
- Flare shape changes (bumps, rebrightenings; Huang et al. 2023, 2025) move peak times by days and act
  as extra noise; a real wobble smaller than about 1-2 days cannot be seen.
- The mapping from N to r_p, and the size of the modulation f, depend on the disc's structure and the
  orbit's inclination, which are not known; f ~ 1 is a reference value, not a prediction.
- Light-travel-time and collision-timing variants (Linial & Quataert sec. 5) predict other signatures;
  they already argued those variants are disfavoured.

## Run

`python -m unittest -v` then `python precession.py --data timings_29.csv --results-json
results/v0.1.0/results.json --out-dir /out`. The output will be committed unchanged as
results/precession-v0.1.0/.
