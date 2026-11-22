# Pre-registration: asassn14ko-clock v0.1.0

Fixed on 2026-11-21, before any ASAS-SN or TESS data after 2021 were downloaded or looked at
by us. Any change after this point is listed under **Deviations** with its reason.

## Question

ASASSN-14ko (nucleus of ESO 253-G003, z = 0.042) flares every ~115 days, and its flares
come a little earlier each cycle. Payne et al. (2021, ApJ 910, 125) found the period
derivative; Payne et al. (2023, ApJ 951, 134; arXiv:2206.11278, posted June 2022) refit it as
Pdot = dP/dt = -0.0026 +- 0.0006 from flares up to November 2021 (cycle n = 24) and forecast
the next five peaks (2022-03-01 to 2023-05-04). In the code this fit is called PAYNE22. Models in which
a star on a ~115-day orbit loses mass at each pericentre (a repeating partial tidal disruption)
struggle to explain a decay this fast: gravitational waves are far too weak (Cufari et al. 2022)
and disc drag or tidal heating has been proposed (Linial & Quataert 2024; Bandopadhyay et al.
2024; Yao et al. 2025). Huang et al. (2025, ApJ) followed the flares to August 2024 with Swift
but did not update the timing model. **No published timing analysis uses the flares after
November 2021.** We ask:

1. Did the 2022-2023 flares arrive when Payne et al. (2023) forecast?
2. Over 2019-2026, is the period still shrinking, and at a constant rate, or is the decay
   speeding up or slowing down?
3. Are the flares fading (the star running out of mass)?

## Data (all public, no login)

- ASAS-SN Sky Patrol v2 g-band light curve of the catalogued source nearest the target
  (within 10", else 30"), through the public `pyasassn` client (Hart et al. 2023).
  Points with quality != "G", non-positive errors, or errors > 3x the median are dropped.
  Each camera's median-minus-global-median flux offset is subtracted (cameras with < 20 points
  dropped). Band = the archive's filter column when present; otherwise g = MJD >= 58450.
  Inverse-variance means in 1-day bins.
- TESS full-frame-image cutouts (MAST TESSCut, 15x15 px), 3x3-px aperture minus the median of
  pixels >= 5 px away, QUALITY == 0, 0.25-day bins. Cross-check only.

## Peak measurement (ko/peaks.py)

Cycle by cycle from n = 0: search +-25 d around the predicted peak (the Payne+2023 ephemeris
until five peaks are measured, then a constant-Pdot fit to the peaks measured so far, with 1 d
added in quadrature). Provisional peak = brightest 3-bin running mean in the window. Fit a
weighted polynomial (degree 4 / 3 / 2 for >= 10 / >= 7 / >= 5 bins) to bins from 12 d before to
15 d after it; peak = polynomial maximum within +-8 d (0.05 d grid). Rejected if: < 8 bins in
the window; < 2 bins on either side; maximum at the grid edge; amplitude above baseline
(median flux 45-15 d before the peak) < 5x the median bin error. Error = standard deviation of
300 bootstrap refits (floor 0.1 d).

## Gates

- **G1 (must pass, else no grading):** our g-band peaks reproduce the published ASAS-SN g peaks
  of cycles 15, 17, 18, 21, 22, 24 (Payne+2021/2023) within max(2 d, 3 sigma combined) for at
  least 5 of the 6.
- **G2 (TESS pipeline only):** our TESS peaks for cycles 14 (sector 4) and 21 (sector 33) are
  within 1.5 d of the published TESS peaks. If G2 fails, TESS results are not used.
- At least 4 measured cycles with n >= 25, else outcome INSUFFICIENT_NEW_DATA.

## Tests

- **T1 (forecast):** O - C of the five Payne+2023 forecasts (n = 25-29); "on time" if
  |O - C| <= 5 d. Also the O - C of every later cycle against the Payne+2023 ephemeris.
- **T2 (primary, period evolution):** all our measured g-band cycles. Models: linear, steady
  (constant Pdot: t = a + b m + c m^2 + d m^3 with d tied to 2c^2/3b, m = n - n_ref;
  Pdot = 2c/b), cubic (d free). One jitter, added in quadrature, chosen so the cubic fit has
  reduced chi^2 = 1 (0 if already <= 1), used for all models. z = (d - 2c^2/3b) / its error.
  Outcome:
  - NO_DETECTED_DECAY if |Pdot| < 3 sigma;
  - PERIOD_GROWING if Pdot > 0 at >= 3 sigma;
  - DECAY_ACCELERATING if Pdot < 0 at >= 3 sigma and z <= -3;
  - DECAY_SLOWING if Pdot < 0 at >= 3 sigma and z >= +3;
  - DECAY_STEADY otherwise.
  Reported with it: Pdot vs Payne+2023 (sigma), period at the last cycle, period change in hours
  per orbit, and the timescale P/|Pdot|.
- **T2b (secondary):** same fit on published peaks for n <= 24 plus ours for n >= 25.
- **T3 (secondary):** Spearman rank correlation of flare amplitude (mJy above baseline) with n,
  two-sided permutation p (10,000 shuffles, seed 20261121): FADING if rho < 0 and p < 0.01,
  BRIGHTENING if rho > 0 and p < 0.01, else NO_TREND. Aperture photometry includes the host
  and its AGN, so this is indicative only.

## Exploratory (not graded)

AT2023uqm: Wang et al. (2025/2026, Science Bulletin; arXiv:2510.26561) forecast its sixth
flare's first peak near MJD 61280 (2026-08-28). If Sky Patrol has a catalogued source within
5", we report whether >= 2 g-band bins in MJD 61250-61330 exceed the MJD 60900-61230 median by
> 4 sigma. ASAS-SN is likely too shallow for this z = 0.238 source (flare 5 peaked near
r = 18.8), so a non-detection means little; this is a pointer for a proper test, not a test.

## What each outcome would mean

- DECAY_STEADY: the orbit has kept shrinking at one rate for 12 years; any model must produce
  a constant Pdot.
- DECAY_ACCELERATING / SLOWING: first evidence that the decay rate itself is changing, which
  separates the models (drag and mass loss depend on the orbit and the star's remaining mass).
- NO_DETECTED_DECAY: the published decay does not hold up on a longer baseline.

## Deviations

**D1 (2026-11-22, data access only; code 0.1.1, `ko/asassn_io.py`).** The v0.1.0 job ran on
2026-11-22 and ended `NO_DATA`: Sky Patrol's cone download returned "Pulled 1 of 2" light
curves and 0 rows for the nearest catalogued source (asas_sn_id 498216284073, 0.9" away), while
Sky Patrol announced on 2026-11-03 that it is rebuilding its catalogues and light curves. No
ASAS-SN flux was read. Change: for each catalogued source within the radius, in order of
separation, use the cone download's rows if there are any, else request that source's light
curve by ID (two tries); use the first source that has rows. The second candidate is a
Swift-catalogue entry 7.9" away, about one ASAS-SN pixel, so it samples the same nucleus. Which
source was used is written to results.json. Gate G1 (reproducing the six published peaks) still
decides whether anything is graded. Nothing else changes.

Note on v0.1.0's TESS output (cross-check only, not graded): G2 did not pass because sector 4
was not measured, so by the rule above TESS results are not used. For full disclosure: before
writing D1 we saw v0.1.0's TESS summary, which reproduced the published sector-33 peak (1.2 d
early) and reported no flare above baseline in sector 98 near the extrapolated ephemeris time
(MJD ~61020). Because G2 failed, that is not a result; D1 does not touch the TESS code.
