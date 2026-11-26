# Grading of asassn14ko-clock v0.1.0 (code 0.1.1)

Graded 2026-11-26 (UTC) against [PREREGISTRATION.md](PREREGISTRATION.md) (fixed 2026-11-21,
deviation D1 of 2026-11-22). Job `asassn14ko-clock-v2` ran 2026-11-26 07:22-07:47 UTC; all files
are in `results/v0.1.0/`. MV2, Munim, Inc. Made with AI agents working in MV2's name; no human
expert has checked it yet.

## Outcome: DECAY_STEADY

## Data

- ASAS-SN Sky Patrol v2, g band. The nearest catalogued source (asas_sn_id 498216284073, 0.9") again
  returned 0 rows (cone and by-ID, two tries), so by deviation D1 the job used the next source,
  asas_sn_id 661431833080 (a Swift-catalogue entry 7.9" away, about one ASAS-SN pixel, same nucleus):
  2,429 rows, 1,667 clean g points in 1,162 daily bins, last g point MJD 60708.97 (2025-01-28).
- Gate G1 (does our pipeline reproduce the published peaks?): **pass, 6 of 6** (cycles 15, 17, 18,
  21, 22, 24 within max(2 d, 3 sigma); differences 0.73-1.93 d). This also validates the fallback
  source.
- New cycles measured (n >= 25): **8** (25, 27, 28, 29, 30, 31, 32, 34). Cycles 26 and 33 fell in the
  seasonal gap (sparse). Gate (>= 4) passed.

## T1: did the Payne et al. (2023) forecasts come true?

| Cycle | Forecast (MJD) | Observed (MJD) | O - C (d) | On time (|O-C| <= 5 d) |
|---|---|---|---|---|
| 25 | 59639.7 | 59635.89 +- 0.84 | -3.81 | yes |
| 26 | 59747.5 | not measured (seasonal gap) | | |
| 27 | 59855.9 | 59855.74 +- 0.67 | -0.16 | yes |
| 28 | 59962.1 | 59963.40 +- 0.63 | +1.30 | yes |
| 29 | 60069.0 | 60082.05 +- 4.16 | +13.05 | **no** |

3 of 4 measured forecasts on time. Against the Payne et al. (2023) ephemeris, every later flare
came late: O - C = +12.5 (n = 29), +13.8 (30), +8.1 (31), +6.7 (32), +14.0 d (34). The ephemeris's
period derivative (-0.0026) shortens the period too fast.

## T2 (primary): period evolution, our g-band peaks only (n = 12-34, 18 cycles)

- Pdot = **-0.00169 +- 0.00041** (4.1 sigma), so the period is still shrinking.
- Cubic excess z = **+1.58** (|z| < 3): no significant acceleration or slowing. -> **DECAY_STEADY**.
- Jitter 3.32 d (added in quadrature so the cubic fit has reduced chi^2 = 1). chi^2: linear 33.1
  (16 dof), steady 16.5 (15), cubic 14.0 (14).
- Versus Payne et al. (2023) Pdot = -0.0026 +- 0.0006: 1.25 sigma.
- Period 108.99 d at n = 23; **106.98 d at n = 34** (late 2024); -4.3 h per orbit; P/|Pdot| = 174 yr.

## T2b (secondary): published peaks n <= 24 plus ours n >= 25 (29 cycles, 2014-2024)

- Pdot = **-0.00197 +- 0.00014** (14.3 sigma); cubic z = +1.01; jitter 2.54 d. -> DECAY_STEADY.
- chi^2: linear 227.7 (27 dof), steady 26.0 (26), cubic 25.0 (25).
- Versus Payne et al. (2023): 1.0 sigma. Versus Payne et al. (2021, -0.0017 +- 0.0003): within 1 sigma.
- Period 110.42 d at n = 17, **106.78 d at n = 34**; -5.0 h per orbit; P/|Pdot| = 148 yr.

## T3 (secondary, indicative): are the flares fading?

Spearman rho = -0.29, permutation p = 0.235 (n = 18): **NO_TREND**.

## TESS (not graded)

Gate G2 failed again because the sector-4 peak (cycle 14) could not be measured, so by the frozen
rule TESS results are not used. For disclosure only: the sector-33 peak was 1.23 d before the
published TESS peak, and sector 98 shows a flare peak at MJD 61042.12 (cycle 38, 2026-01-02), 3.1 d
after the primary steady fit's prediction.

## Forward forecast (new, written 2026-11-26 before any later data were looked at)

From the T2b steady fit in `results/v0.1.0/results.json` (t = a + b m + c m^2 + d m^3, m = n - 17,
a = 58765.7718, b = 110.41875, c = -0.10885594, d = 2c^2/(3b) = 7.1544e-5), with the T2 fit
(m = n - 23, a = 59424.7457, b = 108.98762, c = -0.09196616) as a cross-check:

| Cycle | T2b prediction (MJD, date) | T2 prediction (MJD) | Note |
|---|---|---|---|
| 41 | 61354.1 (2026-11-10) | 61357.0 | just passed |
| 42 | **61459.3 (2027-02-23)** | 61462.7 | next well-placed flare |
| 43 | 61564.3 (2027-06-08) | 61568.1 | near solar conjunction |
| 44 | 61669.1 (2027-09-21) | 61673.4 | |

Expected scatter: the fitted jitter is 2.5-3.3 d, so a peak within about +-7 d of these dates is
consistent with steady decay; a peak 10 d or more early would point to accelerating decay.

## What it means

Over twelve years and 29 flares the period of ASASSN-14ko has shrunk at one constant rate,
dP/dt = -0.0020 +- 0.0001, about four times more precisely measured than the last published value
(-0.0026 +- 0.0006, flares to 2021). The data show no sign that the decay is speeding up or slowing
down, and the newer flares do not follow the steeper published rate: since 2023 they have come 7-14 d
later than the 2023 ephemeris predicted. Any model of the decay (disc drag, tidal heating, mass
transfer) must give a nearly constant dP/dt of about -0.002 over at least 30 orbits.

## Limits

- Aperture photometry includes the host galaxy and its AGN; peak times can be pulled by flare
  shape changes (bumps and rebrightenings, Huang et al. 2023, 2025). The jitter term absorbs this
  scatter but does not model it.
- The ASAS-SN light curve used is the Swift-catalogue source 7.9" from the nucleus (deviation D1);
  G1 shows it reproduces the published peaks.
- ASAS-SN Sky Patrol v2 data end in January 2025; flares from 2025-2026 are not in this test.

## Credit

Payne et al. 2021 (ApJ 910, 125), 2022 (ApJ 926, 142), 2023 (ApJ 951, 134; arXiv:2206.11278);
Huang et al. 2023 (arXiv:2310.03236), 2025 (arXiv:2506.15039). ASAS-SN Sky Patrol v2 (Hart et al.
2023); TESS and TESSCut (Brasseur et al. 2019).
