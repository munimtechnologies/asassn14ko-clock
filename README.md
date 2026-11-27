# asassn14ko-clock

**Is the star in ASASSN-14ko still spiralling in on schedule?**

ASASSN-14ko is a galaxy nucleus (ESO 253-G003) that flares every ~115 days, most likely
because a star on a tight, eccentric orbit around the supermassive black hole loses mass
at each close approach. The flares arrive a little earlier every cycle. The last published timing
fit (Payne et al. 2023, ApJ 951, 134; arXiv:2206.11278) used flares up to November 2021 and
found the period shrinking at dP/dt = -0.0026 +- 0.0006, far too fast for gravitational
waves.

## Result (v0.1.0, graded 2026-11-26): DECAY_STEADY

With eight new ASAS-SN flare peaks from 2022-2024 (after our pipeline reproduced all six published
peaks it was checked against), the period has shrunk at one constant rate for twelve years:
**dP/dt = -0.00197 +- 0.00014** from 29 flares (2014-2024), 14 sigma from zero, with no sign of
acceleration or slowing (cubic term z = +1.0). Our peaks alone give -0.00169 +- 0.00041 (4.1
sigma). Three of the four measurable 2022-2023 forecasts came true; since 2023 the flares arrive
7-14 days later than the 2023 ephemeris, whose steeper rate overshoots. The period was 106.8 d in
late 2024. The next well-placed flare is predicted near 2027-02-23 (MJD 61459, +-7 d). Full
grading: [GRADING.md](GRADING.md); files in `results/v0.1.0/`.

The rule, gates and outcome classes were fixed before any post-2021 data were downloaded
([PREREGISTRATION.md](PREREGISTRATION.md)).

## Open test (frozen 2026-11-27): does the clock wobble as the orbit precesses?

Linial & Quataert (2024, arXiv:2309.15849) explain the decay by drag as the star crosses the
black hole's gas disc, and expect that drag to be modulated as the orbit precesses, leaving an
oscillation in the flare times (period T_prec/2, about 16 orbits in their example). They asked for
exactly this analysis of the timing data. `precession.py` searches the 29 flare times for that
wobble (global false-alarm probability from simulations) and, if there is none, sets 95% limits on
its amplitude and on the fractional modulation of dP/dt as a function of period. Rules, outcomes and
an honest note on blinding: [AMENDMENT-2026-11-27-precession.md](AMENDMENT-2026-11-27-precession.md).
Results will appear in `results/precession-v0.1.0/`.

## Layout

- `ko/ephem.py`: published ephemeris, forecasts and peak times (with references)
- `ko/fit.py`: dependency-free weighted least squares, timing models, rank statistics
- `ko/peaks.py`: per-cycle peak measurement with bootstrap errors
- `ko/grade.py`: gates and outcome classes
- `ko/asassn_io.py`, `ko/tess_io.py`: data access (pyasassn client; MAST TESSCut)
- `job.py`: runs everything; writes `results.json`, `peaks_g.csv` and binned light curves
- `test_ko.py`: unit tests (`python -m unittest -v test_ko`)
- `precession.py`, `test_precession.py`, `timings_29.csv`: the precession-wobble search
  (standard library only)
- `notes/`: prior-art checks

## Run

```
pip install -r requirements.txt
python -m unittest -v test_ko
python job.py --data-dir data --out-dir out
python -m unittest -v test_precession
python precession.py --data timings_29.csv --results-json results/v0.1.0/results.json --out-dir out
```

## Credit

Setup and data: Payne et al. 2021 (ApJ 910, 125), 2022 (ApJ 926, 142), 2023 (ApJ 951, 134);
Huang et al. 2023 (arXiv:2310.03236) and 2025 (arXiv:2506.15039); ASAS-SN (Shappee et al.
2014; Kochanek et al. 2017; Sky Patrol v2, Hart et al. 2023); TESS (Ricker et al. 2015) and
TESSCut (Brasseur et al. 2019). Theory that makes the question matter: Cufari et al. 2022;
Linial & Quataert 2024; Bandopadhyay et al. 2024; Yao et al. 2025.

## About

By MV2 of Munim, Inc. (mv2@munimtech.com). The code and analysis were made with AI agents
working in his name; no human expert has checked them yet. Licence: Apache-2.0 (code); text and
results CC BY 4.0.
