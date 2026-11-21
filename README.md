# asassn14ko-clock

**Is the star in ASASSN-14ko still spiralling in on schedule?**

ASASSN-14ko is a galaxy nucleus (ESO 253-G003) that flares every ~115 days, most likely
because a star on a tight, eccentric orbit around the supermassive black hole loses mass
at each close approach. The flares arrive a little earlier every cycle. The last published timing
fit (Payne et al. 2023, ApJ 951, 134; arXiv:2206.11278) used flares up to November 2021 and
found the period shrinking at dP/dt = -0.0026 +- 0.0006, far too fast for gravitational
waves. No timing analysis since has used the roughly 15 flares from 2022 to 2026.

This repository is an open, **pre-registered** test (see [PREREGISTRATION.md](PREREGISTRATION.md))
of three questions, using only public data (ASAS-SN Sky Patrol v2 and TESS):

1. Did the five forecast peaks of 2022-2023 arrive on time?
2. Is the period still shrinking, and is the decay steady, speeding up or slowing down?
3. Are the flares fading?

The rule, gates and outcome classes were fixed before any post-2021 data were downloaded.

## Layout

- `ko/ephem.py`: published ephemeris, forecasts and peak times (with references)
- `ko/fit.py`: dependency-free weighted least squares, timing models, rank statistics
- `ko/peaks.py`: per-cycle peak measurement with bootstrap errors
- `ko/grade.py`: gates and outcome classes
- `ko/asassn_io.py`, `ko/tess_io.py`: data access (pyasassn client; MAST TESSCut)
- `job.py`: runs everything; writes `results.json`, `peaks_g.csv` and binned light curves
- `test_ko.py`: unit tests (`python -m unittest -v test_ko`)
- `notes/`: prior-art check

## Run

```
pip install -r requirements.txt
python -m unittest -v test_ko
python job.py --data-dir data --out-dir out
```

## Credit

Setup and data: Payne et al. 2021 (ApJ 910, 125), 2022 (ApJ 926, 142), 2023 (ApJ 951, 134);
Huang et al. 2023 (arXiv:2310.03236) and 2025 (arXiv:2506.15039); ASAS-SN (Shappee et al.
2014; Kochanek et al. 2017; Sky Patrol v2, Hart et al. 2023); TESS (Ricker et al. 2015) and
TESSCut (Brasseur et al. 2019). Theory that makes the question matter: Cufari et al. 2022;
Linial & Quataert 2024; Bandopadhyay et al. 2024; Yao et al. 2025.

## About

By MV2 of Munim, Inc. (mv2@munimtech.com). The code and analysis were made with AI agents
working in his name; no human expert has checked them yet. Licence: Apache-2.0.
