# Prior-art check, 2026-11-27 (~18:35 UTC): precession wobble in ASASSN-14ko flare times

Question: has anyone searched ASASSN-14ko's flare arrival times for the oscillation that
Linial & Quataert (2024) expect from apsidal precession of a disc-crossing orbit, or set limits on it?

Searches:
1. arXiv API `all:"ASASSN-14ko"`, newest first (18 results; newest arXiv:2611.21731, 2026-11-16).
   None of the 8 newest abstracts mentions precession or an oscillation in flare times
   (2603.04559, 2510.26561, 2506.15039, 2505.10611, 2410.02146, 2406.03675, 2402.18812 checked).
2. Web search "ASASSN-14ko flare timing apsidal precession oscillation residuals period derivative
   modulation": only Linial & Quataert 2024 (the prediction) and the Payne et al. timing papers
   (2021, 2022, 2023), whose fits stop at flares from 2021 and test a constant or piecewise Pdot.
3. Full HTML of Linial & Quataert 2024 (arXiv:2309.15849v2) read for the prediction: oscillation with
   period T_prec,ap/2 (about 16 orbits in their Fig. 1 example), T_prec/P0 = a(1-e^2)/(3 R_g) (eq. 14),
   and "We encourage further analysis of the flare timing data, and its deviations from constant Pdot".
4. Earlier checks for this repository (notes/scoop-check-2026-11-26*.md): no timing fit with flares
   after 2021 exists apart from ours.

Closest prior work and how this differs:
- Linial & Quataert 2024 (MNRAS 527, 4317): the prediction; no fit to data for the oscillation.
- Payne et al. 2023 (ApJ 951, 134): fits Pdot to flares to late 2021 and remarks the derivative may
  vary; no periodic-term search.

Conclusion: no published search for, or limit on, a periodic wobble in ASASSN-14ko's flare times.
