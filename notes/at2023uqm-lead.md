# Lead: did AT2023uqm's sixth flare arrive in August 2026? (checked 2026-11-21)

AT2023uqm (ZTF23abikzfr, ATLAS23txi; RA 347.66256, Dec +34.90301; z = 0.238) is a periodic
nuclear transient with five flares about 527 days apart whose energy grows by ~1.7-1.8x per
flare. Wang et al. (arXiv:2510.26561; Science Bulletin 71, 4645, 2026) fit the first-peak times
as t = 58657.80 + 524.45 (n - 1) (MJD) and forecast the sixth flare's first peak near
MJD 61280 (2026-08-28). They wrote that the star "may be completely disrupted after only one or
two more flares", and that a giant-star donor would make the flares peak and then decline.

What public data cover the forecast (MJD ~61250-61330), as checked on 2026-11-21:
- **ZTF public alerts (ALeRCE API):** the object's last alert is MJD 60960 (2025-10-12).
  No alert of any object within 1 degree has lastmjd later than MJD 61159 (2026-05-10),
  although ALeRCE carries ZTF alerts elsewhere up to MJD 61362 (2026-11-18). So the field was
  not covered by public ZTF alerts during the forecast window: the missing alerts are
  **not** evidence of a missing flare.
- **Swift:** 51 observations, the last at MJD 61079 (2026-02-10); none in the window.
- **Rubin/LSST:** Dec +35 is outside the main survey footprint.
- **ATLAS forced photometry** covers it (the discovery survey of flare 4) but needs a free
  account (fallingstar-data.com/forcedphot).
- **ASAS-SN Sky Patrol:** probably too shallow; the asassn14ko-clock job looks anyway
  (exploratory only).
- No paper, ATel or AstroNote reporting the sixth flare was found (arXiv API all:AT2023uqm
  shows only the discovery paper, v3 of 2026-10-09).

Next step: with ATLAS access, pre-register a test (flare detected? first-peak time vs 61280;
energy ratio to flare 5 vs the 1.67-1.79 growth) before reading the photometry.
