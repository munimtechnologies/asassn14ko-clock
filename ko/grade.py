# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Pre-registered gates and outcome classes (see PREREGISTRATION.md)."""
import math
import random

from . import ephem as E
from . import fit as F

ON_TIME_DAYS = 5.0       # T1: a Payne+2023 (arXiv 2022) forecast is "on time" if |O - C| <= 5 d
MIN_NEW_CYCLES = 4       # need >= 4 measured cycles with n >= 25
SIG = 3.0                # significance threshold for Pdot and the cubic term
TREND_P = 0.01           # T3 permutation p-value threshold
TESS_GATE_DAYS = 1.5


def gate_g1(measured):
    """measured: {n: (t, err)} from our ASAS-SN g-band pipeline."""
    pub = {n: (t, e) for n, t, e, _, _ in E.PUBLISHED_PEAKS}
    rows, n_meas, n_pass = [], 0, 0
    for n in E.GATE_CYCLES:
        tp, ep = pub[n]
        if n not in measured:
            rows.append({"n": n, "published": tp, "status": "not_measured"})
            continue
        t, e = measured[n]
        n_meas += 1
        tol = max(2.0, 3.0 * math.sqrt(e * e + ep * ep))
        ok = abs(t - tp) <= tol
        n_pass += ok
        rows.append({"n": n, "published": tp, "ours": t, "ours_err": e,
                     "diff": t - tp, "tol": tol, "pass": ok})
    return {"pass": n_pass >= 5, "n_measured": n_meas, "n_pass": n_pass, "rows": rows}


def gate_tess(tess_measured):
    """tess_measured: {n: (t, err)}; both published TESS peaks must be recovered."""
    rows, ok_all = [], True
    for n, tp, ep, sector in E.PUBLISHED_TESS_PEAKS:
        if n not in tess_measured:
            rows.append({"n": n, "sector": sector, "status": "not_measured"})
            ok_all = False
            continue
        t, e = tess_measured[n]
        ok = abs(t - tp) <= TESS_GATE_DAYS
        ok_all = ok_all and ok
        rows.append({"n": n, "sector": sector, "published": tp, "ours": t, "diff": t - tp, "pass": ok})
    return {"pass": ok_all, "rows": rows}


def forecast_test(measured):
    """T1: Payne+2022 forecasts n=25..29, and the Payne ephemeris for every new cycle."""
    rows = []
    for n, tf in sorted(E.PAYNE22_FORECASTS.items()):
        if n in measured:
            t, e = measured[n]
            rows.append({"n": n, "forecast": tf, "observed": t, "err": e,
                         "o_minus_c": t - tf, "on_time": abs(t - tf) <= ON_TIME_DAYS})
        else:
            rows.append({"n": n, "forecast": tf, "status": "not_measured"})
    drift = []
    for n in sorted(measured):
        if n >= E.FIRST_NEW_CYCLE:
            t, e = measured[n]
            drift.append({"n": n, "observed": t, "payne_ephemeris": E.payne_time(n),
                          "o_minus_c": t - E.payne_time(n)})
    meas = [r for r in rows if "observed" in r]
    return {"forecasts": rows, "n_measured": len(meas),
            "n_on_time": sum(1 for r in meas if r["on_time"]), "ephemeris_drift": drift}


def period_evolution(ns, ts, errs):
    """T2: steady (constant Pdot) and cubic fits with a common jitter."""
    n_ref = int(round(sum(ns) / len(ns)))
    jit = F.choose_jitter(ns, ts, errs, model="cubic", n_ref=n_ref)
    lin = F.fit_timing(ns, ts, errs, model="linear", jitter=jit, n_ref=n_ref)
    st = F.fit_timing(ns, ts, errs, model="steady", jitter=jit, n_ref=n_ref)
    cu = F.fit_timing(ns, ts, errs, model="cubic", jitter=jit, n_ref=n_ref)
    n_last = max(ns)
    pdot, perr = st["Pdot"], st["Pdot_err"]
    s = abs(pdot) / perr if perr > 0 else 0.0
    z = cu["cubic_z"]
    if s < SIG:
        outcome = "NO_DETECTED_DECAY"
    elif pdot > 0:
        outcome = "PERIOD_GROWING"
    elif z <= -SIG:
        outcome = "DECAY_ACCELERATING"
    elif z >= SIG:
        outcome = "DECAY_SLOWING"
    else:
        outcome = "DECAY_STEADY"
    p_last = F.period_at(st, n_last)
    dpdn = F.period_change_per_orbit(st, n_last)
    vs_payne = (pdot - E.PAYNE22["Pdot"]) / math.sqrt(perr ** 2 + E.PAYNE22["Pdot_err"] ** 2)
    return {
        "outcome": outcome, "n_cycles": len(ns), "n_first": min(ns), "n_last": n_last,
        "jitter_days": jit, "n_ref": n_ref,
        "Pdot": pdot, "Pdot_err": perr, "Pdot_sigma": s,
        "Pdot_vs_payne22_sigma": vs_payne,
        "period_at_n_ref": st["b"], "period_at_last": p_last,
        "period_change_hours_per_orbit_at_last": dpdn * 24.0,
        "decay_timescale_years_at_last": (p_last / abs(pdot) / 365.25) if pdot < 0 else None,
        "cubic_excess": cu["cubic_excess"], "cubic_excess_err": cu["cubic_excess_err"],
        "cubic_z": z,
        "chi2": {"linear": lin["chi2"], "steady": st["chi2"], "cubic": cu["chi2"]},
        "dof": {"linear": lin["dof"], "steady": st["dof"], "cubic": cu["dof"]},
        "steady_fit": {k: st[k] for k in ("a", "b", "c", "d", "n_ref")},
        "residuals_steady": [{"n": n, "t": t, "o_minus_c": t - F.predict(st, n)}
                             for n, t in zip(ns, ts)],
    }


def amplitude_trend(ns, amps, n_perm=10000, seed=20261121):
    """T3: Spearman rank correlation of flare amplitude with cycle number."""
    if len(ns) < 6:
        return {"outcome": "TOO_FEW", "n": len(ns)}
    rho, p = F.spearman_perm_p(ns, amps, n_perm, random.Random(seed))
    if p < TREND_P and rho < 0:
        outcome = "FADING"
    elif p < TREND_P and rho > 0:
        outcome = "BRIGHTENING"
    else:
        outcome = "NO_TREND"
    return {"outcome": outcome, "rho": rho, "p": p, "n": len(ns)}


def grade(g_results, tess_results=None):
    """g_results: list of measure_cycles dicts for ASAS-SN g. Returns the graded summary."""
    ok = [r for r in g_results if r.get("status") == "ok"]
    measured = {r["n"]: (r["t_peak"], r["err"]) for r in ok}
    out = {"gate_g1": gate_g1(measured)}
    if tess_results is not None:
        tm = {r["n"]: (r["t_peak"], r["err"]) for r in tess_results if r.get("status") == "ok"}
        out["gate_tess"] = gate_tess(tm)
    new = [n for n in measured if n >= E.FIRST_NEW_CYCLE]
    out["n_new_cycles"] = len(new)
    if not out["gate_g1"]["pass"]:
        out["outcome"] = "GATE_FAILED"
        return out
    if len(new) < MIN_NEW_CYCLES:
        out["outcome"] = "INSUFFICIENT_NEW_DATA"
        return out
    out["t1_forecasts"] = forecast_test(measured)
    ns = sorted(measured)
    out["t2_period"] = period_evolution(ns, [measured[n][0] for n in ns], [measured[n][1] for n in ns])
    out["outcome"] = out["t2_period"]["outcome"]
    # Secondary: full baseline = published peaks for n <= 24 plus ours for n >= 25.
    full = {n: (t, e) for n, t, e, _, _ in E.PUBLISHED_PEAKS}
    for n in new:
        full[n] = measured[n]
    fn = sorted(full)
    out["t2_full_baseline"] = period_evolution(fn, [full[n][0] for n in fn], [full[n][1] for n in fn])
    out["t3_amplitude"] = amplitude_trend([r["n"] for r in ok], [r["amplitude"] for r in ok])
    return out


def uqm_window(bins, base=(60900.0, 61230.0), window=(61250.0, 61330.0)):
    """Exploratory AT2023uqm check on ASAS-SN g (not graded): >= 2 bins > baseline + 4 sigma."""
    b = [x for x in bins if base[0] <= x[0] <= base[1]]
    w = [x for x in bins if window[0] <= x[0] <= window[1]]
    if len(b) < 5 or len(w) < 5:
        return {"status": "no_coverage", "n_base": len(b), "n_window": len(w)}
    fl = sorted(x[1] for x in b)
    med = fl[len(fl) // 2] if len(fl) % 2 else 0.5 * (fl[len(fl) // 2 - 1] + fl[len(fl) // 2])
    dev = sorted(abs(x[1] - med) for x in b)
    mad = dev[len(dev) // 2] if len(dev) % 2 else 0.5 * (dev[len(dev) // 2 - 1] + dev[len(dev) // 2])
    rs = 1.4826 * mad
    high = [x for x in w if x[1] - med > 4.0 * max(rs, x[2])]
    return {"status": "detected" if len(high) >= 2 else "not_detected",
            "baseline_flux": med, "robust_sigma": rs, "n_base": len(b), "n_window": len(w),
            "n_high": len(high), "max_window_flux": max(x[1] for x in w)}
