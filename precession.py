#!/usr/bin/env python3
# Copyright 2026 Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
"""Search ASASSN-14ko's flare arrival times for a periodic wobble on top of steady decay.

Linial & Quataert (2024, MNRAS 527, 4317; arXiv:2309.15849) explain the shrinking flare period
of ASASSN-14ko by drag as the star crosses the black hole's gas disc twice per orbit. They
expect the drag, and so the period derivative, to be modulated as the orbit precesses, which
leaves an oscillation in the flare times (their Fig. 1; period T_prec/2), and they encourage
"further analysis of the flare timing data, and its deviations from constant Pdot".

Model (m = n - 17, u = m / 17):
    t_n = a + b m + c m^2 + d m^3 + A cos(2 pi m / N) + B sin(2 pi m / N)
Detection: F statistic of the sinusoid over the cubic-only fit, maximised over a frequency
grid; global false-alarm probability (FAP) from white-noise simulations (the statistic does
not depend on the noise level or the polynomial). Upper limit K95(N): the true amplitude for
which 95% of simulated data sets (random phase, noise = cubic-fit rms) give a fitted amplitude
at least as large as the observed one. Converted to a fractional modulation of Pdot:
    Pdot(m) = Pdot0 (1 + f cos(2 pi m / N))  =>  K = f |Pdot0| P (N / 2 pi)^2.
Standard library only. Rules: AMENDMENT-2026-11-27-precession.md.
"""
import argparse
import csv
import json
import math
import os
import random
import sys
import time

VERSION = "0.1.0"
N_REF = 17
SCALE = 17.0
N_MIN, N_MAX = 4.0, 34.0
OVERSAMPLE = 10
UL_GRID = [6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 25.0, 30.0]
NSIM_FAP = 2000
NSIM_UL = 500
SEED = 20261127
FAP_DETECT = 0.01
FAP_HINT = 0.05
OSC_PER_RP = 6.5  # oscillation period in orbits per (r_p / 2 r_t): T_prec/(2 P0), LQ24 eq. 14 / sec. 4.3


def read_timings(path):
    ns, ts = [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            ns.append(int(row["n"]))
            ts.append(float(row["t_mjd"]))
    return ns, ts


def solve(a, b):
    n = len(b)
    m = [list(map(float, a[i])) + [float(b[i])] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-12:
            raise ValueError("singular system")
        m[col], m[piv] = m[piv], m[col]
        p = m[col][col]
        for r in range(col + 1, n):
            f = m[r][col] / p
            if f:
                row_r, row_c = m[r], m[col]
                for k in range(col, n + 1):
                    row_r[k] -= f * row_c[k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = m[i][n] - sum(m[i][k] * x[k] for k in range(i + 1, n))
        x[i] = s / m[i][i]
    return x


def lstsq(X, y):
    p = len(X[0])
    ata = [[sum(r[i] * r[j] for r in X) for j in range(p)] for i in range(p)]
    aty = [sum(r[i] * yy for r, yy in zip(X, y)) for i in range(p)]
    beta = solve(ata, aty)
    rss = 0.0
    for r, yy in zip(X, y):
        e = yy - sum(bb * v for bb, v in zip(beta, r))
        rss += e * e
    return beta, rss


def base_rows(ns):
    rows = []
    for n in ns:
        u = (n - N_REF) / SCALE
        rows.append([1.0, u, u * u, u * u * u])
    return rows


def sin_rows(ns, period):
    rows = []
    for r, n in zip(base_rows(ns), ns):
        ph = 2.0 * math.pi * (n - N_REF) / period
        rows.append(r + [math.cos(ph), math.sin(ph)])
    return rows


def fit_base(ns, ts):
    y0 = ts[0]
    y = [t - y0 for t in ts]
    beta, rss = lstsq(base_rows(ns), y)
    b = beta[1] / SCALE
    c = beta[2] / SCALE ** 2
    d = beta[3] / SCALE ** 3
    return {"a": beta[0] + y0, "P_days": b, "c": c, "d": d, "Pdot": 2.0 * c / b, "rss": rss,
            "rms": math.sqrt(rss / (len(ns) - 4))}


def freq_grid(ns):
    span = max(ns) - min(ns)
    df = 1.0 / (span * OVERSAMPLE)
    f, fmax = 1.0 / N_MAX, 1.0 / N_MIN
    out = []
    while f <= fmax + 1e-12:
        out.append(1.0 / f)
        f += df
    return out


class Scanner:
    """Precomputes design matrices so repeated fits on the same epochs are cheap."""

    def __init__(self, ns, periods):
        self.ns = list(ns)
        self.periods = list(periods)
        self.X0 = base_rows(ns)
        self.Xs = [sin_rows(ns, p) for p in periods]
        self.dof1 = len(ns) - 6

    def fit_one(self, y, i):
        beta, rss1 = lstsq(self.Xs[i], y)
        return math.hypot(beta[4], beta[5]), math.atan2(-beta[5], beta[4]), rss1

    def scan(self, y):
        _, rss0 = lstsq(self.X0, y)
        out = []
        for i, p in enumerate(self.periods):
            K, phase, rss1 = self.fit_one(y, i)
            F = ((rss0 - rss1) / 2.0) / (rss1 / self.dof1) if rss1 > 0 else float("inf")
            out.append({"N": p, "F": F, "K": K, "phase": phase})
        return out


def fap_of(scanner, f_obs, nsim, seed):
    rng = random.Random(seed)
    hits = 0
    for _ in range(nsim):
        y = [rng.gauss(0.0, 1.0) for _ in scanner.ns]
        if max(r["F"] for r in scanner.scan(y)) >= f_obs:
            hits += 1
    return (hits + 1.0) / (nsim + 1.0)


def upper_limit(ns, period, k_obs, sigma, nsim, seed, cl=0.95):
    sc = Scanner(ns, [period])
    rng = random.Random(seed)
    noise = [[rng.gauss(0.0, sigma) for _ in ns] for _ in range(nsim)]
    phases = [rng.uniform(0.0, 2.0 * math.pi) for _ in range(nsim)]
    ms = [n - N_REF for n in ns]

    def frac(k):
        hits = 0
        for e, ph in zip(noise, phases):
            y = [k * math.cos(2.0 * math.pi * m / period + ph) + ee for m, ee in zip(ms, e)]
            kf, _, _ = sc.fit_one(y, 0)
            if kf >= k_obs:
                hits += 1
        return hits / float(nsim)

    lo, hi = 0.0, max(4.0 * k_obs, 10.0 * sigma)
    while frac(hi) < cl and hi < 1e4:
        hi *= 2.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if frac(mid) >= cl:
            hi = mid
        else:
            lo = mid
    return hi


def analyse(ns, ts, nsim_fap=NSIM_FAP, nsim_ul=NSIM_UL, ul_grid=None, seed=SEED, jackknife=True):
    ul_grid = UL_GRID if ul_grid is None else ul_grid
    base = fit_base(ns, ts)
    y = [t - ts[0] for t in ts]
    periods = freq_grid(ns)
    sc = Scanner(ns, periods)
    rows = sc.scan(y)
    best = max(rows, key=lambda r: r["F"])
    fap = fap_of(sc, best["F"], nsim_fap, seed)
    if fap < FAP_DETECT:
        outcome = "WOBBLE_DETECTED"
    elif fap < FAP_HINT:
        outcome = "WOBBLE_HINT"
    else:
        outcome = "NO_WOBBLE"
    pdp = abs(base["Pdot"]) * base["P_days"]  # period change per orbit, days
    limits = []
    for j, N in enumerate(ul_grid):
        k_obs, _, _ = Scanner(ns, [N]).fit_one(y, 0)
        k95 = upper_limit(ns, N, k_obs, base["rms"], nsim_ul, seed + 1 + j)
        f95 = k95 * (2.0 * math.pi / N) ** 2 / pdp
        limits.append({"N": N, "K_obs": k_obs, "K95": k95, "f95": f95, "rp_over_2rt": N / OSC_PER_RP})
    n_star = None
    for lim in sorted(limits, key=lambda r: -r["N"]):
        if lim["f95"] < 1.0:
            n_star = lim["N"]
        else:
            break
    out = {"version": VERSION, "n_flares": len(ns), "base": base, "n_periods": len(periods),
           "best": best, "fap": fap, "nsim_fap": nsim_fap, "outcome": outcome,
           "limits": limits, "N_star_f_below_1": n_star,
           "rp_over_2rt_at_N_star": (n_star / OSC_PER_RP) if n_star else None,
           "scan": rows}
    if jackknife and 23 in ns:
        keep = [i for i, n in enumerate(ns) if n != 23]
        ns2 = [ns[i] for i in keep]
        y2 = [y[i] for i in keep]
        sc2 = Scanner(ns2, periods)
        b2 = max(sc2.scan(y2), key=lambda r: r["F"])
        out["jackknife_drop_n23"] = {"best": b2, "fap": fap_of(sc2, b2["F"], max(200, nsim_fap // 2), seed + 99)}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="timings_29.csv")
    ap.add_argument("--results-json", default="", help="optional cross-check against results/v0.1.0/results.json")
    ap.add_argument("--out-dir", default="out")
    a = ap.parse_args(argv)
    t0 = time.time()
    os.makedirs(a.out_dir, exist_ok=True)
    ns, ts = read_timings(a.data)
    check = "not run"
    if a.results_json and os.path.exists(a.results_json):
        with open(a.results_json) as fh:
            rj = json.load(fh)
        ref = {r["n"]: r["t"] for r in rj["graded"]["t2_full_baseline"]["residuals_steady"]}
        bad = [n for n, t in zip(ns, ts) if n not in ref or abs(ref[n] - t) > 1e-6]
        check = "ok" if not bad and len(ref) == len(ns) else "MISMATCH %s" % bad
    print("[%6.1fs] %d flares, cross-check vs results.json: %s" % (time.time() - t0, len(ns), check), flush=True)
    res = analyse(ns, ts)
    res["crosscheck_results_json"] = check
    with open(os.path.join(a.out_dir, "precession_results.json"), "w") as fh:
        json.dump({k: v for k, v in res.items() if k != "scan"}, fh, indent=1)
    with open(os.path.join(a.out_dir, "scan.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N_orbits", "F", "K_days", "phase_rad"])
        for r in res["scan"]:
            w.writerow(["%.4f" % r["N"], "%.5f" % r["F"], "%.4f" % r["K"], "%.4f" % r["phase"]])
    with open(os.path.join(a.out_dir, "limits.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N_orbits", "K_obs_days", "K95_days", "f95", "rp_over_2rt"])
        for r in res["limits"]:
            w.writerow(["%.1f" % r["N"], "%.3f" % r["K_obs"], "%.3f" % r["K95"], "%.3f" % r["f95"], "%.2f" % r["rp_over_2rt"]])
    b = res["base"]
    lines = ["# ASASSN-14ko precession-wobble search v%s" % VERSION, "",
             "- flares: %d; cubic-only fit: P = %.3f d at n=17, Pdot = %.5f, rms = %.2f d" % (len(ns), b["P_days"], b["Pdot"], b["rms"]),
             "- best period %.2f orbits: F = %.2f, K = %.2f d; global FAP = %.4f (%d sims)" % (res["best"]["N"], res["best"]["F"], res["best"]["K"], res["fap"], res["nsim_fap"]),
             "- OUTCOME %s" % res["outcome"],
             "- N* (f95 < 1 for all grid N >= N*): %s" % res["N_star_f_below_1"]]
    if "jackknife_drop_n23" in res:
        jk = res["jackknife_drop_n23"]
        lines.append("- without flare 23: best %.2f orbits, F = %.2f, FAP = %.4f" % (jk["best"]["N"], jk["best"]["F"], jk["fap"]))
    lines += ["", "| N (orbits) | K_obs (d) | K95 (d) | f95 | r_p/2r_t |", "|---|---|---|---|---|"]
    for r in res["limits"]:
        lines.append("| %.0f | %.2f | %.2f | %.2f | %.2f |" % (r["N"], r["K_obs"], r["K95"], r["f95"], r["rp_over_2rt"]))
    with open(os.path.join(a.out_dir, "summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    print("[%6.1fs] done" % (time.time() - t0), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
