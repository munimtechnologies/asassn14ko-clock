# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Flare peak times from a sparse light curve (pure Python, pre-registered settings).

Per cycle: look in a window around the predicted peak, take the brightest 3-bin running
mean as a provisional peak, fit a weighted polynomial (degree 4 / 3 / 2 for >=10 / >=7 / >=5
bins) to the bins from 12 d before to 15 d after it, and take the maximum of the polynomial
within +-8 d. The error is the bootstrap standard deviation (300 resamples of the bins).
"""
import math
import random

from . import fit as F

HALF_WINDOW = 25.0     # days around the predicted peak
PRE, POST = 12.0, 15.0  # fit range around the provisional peak
GRID_HALF = 8.0        # peak searched within +-8 d of the provisional peak
GRID_STEP = 0.05
EDGE = 0.1
MIN_WINDOW_BINS = 8
MIN_SNR = 5.0          # amplitude above baseline in units of the median bin error
N_BOOT = 300
ERR_FLOOR = 0.1


def clean_points(points, max_err_factor=3.0):
    """Drop non-finite points and those with errors > 3x the median error."""
    good = [(t, f, e) for t, f, e in points
            if all(math.isfinite(v) for v in (t, f, e)) and e > 0]
    if not good:
        return []
    errs = sorted(e for _, _, e in good)
    med = errs[len(errs) // 2]
    return sorted([p for p in good if p[2] <= max_err_factor * med])


def normalize_cameras(points, min_points=20):
    """points: (t, flux, err, camera). Subtract each camera's median-minus-global-median
    offset; cameras with fewer than min_points points are dropped. Returns ((t, f, e) list,
    {camera: offset})."""
    by = {}
    for p in points:
        by.setdefault(p[3], []).append(p)
    by = {cam: ps for cam, ps in by.items() if len(ps) >= min_points}
    if not by:
        return [], {}
    glob = median([p[1] for ps in by.values() for p in ps])
    offsets, out = {}, []
    for cam, ps in sorted(by.items()):
        off = median([p[1] for p in ps]) - glob
        offsets[cam] = off
        out.extend((t, f - off, e) for t, f, e, _ in ps)
    return sorted(out), offsets


def bin_points(points, width=1.0):
    """Inverse-variance weighted means in bins of floor(t / width)."""
    groups = {}
    for t, f, e in points:
        groups.setdefault(math.floor(t / width), []).append((t, f, e))
    out = []
    for key in sorted(groups):
        g = groups[key]
        ws = [1.0 / (e * e) for _, _, e in g]
        wsum = sum(ws)
        tm = sum(w * t for w, (t, _, _) in zip(ws, g)) / wsum
        fm = sum(w * f for w, (_, f, _) in zip(ws, g)) / wsum
        out.append((tm, fm, 1.0 / math.sqrt(wsum)))
    return out


def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def _poly_peak(fit_bins, degree, t_ref, lo, hi):
    xs = [(t - t_ref) / 10.0 for t, _, _ in fit_bins]
    try:
        coef, _, _ = F.wls_poly(xs, [f for _, f, _ in fit_bins], [e for _, _, e in fit_bins], degree)
    except ValueError:
        return None, None
    steps = int(round((hi - lo) / GRID_STEP))
    best_t, best_v = None, -float("inf")
    for i in range(steps + 1):
        t = lo + i * GRID_STEP
        v = F.polyval(coef, (t - t_ref) / 10.0)
        if v > best_v:
            best_t, best_v = t, v
    return best_t, best_v


def find_peak(bins, t_pred, seed, half_window=HALF_WINDOW, n_boot=N_BOOT,
              min_window_bins=MIN_WINDOW_BINS, min_snr=MIN_SNR):
    """Measure one flare peak near t_pred. bins: sorted list of (t, flux, err)."""
    res = {"t_pred": t_pred}
    win = [b for b in bins if abs(b[0] - t_pred) <= half_window]
    res["n_window"] = len(win)
    if len(win) < min_window_bins:
        res["status"] = "sparse"
        return res
    best_i, best_v = 0, -float("inf")
    for i in range(len(win)):
        lo, hi = max(0, i - 1), min(len(win), i + 2)
        v = sum(b[1] for b in win[lo:hi]) / (hi - lo)
        if v > best_v:
            best_i, best_v = i, v
    t_max = win[best_i][0]
    res["t_max"] = t_max
    fit_bins = [b for b in bins if t_max - PRE <= b[0] <= t_max + POST]
    n_pre = sum(1 for b in fit_bins if b[0] < t_max - 0.5)
    n_post = sum(1 for b in fit_bins if b[0] > t_max + 0.5)
    res["n_fit"] = len(fit_bins)
    if len(fit_bins) < 5 or n_pre < 2 or n_post < 2:
        res["status"] = "unbracketed"
        return res
    deg = 4 if len(fit_bins) >= 10 else (3 if len(fit_bins) >= 7 else 2)
    res["degree"] = deg
    lo = max(t_max - GRID_HALF, fit_bins[0][0])
    hi = min(t_max + GRID_HALF, fit_bins[-1][0])
    peak, val = _poly_peak(fit_bins, deg, t_max, lo, hi)
    if peak is None:
        res["status"] = "fit_failed"
        return res
    res["t_peak"] = peak
    res["peak_flux"] = val
    if peak - lo < EDGE or hi - peak < EDGE:
        res["status"] = "edge"
        return res
    base = [b[1] for b in bins if peak - 45.0 <= b[0] <= peak - 15.0]
    if len(base) >= 3:
        baseline = median(base)
    else:
        wf = sorted(b[1] for b in win)
        baseline = wf[int(0.1 * (len(wf) - 1))]
    noise = median([b[2] for b in fit_bins])
    amp = val - baseline
    res["baseline"] = baseline
    res["amplitude"] = amp
    res["snr"] = amp / noise if noise > 0 else 0.0
    if amp < min_snr * noise:
        res["status"] = "no_flare"
        return res
    rng = random.Random(seed)
    boots = []
    for _ in range(n_boot):
        sample = sorted(fit_bins[rng.randrange(len(fit_bins))] for _ in fit_bins)
        if len(set(b[0] for b in sample)) < deg + 2:
            continue
        p, _ = _poly_peak(sample, deg, t_max, lo, hi)
        if p is not None:
            boots.append(p)
    if len(boots) < 20:
        res["status"] = "unstable"
        return res
    mb = sum(boots) / len(boots)
    sd = math.sqrt(sum((b - mb) ** 2 for b in boots) / (len(boots) - 1))
    res["err"] = max(ERR_FLOOR, sd)
    res["n_boot"] = len(boots)
    res["status"] = "ok"
    return res


def measure_cycles(bins, first_guess, seed=20261121, n_start=0, n_max=200, **kw):
    """Walk cycle by cycle. The window is centred on first_guess(n) until 5 peaks are
    measured, then on a constant-Pdot fit to the peaks measured so far."""
    if not bins:
        return []
    last_t = bins[-1][0]
    half = kw.get("half_window", HALF_WINDOW)
    found, out = [], []
    for n in range(n_start, n_max + 1):
        if len(found) >= 5:
            f = F.fit_timing([x[0] for x in found], [x[1] for x in found],
                             [x[2] for x in found], model="steady", jitter=1.0)
            t_pred = F.predict(f, n)
        else:
            t_pred = first_guess(n)
        if t_pred + half > last_t:
            break
        r = find_peak(bins, t_pred, seed + n, **kw)
        r["n"] = n
        out.append(r)
        if r["status"] == "ok":
            found.append((n, r["t_peak"], r["err"]))
    return out
