# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Dependency-free weighted polynomial least squares and flare-timing models.

Timing models, with m = n - n_ref:
  linear : t = a + b m
  steady : t = a + b m + c m^2 + d m^3 with d tied to 2 c^2 / (3 b)  (constant dP/dt,
           Payne et al. 2022 eq. 1 expanded about n_ref; P = b, Pdot = 2c / b)
  cubic  : t = a + b m + c m^2 + d m^3 with d free (tests whether dP/dt itself changes)
"""
import math


def solve_inverse(matrix):
    """Inverse of a small square matrix (Gauss-Jordan, partial pivoting)."""
    n = len(matrix)
    scale = max(abs(v) for row in matrix for v in row) or 1.0
    a = [[float(v) for v in row] + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(matrix)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[piv][col]) <= 1e-13 * scale:
            raise ValueError("singular matrix")
        a[col], a[piv] = a[piv], a[col]
        p = a[col][col]
        a[col] = [v / p for v in a[col]]
        for r in range(n):
            if r != col:
                f = a[r][col]
                if f != 0.0:
                    a[r] = [vr - f * vc for vr, vc in zip(a[r], a[col])]
    return [row[n:] for row in a]


def polyval(coef, x):
    """sum_j coef[j] * x**j (Horner)."""
    y = 0.0
    for c in reversed(coef):
        y = y * x + c
    return y


def wls_poly(xs, ys, sigmas, degree):
    """Weighted least-squares polynomial. Returns (coef, covariance, chi2)."""
    k = degree + 1
    if len(xs) < k:
        raise ValueError("not enough points for the polynomial degree")
    ata = [[0.0] * k for _ in range(k)]
    aty = [0.0] * k
    for x, y, s in zip(xs, ys, sigmas):
        w = 1.0 / (s * s)
        pw = [x ** j for j in range(k)]
        for i in range(k):
            aty[i] += w * pw[i] * y
            for j in range(k):
                ata[i][j] += w * pw[i] * pw[j]
    cov = solve_inverse(ata)
    coef = [sum(cov[i][j] * aty[j] for j in range(k)) for i in range(k)]
    chi2 = sum(((y - polyval(coef, x)) / s) ** 2 for x, y, s in zip(xs, ys, sigmas))
    return coef, cov, chi2


_U = 10.0  # m is scaled to u = m / 10 for conditioning


def fit_timing(ns, ts, errs, model="steady", jitter=0.0, n_ref=None):
    """Fit peak times t(n). Returns a dict with a, b (=P at n_ref), c, d, covariance of
    (a, b, c[, d]) in m-units, chi2, dof, Pdot and its error (steady/cubic)."""
    if n_ref is None:
        n_ref = int(round(sum(ns) / len(ns)))
    sig = [math.sqrt(e * e + jitter * jitter) for e in errs]
    us = [(n - n_ref) / _U for n in ns]
    ms = [n - n_ref for n in ns]
    if model == "linear":
        deg = 1
    elif model == "steady":
        deg = 2
    elif model == "cubic":
        deg = 3
    else:
        raise ValueError(model)
    d_fix = 0.0
    for _ in range(100):
        ys = [t - d_fix * m ** 3 for t, m in zip(ts, ms)] if model == "steady" else list(ts)
        coef, cov_u, chi2 = wls_poly(us, ys, sig, deg)
        if model != "steady":
            break
        b = coef[1] / _U
        c = coef[2] / _U ** 2
        new = 2.0 * c * c / (3.0 * b)
        if abs(new - d_fix) < 1e-14:
            d_fix = new
            break
        d_fix = new
    # Convert coefficients and covariance from u-units to m-units.
    k = deg + 1
    fac = [1.0 / _U ** j for j in range(k)]
    coef_m = [coef[j] * fac[j] for j in range(k)]
    cov_m = [[cov_u[i][j] * fac[i] * fac[j] for j in range(k)] for i in range(k)]
    out = {
        "model": model, "n_ref": n_ref, "jitter": jitter, "npts": len(ns),
        "a": coef_m[0], "b": coef_m[1],
        "c": coef_m[2] if k > 2 else 0.0,
        "d": coef_m[3] if k > 3 else (d_fix if model == "steady" else 0.0),
        "cov": cov_m, "chi2": chi2, "dof": len(ns) - k,
    }
    if model in ("steady", "cubic"):
        b, c = out["b"], out["c"]
        vb, vc, vbc = cov_m[1][1], cov_m[2][2], cov_m[1][2]
        pdot = 2.0 * c / b
        gb, gc = -2.0 * c / (b * b), 2.0 / b
        var = gb * gb * vb + gc * gc * vc + 2.0 * gb * gc * vbc
        out["Pdot"] = pdot
        out["Pdot_err"] = math.sqrt(max(var, 0.0))
    if model == "cubic":
        b, c, d = out["b"], out["c"], out["d"]
        g = [0.0, 2.0 * c * c / (3.0 * b * b), -4.0 * c / (3.0 * b), 1.0]
        var = sum(g[i] * g[j] * cov_m[i][j] for i in range(4) for j in range(4))
        excess = d - 2.0 * c * c / (3.0 * b)
        out["cubic_excess"] = excess
        out["cubic_excess_err"] = math.sqrt(max(var, 0.0))
        out["cubic_z"] = excess / out["cubic_excess_err"] if out["cubic_excess_err"] > 0 else 0.0
    return out


def predict(fit, n):
    m = n - fit["n_ref"]
    return fit["a"] + fit["b"] * m + fit["c"] * m * m + fit["d"] * m ** 3


def period_at(fit, n):
    """Instantaneous period dt/dn at cycle n (days)."""
    m = n - fit["n_ref"]
    return fit["b"] + 2.0 * fit["c"] * m + 3.0 * fit["d"] * m * m


def period_change_per_orbit(fit, n):
    """d(period)/dn at cycle n (days per orbit)."""
    m = n - fit["n_ref"]
    return 2.0 * fit["c"] + 6.0 * fit["d"] * m


def choose_jitter(ns, ts, errs, model="cubic", n_ref=None, hi=50.0):
    """Smallest extra error (days, added in quadrature) giving reduced chi2 <= 1."""
    f0 = fit_timing(ns, ts, errs, model=model, jitter=0.0, n_ref=n_ref)
    if f0["dof"] <= 0 or f0["chi2"] / f0["dof"] <= 1.0:
        return 0.0
    lo = 0.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = fit_timing(ns, ts, errs, model=model, jitter=mid, n_ref=n_ref)
        if f["chi2"] / f["dof"] > 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ranks(values):
    """Average ranks (1-based), ties share the mean rank."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    r = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return 0.0
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def spearman_perm_p(x, y, n_perm, rng):
    """Two-sided permutation p-value of Spearman's rho."""
    rho = spearman(x, y)
    ry = ranks(y)
    rx = ranks(x)
    yy = list(ry)
    hits = 0
    for _ in range(n_perm):
        rng.shuffle(yy)
        if abs(pearson(rx, yy)) >= abs(rho) - 1e-12:
            hits += 1
    return rho, (hits + 1.0) / (n_perm + 1.0)
