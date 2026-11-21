# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""ASAS-SN Sky Patrol v2 light curves through the public pyasassn client (network I/O only;
not imported by the tests)."""
import csv
import math
import subprocess
import sys


def ensure_client():
    try:
        from pyasassn.client import SkyPatrolClient  # noqa: F401
        return "already installed"
    except ImportError:
        pass
    for spec in ("skypatrol", "git+https://github.com/asas-sn/skypatrol.git"):
        r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", spec],
                           capture_output=True, text=True)
        try:
            from pyasassn.client import SkyPatrolClient  # noqa: F401
            return "installed " + spec
        except ImportError:
            last = (r.stdout or "")[-500:] + (r.stderr or "")[-1500:]
    raise RuntimeError("could not install pyasassn: " + last)


def _sep_arcsec(ra1, dec1, ra2, dec2):
    d1, d2 = math.radians(dec1), math.radians(dec2)
    dra = math.radians(ra2 - ra1)
    c = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(dra)
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 3600.0


def fetch(ra, dec, out_csv, radius_arcsec=10.0, log=print):
    """Download the light curve of the catalogued source nearest (ra, dec).
    Returns a dict describing what was found (or why not)."""
    from pyasassn.client import SkyPatrolClient
    client = SkyPatrolClient()
    rad = radius_arcsec / 3600.0
    cat = client.cone_search(ra_deg=ra, dec_deg=dec, radius=rad, catalog="master_list")
    info = {"radius_arcsec": radius_arcsec}
    if cat is None or len(cat) == 0:
        info["status"] = "not_in_catalog"
        return info
    rows = cat.to_dict("records")
    for r in rows:
        r["sep_arcsec"] = _sep_arcsec(ra, dec, float(r["ra_deg"]), float(r["dec_deg"]))
    rows.sort(key=lambda r: r["sep_arcsec"])
    best = rows[0]
    info["candidates"] = [{"asas_sn_id": str(r["asas_sn_id"]), "sep_arcsec": r["sep_arcsec"],
                           "catalog_sources": str(r.get("catalog_sources"))} for r in rows[:5]]
    info["asas_sn_id"] = str(best["asas_sn_id"])
    info["sep_arcsec"] = best["sep_arcsec"]
    lcs = client.cone_search(ra_deg=ra, dec_deg=dec, radius=rad, catalog="master_list",
                             download=True, threads=1)
    df = lcs.data
    df = df[df["asas_sn_id"].astype(str) == info["asas_sn_id"]]
    info["columns"] = [str(c) for c in df.columns]
    info["n_rows"] = int(len(df))
    df.to_csv(out_csv, index=False)
    info["status"] = "ok"
    log("ASAS-SN %s: %d rows, columns %s" % (info["asas_sn_id"], info["n_rows"], info["columns"]))
    return info


def load_points(csv_path, g_cut_mjd=58450.0):
    """Read the saved CSV into {'g': [(mjd, flux, err, camera)], 'V': [...]}.
    Band = phot_filter column when present, else by date: g from MJD 58450 (after the
    V-band cameras stopped), V before (see PREREGISTRATION.md)."""
    out = {"g": [], "V": []}
    meta = {"band_source": None, "n_total": 0, "n_bad_quality": 0, "last_mjd": None}
    with open(csv_path, newline="") as fh:
        rd = csv.DictReader(fh)
        cols = rd.fieldnames or []
        fcol = next((c for c in ("phot_filter", "filter", "band") if c in cols), None)
        meta["band_source"] = fcol or "date"
        for row in rd:
            meta["n_total"] += 1
            try:
                jd = float(row["jd"])
                flux = float(row["flux"])
                err = float(row["flux_err"])
            except (KeyError, ValueError):
                continue
            q = (row.get("quality") or "G").strip()
            if q and q != "G":
                meta["n_bad_quality"] += 1
                continue
            mjd = jd - 2400000.5 if jd > 2400000 else jd
            if fcol:
                band = (row.get(fcol) or "").strip()
                band = "g" if band.lower().startswith("g") else ("V" if band.upper().startswith("V") else band)
            else:
                band = "g" if mjd >= g_cut_mjd else "V"
            if band in out:
                out[band].append((mjd, flux, err, (row.get("camera") or "?").strip()))
            meta["last_mjd"] = mjd if meta["last_mjd"] is None else max(meta["last_mjd"], mjd)
    return out, meta
