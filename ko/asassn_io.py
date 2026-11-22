# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""ASAS-SN Sky Patrol v2 light curves through the public pyasassn client (network I/O only;
not imported by the tests).

0.1.1 (2026-11-22, deviation D1 in PREREGISTRATION.md): v0.1.0 took the light curve of the
nearest catalogued source from one cone download. On 2026-11-22 that download returned
"Pulled 1 of 2" and 0 rows for the nearest source (Sky Patrol is rebuilding its catalogues,
Hart, 2026-11-03). Now, for each candidate in order of separation, the code uses the cone
download's rows if there are any, otherwise asks for that source's light curve by ID (two
tries), and takes the first candidate within the radius that has rows. Which candidate was used,
and every count, is written to results.json.
"""
import csv
import math
import subprocess
import sys
import time


def ensure_client():
    try:
        from pyasassn.client import SkyPatrolClient  # noqa: F401
        return "already installed"
    except ImportError:
        pass
    last = ""
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


def _rows_for(df, asas_sn_id):
    if df is None or len(df) == 0 or "asas_sn_id" not in df.columns:
        return None
    return df[df["asas_sn_id"].astype(str) == str(asas_sn_id)]


def _by_id(client, asas_sn_id, log, tries=2):
    """Light curve of one source by ID; returns (dataframe or None, note)."""
    note = []
    for k in range(tries):
        try:
            lcs = client.query_list([int(asas_sn_id)], catalog="master_list",
                                    id_col="asas_sn_id", download=True, threads=1)
            df = getattr(lcs, "data", None)
            sub = _rows_for(df, asas_sn_id)
            n = 0 if sub is None else int(len(sub))
            note.append("try %d: %d rows" % (k + 1, n))
            if n > 0:
                return sub, "; ".join(note)
        except Exception as exc:  # network or API hiccup
            note.append("try %d: %r" % (k + 1, exc))
        time.sleep(5)
    log("ASAS-SN by-ID %s: %s" % (asas_sn_id, "; ".join(note)))
    return None, "; ".join(note)


def fetch(ra, dec, out_csv, radius_arcsec=10.0, log=print):
    """Download the light curve of the nearest catalogued source (ra, dec) that has rows.
    Returns a dict describing what was found (or why not)."""
    from pyasassn.client import SkyPatrolClient
    client = SkyPatrolClient()
    rad = radius_arcsec / 3600.0
    cat = client.cone_search(ra_deg=ra, dec_deg=dec, radius=rad, catalog="master_list")
    info = {"radius_arcsec": radius_arcsec, "io_version": "0.1.1"}
    if cat is None or len(cat) == 0:
        info["status"] = "not_in_catalog"
        return info
    rows = cat.to_dict("records")
    for r in rows:
        r["sep_arcsec"] = _sep_arcsec(ra, dec, float(r["ra_deg"]), float(r["dec_deg"]))
    rows.sort(key=lambda r: r["sep_arcsec"])
    rows = [r for r in rows if r["sep_arcsec"] <= radius_arcsec]
    info["candidates"] = [{"asas_sn_id": str(r["asas_sn_id"]), "sep_arcsec": r["sep_arcsec"],
                           "catalog_sources": str(r.get("catalog_sources"))} for r in rows[:5]]
    df_all = None
    try:
        lcs = client.cone_search(ra_deg=ra, dec_deg=dec, radius=rad, catalog="master_list",
                                 download=True, threads=1)
        df_all = getattr(lcs, "data", None)
    except Exception as exc:
        info["cone_download_error"] = repr(exc)
    tried = []
    chosen = None
    for rank, r in enumerate(rows[:5]):
        sid = str(r["asas_sn_id"])
        sub = _rows_for(df_all, sid)
        n_cone = 0 if sub is None else int(len(sub))
        entry = {"rank": rank, "asas_sn_id": sid, "sep_arcsec": r["sep_arcsec"], "n_cone": n_cone}
        if n_cone == 0:
            sub, note = _by_id(client, sid, log)
            entry["by_id"] = note
            entry["n_by_id"] = 0 if sub is None else int(len(sub))
        tried.append(entry)
        if sub is not None and len(sub) > 0:
            chosen = (rank, r, sub)
            break
    info["tried"] = tried
    if chosen is None:
        info["asas_sn_id"] = str(rows[0]["asas_sn_id"]) if rows else None
        info["n_rows"] = 0
        info["status"] = "no_rows"
        with open(out_csv, "w", newline="") as fh:
            fh.write("asas_sn_id,jd,flux,flux_err,mag,mag_err,limit,fwhm,image_id,camera,quality,phot_filter\n")
        log("ASAS-SN: no candidate within %.0f\" returned rows: %s" % (radius_arcsec, tried))
        return info
    rank, best, df = chosen
    info["asas_sn_id"] = str(best["asas_sn_id"])
    info["sep_arcsec"] = best["sep_arcsec"]
    info["chosen_rank"] = rank
    info["deviation_D1_fallback_used"] = rank > 0
    info["columns"] = [str(c) for c in df.columns]
    info["n_rows"] = int(len(df))
    df.to_csv(out_csv, index=False)
    info["status"] = "ok"
    log("ASAS-SN %s (rank %d, %.1f\"): %d rows" % (info["asas_sn_id"], rank, best["sep_arcsec"],
                                                   info["n_rows"]))
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
