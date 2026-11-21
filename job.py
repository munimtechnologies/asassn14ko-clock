# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""asassn14ko-clock: is the orbit of the star in ASASSN-14ko still shrinking on schedule?

Runs the pre-registered analysis in PREREGISTRATION.md:
  1. ASAS-SN Sky Patrol v2 g-band light curve of ASASSN-14ko (public, pyasassn client)
  2. peak time of every flare, cycle by cycle (ko.peaks)
  3. gate G1 against published peaks, then tests T1-T3 (ko.grade)
  4. TESS FFI cross-check (gate G2), and an exploratory look at AT2023uqm's predicted
     August 2026 flare in the same archive.
Usage: python job.py --data-dir /data --out-dir /out
"""
import argparse
import csv
import json
import os
import sys
import time
import traceback

from ko import ephem as E
from ko import fit as F
from ko import grade as G
from ko import peaks as P

UQM = {"name": "AT2023uqm", "ra_deg": 347.66256, "dec_deg": 34.90301}


def log(*a):
    print(*a, flush=True)


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)


def write_csv(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols])


def run_asassn(data_dir, out_dir, results):
    from ko import asassn_io as A
    results["asassn_client"] = A.ensure_client()
    raw_csv = os.path.join(out_dir, "asassn_raw_14ko.csv")
    info = None
    for attempt in range(3):
        try:
            info = A.fetch(E.TARGET["ra_deg"], E.TARGET["dec_deg"], raw_csv, radius_arcsec=10.0, log=log)
            if info.get("status") == "not_in_catalog":
                info = A.fetch(E.TARGET["ra_deg"], E.TARGET["dec_deg"], raw_csv, radius_arcsec=30.0, log=log)
            break
        except Exception as exc:  # network hiccups: retry
            log("ASAS-SN fetch attempt %d failed: %r" % (attempt + 1, exc))
            info = {"status": "error", "error": repr(exc), "trace": traceback.format_exc()[-2000:]}
            time.sleep(30)
    results["asassn_fetch"] = info
    if not info or info.get("status") != "ok":
        return None
    pts, meta = A.load_points(raw_csv)
    results["asassn_meta"] = meta
    out = {}
    for band in ("g", "V"):
        cleaned4 = [p for p in pts[band] if p[2] > 0]
        norm, offsets = P.normalize_cameras(cleaned4)
        clean = P.clean_points(norm)
        bins = P.bin_points(clean, 1.0)
        out[band] = bins
        results.setdefault("camera_offsets", {})[band] = offsets
        results.setdefault("n_points", {})[band] = {"raw": len(pts[band]), "clean": len(clean), "bins": len(bins)}
        with open(os.path.join(out_dir, "lightcurve_%s_binned.csv" % band), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mjd", "flux_mjy", "err_mjy"])
            for t, f, e in bins:
                w.writerow(["%.4f" % t, "%.5f" % f, "%.5f" % e])
    return out


def run_tess(data_dir, out_dir, final_fit, results):
    from ko import tess_io as T
    secs = T.sectors(E.TARGET["ra_deg"], E.TARGET["dec_deg"])
    results["tess_sectors"] = secs
    tdir = os.path.join(data_dir, "tess")
    os.makedirs(tdir, exist_ok=True)
    rows = []
    for s in secs:
        try:
            path = T.cutout(E.TARGET["ra_deg"], E.TARGET["dec_deg"], s, tdir)
            lc = T.lightcurve(path, E.TARGET["ra_deg"], E.TARGET["dec_deg"])
        except Exception as exc:
            rows.append({"sector": s, "status": "error", "error": repr(exc)})
            continue
        if len(lc) < 50:
            rows.append({"sector": s, "status": "too_few_cadences", "n": len(lc)})
            continue
        bins = P.bin_points(P.clean_points(lc), 0.25)
        t0, t1 = bins[0][0], bins[-1][0]
        with open(os.path.join(out_dir, "tess_s%04d_binned.csv" % s), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mjd", "flux_e_per_s", "err"])
            for t, f, e in bins:
                w.writerow(["%.4f" % t, "%.3f" % f, "%.3f" % e])
        for n in range(0, 200):
            tp = F.predict(final_fit, n) if final_fit else E.payne_time(n)
            if tp < t0 + 3.0 or tp > t1 - 3.0:
                continue
            r = P.find_peak(bins, tp, 1000 + n, half_window=12.0)
            r.update({"n": n, "sector": s, "sector_start": t0, "sector_end": t1})
            rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--skip-tess", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.data_dir, exist_ok=True)
    results = {"project": "asassn14ko-clock", "version": "0.1.0",
               "run_started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "target": E.TARGET, "prereg": "PREREGISTRATION.md"}
    bins = None
    try:
        bins = run_asassn(args.data_dir, args.out_dir, results)
    except Exception as exc:
        results["asassn_error"] = repr(exc)
        results["asassn_trace"] = traceback.format_exc()[-3000:]
    final_fit = None
    g_cycles = []
    if bins and bins.get("g"):
        g_cycles = P.measure_cycles(bins["g"], E.payne_time, seed=20261121)
        cols = ["n", "status", "t_pred", "t_max", "t_peak", "err", "amplitude", "baseline",
                "snr", "n_window", "n_fit", "degree", "n_boot"]
        write_csv(os.path.join(args.out_dir, "peaks_g.csv"), g_cycles, cols)
        results["graded"] = G.grade(g_cycles)
        tp = results["graded"].get("t2_period")
        if tp:
            final_fit = dict(tp["steady_fit"])
        if bins.get("V"):
            v_cycles = P.measure_cycles(bins["V"], E.payne_time, seed=20261121)
            write_csv(os.path.join(args.out_dir, "peaks_V.csv"), v_cycles, cols)
            results["v_band_peaks_secondary"] = [
                {k: r.get(k) for k in ("n", "status", "t_peak", "err", "amplitude")} for r in v_cycles
                if r.get("status") == "ok"]
        results["last_g_mjd"] = bins["g"][-1][0] if bins["g"] else None
    else:
        results["outcome"] = "NO_DATA"
    if not args.skip_tess:
        try:
            tess_rows = run_tess(args.data_dir, args.out_dir, final_fit, results)
            results["tess_peaks"] = tess_rows
            tm = {r["n"]: (r["t_peak"], r["err"]) for r in tess_rows if r.get("status") == "ok"}
            results["gate_tess"] = G.gate_tess(tm)
            if final_fit:
                results["tess_vs_model"] = [
                    {"n": n, "sector": r["sector"], "t_tess": r["t_peak"],
                     "o_minus_c_steady": r["t_peak"] - F.predict(final_fit, n)}
                    for r in tess_rows if r.get("status") == "ok" for n in [r["n"]]]
        except Exception as exc:
            results["tess_error"] = repr(exc)
            results["tess_trace"] = traceback.format_exc()[-3000:]
    # Exploratory: AT2023uqm sixth-flare window in the same archive (not graded).
    try:
        from ko import asassn_io as A
        uqm_csv = os.path.join(args.out_dir, "asassn_raw_at2023uqm.csv")
        info = A.fetch(UQM["ra_deg"], UQM["dec_deg"], uqm_csv, radius_arcsec=5.0, log=log)
        results["at2023uqm_fetch"] = info
        if info.get("status") == "ok":
            pts, meta = A.load_points(uqm_csv)
            g = P.bin_points([(t, f, e) for t, f, e, _ in pts["g"] if e > 0], 1.0)
            results["at2023uqm_window"] = G.uqm_window(g)
    except Exception as exc:
        results["at2023uqm_error"] = repr(exc)
    if "graded" in results:
        results["outcome"] = results["graded"].get("outcome")
    results["run_finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    write_json(os.path.join(args.out_dir, "results.json"), results)
    log("OUTCOME:", results.get("outcome"))
    gr = results.get("graded", {})
    log("gate G1:", json.dumps(gr.get("gate_g1", {}), default=str)[:2000])
    if "t2_period" in gr:
        t2 = gr["t2_period"]
        log("Pdot = %.5f +- %.5f (%.1f sigma); cubic z = %.2f; P(last) = %.2f d; jitter %.2f d"
            % (t2["Pdot"], t2["Pdot_err"], t2["Pdot_sigma"], t2["cubic_z"], t2["period_at_last"],
               t2["jitter_days"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
