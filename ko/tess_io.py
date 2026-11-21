# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""TESS full-frame-image cutouts from MAST TESSCut and simple aperture photometry
(network and numpy/astropy I/O only; not imported by the tests)."""
import io
import json
import os
import urllib.request
import zipfile

BASE = "https://mast.stsci.edu/tesscut/api/v0.1"
UA = {"User-Agent": "asassn14ko-clock/0.1 (open research script)"}


def _get(url, timeout=600):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def sectors(ra, dec):
    data = json.loads(_get("%s/sector?ra=%.6f&dec=%.6f&radius=1m" % (BASE, ra, dec)).decode())
    return sorted(int(r["sector"]) for r in data.get("results", []))


def cutout(ra, dec, sector, out_dir, size=15):
    path = os.path.join(out_dir, "tess_s%04d.fits" % sector)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    raw = _get("%s/astrocut?ra=%.6f&dec=%.6f&y=%d&x=%d&sector=%d" % (BASE, ra, dec, size, size, sector))
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names = [n for n in z.namelist() if n.endswith(".fits")]
        if not names:
            raise RuntimeError("no FITS in TESSCut response for sector %d" % sector)
        with open(path, "wb") as fh:
            fh.write(z.read(names[0]))
    return path


def lightcurve(path, ra, dec):
    """3x3-pixel aperture on the target minus a per-cadence background (median of pixels
    >= 5 px from the target). Returns [(mjd, flux, err)] for QUALITY == 0 cadences."""
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path) as h:
        d = h[1].data
        time = np.array(d["TIME"], dtype=float)
        flux = np.array(d["FLUX"], dtype=float)
        ferr = np.array(d["FLUX_ERR"], dtype=float)
        qual = np.array(d["QUALITY"], dtype=int)
        ny, nx = flux.shape[1], flux.shape[2]
        try:
            w = WCS(h[2].header)
            px, py = w.all_world2pix(ra, dec, 0)
            cx, cy = int(round(float(px))), int(round(float(py)))
        except Exception:
            cx, cy = nx // 2, ny // 2
    cx = min(max(cx, 1), nx - 2)
    cy = min(max(cy, 1), ny - 2)
    yy, xx = np.mgrid[0:ny, 0:nx]
    far = (np.hypot(xx - cx, yy - cy) >= 5.0)
    out = []
    for i in range(len(time)):
        if qual[i] != 0 or not np.isfinite(time[i]):
            continue
        img = flux[i]
        if not np.all(np.isfinite(img[cy - 1:cy + 2, cx - 1:cx + 2])):
            continue
        bk = np.nanmedian(img[far])
        f = float(np.sum(img[cy - 1:cy + 2, cx - 1:cx + 2]) - 9.0 * bk)
        e = float(np.sqrt(np.nansum(ferr[i][cy - 1:cy + 2, cx - 1:cx + 2] ** 2)))
        if np.isfinite(f) and np.isfinite(e) and e > 0:
            out.append((float(time[i]) + 56999.5, f, e))   # BTJD -> MJD
    return out
