# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests for the pure-Python parts of asassn14ko-clock (no network, no numpy)."""
import math
import unittest

from ko import ephem as E
from ko import fit as F
from ko import grade as G
from ko import peaks as P


def gauss_lc(centers, t0, t1, amp=0.5, sigma=6.0, base=1.0, err=0.01, step=1.0):
    bins = []
    t = t0
    while t <= t1:
        f = base + sum(amp * math.exp(-((t - c) / sigma) ** 2) for c in centers)
        bins.append((t, f, err))
        t += step
    return bins


class TestEphem(unittest.TestCase):
    def test_payne_forecasts_reproduced(self):
        for n, tf in E.PAYNE22_FORECASTS.items():
            self.assertLess(abs(E.payne_time(n) - tf), 1.0, n)

    def test_published_mjd_conversion(self):
        pub = {n: t for n, t, _, _, _ in E.PUBLISHED_PEAKS}
        self.assertAlmostEqual(pub[0], 56860.7, places=3)
        self.assertAlmostEqual(pub[15], 58546.6, places=3)
        self.assertAlmostEqual(E.PUBLISHED_TESS_PEAKS[1][1], 59205.63, places=3)


class TestFit(unittest.TestCase):
    def test_wls_exact_quadratic(self):
        xs = [-2.0, -1.0, 0.0, 1.0, 2.0]
        ys = [1 + 2 * x - 0.5 * x * x for x in xs]
        coef, cov, chi2 = F.wls_poly(xs, ys, [1.0] * 5, 2)
        for got, want in zip(coef, [1.0, 2.0, -0.5]):
            self.assertAlmostEqual(got, want, places=9)
        self.assertLess(chi2, 1e-15)

    def test_wls_constant_variance(self):
        coef, cov, _ = F.wls_poly([0, 1, 2, 3], [5, 5, 5, 5], [2.0] * 4, 0)
        self.assertAlmostEqual(coef[0], 5.0)
        self.assertAlmostEqual(cov[0][0], 1.0)

    def test_steady_recovers_constant_pdot(self):
        k = -0.0026
        ns = list(range(41))
        ts = [56851.5 + 115.2 * (math.exp(k * n) - 1.0) / k for n in ns]
        f = F.fit_timing(ns, ts, [1.0] * len(ns), model="steady", n_ref=20)
        self.assertLess(abs(f["Pdot"] - k), 1e-5)
        self.assertLess(abs(F.period_at(f, 20) - 115.2 * math.exp(k * 20)), 0.05)
        c = F.fit_timing(ns, ts, [1.0] * len(ns), model="cubic", n_ref=20)
        self.assertLess(abs(c["cubic_z"]), 1.0)

    def test_ranks_and_spearman(self):
        self.assertEqual(F.ranks([10, 20, 20, 30]), [1.0, 2.5, 2.5, 4.0])
        self.assertAlmostEqual(F.spearman([1, 2, 3], [3, 2, 1]), -1.0)


class TestPeriodEvolution(unittest.TestCase):
    def _cubic(self, d):
        ns = list(range(41))
        ts = [1000 + 115 * n - 0.15 * n * n + d * n ** 3 for n in ns]
        return G.period_evolution(ns, ts, [0.5] * len(ns))

    def test_accelerating(self):
        self.assertEqual(self._cubic(-0.002)["outcome"], "DECAY_ACCELERATING")

    def test_slowing(self):
        self.assertEqual(self._cubic(0.002)["outcome"], "DECAY_SLOWING")

    def test_no_decay(self):
        ns = list(range(30))
        r = G.period_evolution(ns, [100.0 + 110.0 * n for n in ns], [1.0] * len(ns))
        self.assertEqual(r["outcome"], "NO_DETECTED_DECAY")


class TestPeaks(unittest.TestCase):
    def test_bin_points(self):
        b = P.bin_points([(10.2, 1.0, 1.0), (10.8, 3.0, 1.0), (11.5, 5.0, 2.0)], 1.0)
        self.assertEqual(len(b), 2)
        self.assertAlmostEqual(b[0][0], 10.5)
        self.assertAlmostEqual(b[0][1], 2.0)
        self.assertAlmostEqual(b[0][2], 1 / math.sqrt(2))
        self.assertEqual(b[1], (11.5, 5.0, 2.0))

    def test_normalize_cameras(self):
        pts = [(float(i), 1.0, 0.1, "A") for i in range(25)]
        pts += [(float(i) + 0.5, 1.2, 0.1, "B") for i in range(25)]
        pts += [(100.0 + i, 9.0, 0.1, "C") for i in range(3)]
        out, off = P.normalize_cameras(pts)
        self.assertEqual(len(out), 50)
        self.assertNotIn("C", off)
        for _, f, _ in out:
            self.assertAlmostEqual(f, 1.1)

    def test_find_peak_gaussian(self):
        bins = gauss_lc([1050.3], 1000, 1100)
        r = P.find_peak(bins, 1048.0, seed=1, n_boot=60)
        self.assertEqual(r["status"], "ok")
        self.assertLess(abs(r["t_peak"] - 1050.3), 1.0)
        self.assertGreater(r["snr"], 20)

    def test_find_peak_sparse_and_flat(self):
        bins = gauss_lc([], 1000, 1100)
        self.assertEqual(P.find_peak(bins[:30:6], 1010.0, seed=1)["status"], "sparse")
        self.assertNotEqual(P.find_peak(bins, 1050.0, seed=1, n_boot=30)["status"], "ok")

    def test_measure_cycles_tracks_ephemeris(self):
        true = {n: E.payne_time(n) for n in range(20, 31)}
        bins = gauss_lc(list(true.values()), true[20] - 60, true[30] + 40)
        res = P.measure_cycles(bins, E.payne_time, n_boot=40)
        ok = {r["n"]: r for r in res if r["status"] == "ok"}
        for n, t in true.items():
            self.assertIn(n, ok)
            self.assertLess(abs(ok[n]["t_peak"] - t), 1.0, n)


class TestGrade(unittest.TestCase):
    def test_gate_g1(self):
        pub = {n: t for n, t, _, _, _ in E.PUBLISHED_PEAKS}
        m = {n: (pub[n], 0.5) for n in E.GATE_CYCLES}
        self.assertTrue(G.gate_g1(m)["pass"])
        m[15] = (pub[15] + 10.0, 0.5)
        m[17] = (pub[17] - 10.0, 0.5)
        g = G.gate_g1(m)
        self.assertFalse(g["pass"])
        self.assertEqual(g["n_pass"], 4)

    def test_grade_end_to_end_steady(self):
        rows = [{"n": n, "status": "ok", "t_peak": E.payne_time(n), "err": 0.5,
                 "amplitude": 1.0 - 0.01 * n} for n in range(11, 41)]
        out = G.grade(rows)
        self.assertTrue(out["gate_g1"]["pass"])
        self.assertEqual(out["outcome"], "DECAY_STEADY")
        self.assertLess(abs(out["t2_period"]["Pdot"] + 0.0026), 3e-5)
        self.assertEqual(out["t1_forecasts"]["n_on_time"], 5)
        self.assertEqual(out["t3_amplitude"]["outcome"], "FADING")
        self.assertIn("t2_full_baseline", out)

    def test_grade_needs_new_cycles(self):
        rows = [{"n": n, "status": "ok", "t_peak": E.payne_time(n), "err": 0.5,
                 "amplitude": 1.0} for n in range(11, 27)]
        self.assertEqual(G.grade(rows)["outcome"], "INSUFFICIENT_NEW_DATA")

    def test_uqm_window(self):
        base = [(60900.0 + 10 * i, 0.95 if i % 2 else 1.05, 0.1) for i in range(10)]
        win = [(61250.0 + 10 * i, 2.0 if i < 2 else 1.0, 0.1) for i in range(6)]
        self.assertEqual(G.uqm_window(base + win)["status"], "detected")
        self.assertEqual(G.uqm_window(base)["status"], "no_coverage")


if __name__ == "__main__":
    unittest.main()
