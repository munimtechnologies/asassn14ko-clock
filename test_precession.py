# Copyright 2026 Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
"""Unit tests for precession.py (standard library only)."""
import math
import os
import random
import unittest

import precession as P

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "timings_29.csv")


class TestPrecession(unittest.TestCase):
    def test_solve(self):
        x = P.solve([[2.0, 1.0], [1.0, 3.0]], [3.0, 5.0])
        self.assertAlmostEqual(x[0], 0.8, places=10)
        self.assertAlmostEqual(x[1], 1.4, places=10)

    def test_cubic_fit_exact(self):
        ns = list(range(0, 35))
        ts = [58000.0 + 110.0 * (n - 17) - 0.1 * (n - 17) ** 2 + 1e-4 * (n - 17) ** 3 for n in ns]
        b = P.fit_base(ns, ts)
        self.assertAlmostEqual(b["P_days"], 110.0, places=6)
        self.assertAlmostEqual(b["c"], -0.1, places=8)
        self.assertAlmostEqual(b["d"], 1e-4, places=9)
        self.assertAlmostEqual(b["Pdot"], -0.2 / 110.0, places=9)
        self.assertLess(b["rss"], 1e-12)

    def test_freq_grid(self):
        g = P.freq_grid(list(range(0, 35)))
        self.assertAlmostEqual(g[0], 34.0, places=6)
        self.assertGreaterEqual(min(g), 4.0 - 1e-6)
        self.assertGreater(len(g), 60)

    def test_injected_wobble_recovered(self):
        ns, _ = P.read_timings(DATA)
        rng = random.Random(7)
        ts = [58765.0 + 110.4 * (n - 17) - 0.109 * (n - 17) ** 2
              + 5.0 * math.cos(2 * math.pi * (n - 17) / 12.0 + 0.7) + rng.gauss(0.0, 1.0) for n in ns]
        r = P.analyse(ns, ts, nsim_fap=200, nsim_ul=60, ul_grid=[12.0], jackknife=False)
        self.assertEqual(r["outcome"], "WOBBLE_DETECTED")
        self.assertLess(abs(r["best"]["N"] - 12.0), 1.0)
        self.assertLess(abs(r["best"]["K"] - 5.0), 1.5)
        self.assertGreaterEqual(r["limits"][0]["K95"], r["limits"][0]["K_obs"])

    def test_upper_limit_covers_truth(self):
        # With no wobble, the 95% limit must sit above the fitted amplitude and be of order the noise.
        ns, _ = P.read_timings(DATA)
        rng = random.Random(11)
        ts = [58765.0 + 110.4 * (n - 17) - 0.109 * (n - 17) ** 2 + rng.gauss(0.0, 2.5) for n in ns]
        r = P.analyse(ns, ts, nsim_fap=100, nsim_ul=60, ul_grid=[16.0], jackknife=False)
        lim = r["limits"][0]
        self.assertGreater(lim["K95"], lim["K_obs"])
        self.assertLess(lim["K95"], 6.0)
        expect_f = lim["K95"] * (2 * math.pi / 16.0) ** 2 / (abs(r["base"]["Pdot"]) * r["base"]["P_days"])
        self.assertAlmostEqual(lim["f95"], expect_f, places=9)

    def test_data_file(self):
        ns, ts = P.read_timings(DATA)
        self.assertEqual(len(ns), 29)
        self.assertEqual(len(set(ns)), 29)
        self.assertEqual(ns, sorted(ns))
        self.assertTrue(all(t2 > t1 for t1, t2 in zip(ts, ts[1:])))

    @unittest.skipUnless(os.path.exists(os.path.join(HERE, "results", "v0.1.0", "results.json")),
                         "repo results file not present")
    def test_timings_match_v010_results(self):
        import json
        with open(os.path.join(HERE, "results", "v0.1.0", "results.json")) as fh:
            rj = json.load(fh)
        ref = {r["n"]: r["t"] for r in rj["graded"]["t2_full_baseline"]["residuals_steady"]}
        ns, ts = P.read_timings(DATA)
        self.assertEqual(sorted(ref), ns)
        for n, t in zip(ns, ts):
            self.assertAlmostEqual(ref[n], t, places=6)

    def test_real_pdot_matches_published_fit(self):
        # Equal-weight cubic fit should land near the weighted v0.1.0 fit (-0.00197 +- 0.00014).
        ns, ts = P.read_timings(DATA)
        b = P.fit_base(ns, ts)
        self.assertLess(b["Pdot"], -0.0014)
        self.assertGreater(b["Pdot"], -0.0026)
        self.assertGreater(b["P_days"], 109.5)
        self.assertLess(b["P_days"], 111.5)


if __name__ == "__main__":
    unittest.main()
