"""基本检查：python -m unittest discover -s tests（在 rice-crf-map 目录下运行）。"""
import json
import unittest
from pathlib import Path

import numpy as np

from crfmap.data import COUNTRIES, ZONES
from crfmap.demand import season_setup
from crfmap.optimize import best_recipe
from crfmap.params import CROP, EST, P, PRODUCTS
from crfmap.release import arrhenius, release

HERE = Path(__file__).resolve().parent.parent


class ReleaseModel(unittest.TestCase):
    def test_nominal_period_gives_80_percent(self):
        for D in PRODUCTS:
            self.assertAlmostEqual(float(release(np.array([D]), D)[0]), 0.80, places=9)

    def test_arrhenius_reference_and_q10(self):
        self.assertAlmostEqual(float(arrhenius(25.0, P["Ea"])), 1.0, places=12)
        q10 = float(arrhenius(30.0, P["Ea"]) / arrhenius(20.0, P["Ea"]))
        self.assertTrue(1.9 < q10 < 2.1, q10)

    def test_release_is_monotonic(self):
        tau = np.linspace(0, 400, 401)
        for D in (30, 90, 360):
            self.assertTrue(np.all(np.diff(release(tau, D)) >= 0))


class DataFiles(unittest.TestCase):
    def test_zone_and_season_counts(self):
        self.assertEqual(len(ZONES), 50)
        self.assertEqual(sum(len(z["seasons"]) for z in ZONES), 73)

    def test_codes_are_known(self):
        ids = [z["id"] for z in ZONES]
        self.assertEqual(len(ids), len(set(ids)))
        for z in ZONES:
            self.assertIn(z["country"], COUNTRIES)
            self.assertEqual(len(z["T"]), 12)
            self.assertIn(z["crop"], {"single", "double", "triple"})
            self.assertTrue(z["seasons"], z["id"])
            for s in z["seasons"]:
                self.assertIn(s["type"], CROP)
                self.assertIn(s["est"], EST)


class Regression(unittest.TestCase):
    """与已发布 data.json 对照：长江中游一季中稻（武汉）。"""

    def test_wuhan_single_rice(self):
        z = next(z for z in ZONES if z["id"] == "cn-myz")
        main, *_ = best_recipe(season_setup(z, z["seasons"][0]))
        pub = json.loads((HERE / "data.json").read_text())
        s = next(z for z in pub["zones"] if z["id"] == "cn-myz")["seasons"][0]
        self.assertEqual(main["total"], s["total"])
        self.assertEqual([list(r) for r in main["rec"]], s["rec"])


if __name__ == "__main__":
    unittest.main()
