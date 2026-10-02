"""验证分析流程的基本检查：python -m unittest validation/test_validation.py（在仓库根目录运行）。"""
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyze  # noqa: E402
import vmodel as VM  # noqa: E402


class TestVModel(unittest.TestCase):
    def test_static_water_d80(self):
        # 25 °C、混施、壤土、f_soil=1 时，第 D 天应释放 80%
        F = VM.release_curve(np.ones(200), 90, 1.3, 0.03, 0.0)
        self.assertAlmostEqual(F[90], 0.8, places=6)

    def test_vg_roundtrip(self):
        vg = VM.VG_DEFAULT["loam"]
        for psi in (-10, -33, -300, -1500):
            self.assertAlmostEqual(float(VM.vg_psi(VM.vg_theta(psi, *vg), *vg)), psi, delta=abs(psi) * 1e-6)

    def test_dryness_matches_field_release(self):
        psi = np.array([0, -33, -200, -1500, -3000])
        np.testing.assert_allclose(1 - 0.25 * VM.dryness(psi), VM.FR.moisture_factor(psi, 0.25))


class TestAnalysis(unittest.TestCase):
    def test_example_recovers_parameters(self):
        out = Path(tempfile.mkdtemp())
        try:
            sites, products, obs, weather, soilT, moist, mgmt = analyze.load(HERE / "example")
            groups, _ = analyze.build_groups(sites, products, obs, weather, soilT, moist, mgmt, [])
            P, sup = analyze.calibrate(groups)
            # 示例数据的“真实”值：f_soil 0.80、条施 0.80、表施 2.0（见 make_example.py）
            self.assertAlmostEqual(P["f_soil"], 0.80, delta=0.08)
            self.assertAlmostEqual(P["band"], 0.80, delta=0.08)
            self.assertAlmostEqual(P["surface"], 2.0, delta=0.3)
            self.assertEqual(analyze.analyse(HERE / "example", out), 0)
            self.assertTrue((out / "report.md").exists())
        finally:
            shutil.rmtree(out)

    def test_empty_templates_rejected(self):
        with self.assertRaises(analyze.DataError):
            analyze.analyse(HERE / "templates", Path(tempfile.mkdtemp()))


if __name__ == "__main__":
    unittest.main()
