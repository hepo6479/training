"""GetDP 版の検証 (gmsh と getdp が無ければスキップ)

    python -m unittest discover -s tests
"""
import math
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import femm_model as fm  # noqa: E402
import getdp_model as gm  # noqa: E402

SKIP = gm.available()


@unittest.skipIf(SKIP, SKIP or "")
class GetDPTest(unittest.TestCase):
    def test_magnet_in_air_matches_analytic(self):
        """鉄を空気、mu_r = 1 にした円柱磁石の軸上 Bz を解析解と比べる (定式化の検証)"""
        p = fm.normalize({})
        G = fm.geometry(p)
        wd = tempfile.mkdtemp()
        try:
            nodes, tris = gm.build_mesh(p, G, os.path.join(wd, "speaker.msh"))
            gm.write_pro(p, os.path.join(wd, "speaker.pro"))
            path = os.path.join(wd, "speaker.pro")
            s = open(path).read()
            s = s.replace("nu[Iron] = nu_fe[$1];", "nu[Iron] = 1/mu0;")
            s = s.replace("dhdb_NL[Iron] = 2*dnudb2_fe[$1#1]*SquDyadicProduct[#1];",
                          "dhdb_NL[Iron] = 0*TensorDiag[1,1,1];")
            s = s.replace("nu[Magnet] = 1/(mu0*1.05);", "nu[Magnet] = 1/mu0;")
            open(path, "w").write(s)
            gm.run_getdp(wd)
            f = gm.Field(nodes, tris, gm.read_b_elements(os.path.join(wd, "b.txt"), tris, nodes))
        finally:
            shutil.rmtree(wd, ignore_errors=True)
        Br, a, L, zc = p["br"], p["magnet_r"], p["magnet_t"], G["tb"] + p["magnet_t"] / 2
        for dz in (0, 2, 4, 8, -6):
            u = dz
            exact = Br / 2 * ((u + L / 2) / math.hypot(u + L / 2, a) - (u - L / 2) / math.hypot(u - L / 2, a))
            self.assertAlmostEqual(f.getb(0.05, zc + dz)[1], exact, delta=0.02 * Br)

    def test_reference_case(self):
        r = gm.analyze(fm.DEFAULTS)
        self.assertTrue(r["ok"])
        s = r["summary"]
        # 基準値: φ20×5, Br 1.3 T, プレート 3, ギャップ 1, 側壁 2, 底 2.5 mm
        # GetDP 結果 0.97 T (仕様の目安 0.92 T と約 5% 差。README 参照)
        self.assertTrue(0.90 <= s["b_gap_center"] <= 1.0, s["b_gap_center"])
        self.assertIsNotNone(s["xmax"])

    def test_mesh_convergence(self):
        a = gm.analyze({"gap_mesh": 0.1})["summary"]["b_gap_center"]
        b = gm.analyze({"gap_mesh": 0.05})["summary"]["b_gap_center"]
        self.assertLess(abs(a - b) / b, 0.01)

    def test_other_cases_run(self):
        for c in ({"magnet_r": 8, "plate_r": 10}, {"magnet_type": "Ferrite", "br": 0.4, "magnet_r": 15,
                  "plate_r": 15, "magnet_t": 6}, {"iron": "1018 Steel", "coil_offset": -1}):
            r = gm.analyze(c)
            self.assertTrue(r["ok"], c)


if __name__ == "__main__":
    unittest.main()
