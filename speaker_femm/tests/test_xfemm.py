"""xfemm 版の検査。femmcli が無ければ実解析のテストはスキップ。

    FEMMCLI=/path/to/femmcli python -m unittest discover -s tests
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import femm_model as fm  # noqa: E402
import getdp_model as gm  # noqa: E402
import xfemm_model as xm  # noqa: E402

SKIP = xm.available()


class LuaTest(unittest.TestCase):
    def test_lua_script(self):
        p = fm.normalize({})
        G = fm.geometry(p)
        lua = xm.LuaFemm()
        fm.build_and_solve(lua, p, G, "/tmp")
        text = "\n".join(lua.lines)
        self.assertNotIn("mi_getmaterial", text)  # 材料ライブラリに頼らない
        self.assertEqual(text.count('mi_addbhpoint("1010 Steel"'), len(xm.BH_CURVES["1010 Steel"]))
        self.assertIn('mi_probdef(0, "millimeters", "axi"', text)
        self.assertIn("mi_analyze(1)", text)
        self.assertEqual(len(xm.collect_points(p, G)), len(xm.collect_points(p, G)))


@unittest.skipIf(SKIP, SKIP or "")
class XfemmTest(unittest.TestCase):
    def test_reference_case(self):
        r = xm.analyze(fm.DEFAULTS)
        self.assertTrue(r["ok"])
        self.assertAlmostEqual(r["summary"]["b_gap_center"], 0.971, delta=0.005)

    @unittest.skipIf(gm.available(), "getdp なし")
    def test_agrees_with_getdp(self):
        for c in ({}, {"magnet_r": 8, "plate_r": 10},
                  {"magnet_type": "Ferrite", "br": 0.4, "magnet_r": 15, "plate_r": 15, "magnet_t": 6}):
            x = xm.analyze(c)["summary"]
            g = gm.analyze(c)["summary"]
            for k in ("b_gap_center", "b_gap_avg", "bl0"):
                self.assertLess(abs(x[k] - g[k]) / abs(x[k]), 0.01, (c, k, x[k], g[k]))


if __name__ == "__main__":
    unittest.main()
