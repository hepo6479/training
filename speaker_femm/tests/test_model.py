"""モックで femm_model の形状生成と後処理を検査する (FEMM 不要)。

    python -m unittest discover -s tests
"""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import femm_model as fm  # noqa: E402
import mock_femm  # noqa: E402

CASES = [
    {},
    {"magnet_r": 8, "plate_r": 10},
    {"magnet_type": "Ferrite", "br": 0.4, "magnet_r": 15, "plate_r": 15, "magnet_t": 6},
    {"coil_h": 3, "coil_offset": -1.5, "gap": 0.6, "gap_mesh": 0.05},
]


def regions_ok(p):
    """線分をラスタ化し、ブロックラベルごとに閉領域がちょうど 1 つずつあるか調べる。"""
    G = fm.geometry(p)
    nodes, segs = fm._segments(G, with_axis=False)
    h = 0.02
    W, H = G["Ro"] + 3, G["zg1"] + 3
    nx, ny = int(W / h) + 1, int((H + 2) / h) + 1
    wall = np.zeros((ny, nx), bool)
    wall[:, 0] = True  # 対称軸
    def ij(x, y):
        return int(round((y + 2) / h)), int(round(x / h))
    for (x1, y1), (x2, y2) in segs:
        n = int(max(abs(x2 - x1), abs(y2 - y1)) / h * 2) + 2
        for t in np.linspace(0, 1, n):
            i, j = ij(x1 + (x2 - x1) * t, y1 + (y2 - y1) * t)
            wall[i, j] = True
    from scipy_free_label import label
    lab = label(~wall)
    labels = []
    mock_femm.calls.clear()
    fm.build_and_solve(mock_femm, p, G, "/nonexistent")
    for name, args in mock_femm.calls:
        if name == "mi_addblocklabel":
            x, y = args
            if x < W - 0.1:
                labels.append(lab[ij(x, y)])
    outer = lab[ij(W - 0.05, H / 2)]
    # 全領域 = 外側空気 + 内部ラベル。ラベル同士が重複せず、ラベルの無い閉領域がないこと
    regions = set(np.unique(lab)) - {0}
    return labels, outer, regions


class T(unittest.TestCase):
    def test_regions(self):
        for c in CASES:
            p = fm.normalize(c)
            self.assertEqual(fm.validate(p), [])
            fm.analyze(p, femm=mock_femm)  # configure
            labels, outer, regions = regions_ok(p)
            self.assertNotIn(0, labels, c)
            self.assertEqual(len(labels), len(set(labels)), f"ラベル重複 {c}")
            self.assertEqual(set(labels) | {outer}, regions, f"ラベル無し領域 {c}")

    def test_analyze(self):
        for c in CASES:
            r = fm.analyze(c, femm=mock_femm)
            self.assertTrue(r["ok"])
            m = r["map"]
            self.assertEqual(len(m["br"]), m["nr"] * m["nz"])
            self.assertEqual(len(r["bl_x"]["x"]), len(r["bl_x"]["bl"]))

    def test_validate(self):
        self.assertTrue(fm.validate(fm.normalize({"magnet_r": 12, "plate_r": 10})))
        self.assertTrue(fm.validate(fm.normalize({"gap": -1})))
        self.assertTrue(fm.validate(fm.normalize({"gap_mesh": 0.8})))
        r = fm.analyze({"gap": 0}, femm=mock_femm)
        self.assertFalse(r["ok"])


if __name__ == "__main__":
    unittest.main()
