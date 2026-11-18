"""FEMM が無い環境 (Linux / CI) で UI と後処理を動かすためのモック。

pyfemm と同じ関数名を持ち、呼び出しを記録する。磁界は簡易磁気回路モデルから
それらしい値を合成するだけで、FEM 解析ではない。数値は目安にもしないこと。
"""

import math

BACKEND_NAME = "mock"
calls = []
_state = {}


def _rec(name):
    def f(*args):
        calls.append((name, args))
    return f


for _n in ["openfemm", "closefemm", "newdocument", "mi_probdef", "mi_getmaterial",
           "mi_addmaterial", "mi_addboundprop", "mi_addnode", "mi_addsegment", "mi_addarc",
           "mi_selectarcsegment", "mi_setarcsegmentprop", "mi_clearselected",
           "mi_addblocklabel", "mi_selectlabel", "mi_setblockprop", "mi_saveas",
           "mi_analyze", "mi_loadsolution", "mo_close", "mi_close"]:
    globals()[_n] = _rec(_n)


def configure(p, G):
    calls.clear()
    Am = math.pi * G["Rm"] ** 2
    Ag = 2 * math.pi * G["rgap"] * p["plate_t"]
    mur = 1.05 if p["magnet_type"] == "NdFeB" else 1.1
    sigma = 1.9  # 漏れ係数 (仮定)
    Bm = p["br"] / (1 + mur * p["gap"] * Am / (sigma * Ag * p["magnet_t"]))
    phi = Bm * Am
    _state.update(p=p, G=G, Bm=Bm, phi=phi, Bg=phi / (sigma * Ag), sigma=sigma)


def _fringe(z, G, w):
    a, b = G["zpb"], G["ztop"]
    norm = 2 * math.atan((b - a) / (2 * w))
    return (math.atan((z - a) / w) - math.atan((z - b) / w)) / norm


def mo_getb(r, z):
    s, G, p = _state, _state["G"], _state["p"]
    r = abs(r)
    Bm, phi, Bg = s["Bm"], s["phi"], s["Bg"]
    tb, zpb, ztop = G["tb"], G["zpb"], G["ztop"]
    Rm, Rp, Ri, Ro = G["Rm"], G["Rp"], G["Ri"], G["Ro"]
    tp = p["plate_t"]
    phi_w = phi / 1.3
    Aw = math.pi * (Ro ** 2 - Ri ** 2)
    re = max(r, 1e-6)
    if 0 <= z <= tb and r <= Ro:                       # ヨーク底
        rr = min(r, Rm)
        flow = phi_w * rr * rr / (Rm * Rm) if r <= Ri else phi_w * (Ro - r) / (Ro - Ri)
        return -flow / (2 * math.pi * re * tb), Bm * (z / tb) * (1 if r <= Rm else 0.3)
    if tb < z <= ztop and Ri <= r <= Ro:               # 側壁
        k = min(1.0, max(0.2, (ztop - z) / tp))
        return -0.15 * Bg * k, -phi_w / Aw * k
    if tb < z <= zpb and r <= Rm:                      # マグネット
        return 0.0, Bm
    if zpb < z <= ztop and r <= Rp:                    # プレート
        rr = min(r, Rm)
        br = Bm * rr * rr / (2 * re * tp) / math.sqrt(s["sigma"]) * (z - zpb + 0.5) / (tp + 0.5) * 1.6
        bz = Bm * max(0.0, (ztop - z) / tp) * (1 if r <= Rm else 0.5)
        return br, bz
    # 空気
    w = 0.6 * p["gap"]
    f = _fringe(z, G, w)
    d = max(0.0, abs(r - G["rgap"]) - p["gap"] / 2)
    br = Bg * f * (G["rgap"] / re if r >= Rp else r / Rp) * math.exp(-d / 1.5)
    bz = 0.0
    if z > ztop:
        bz = 0.25 * Bg * math.exp(-(z - ztop) / 2) * (1 - min(1, r / Ro))
    if z < 0:
        bz = -0.1 * Bg * math.exp(z / 2) * (1 - min(1, r / Ro))
    return br, bz
