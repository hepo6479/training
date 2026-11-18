"""内磁型スピーカ磁気回路の軸対称静磁場解析 (FEMM 4.2 + pyfemm)

analyze(params) -> dict (JSON 化可能)

座標系: r = 半径方向 (FEMM の x), z = 軸方向 (FEMM の y)。単位 mm。
    z = 0              ヨーク底面
    z = tb             ヨーク底の上面 = マグネット下面
    z = tb + hm        マグネット上面 = プレート下面
    z = tb + hm + tp   プレート上面 = ヨーク側壁上端 (面一)
ヨーク内径 Ri = Rp + g, 外径 Ro = Ri + tw。マグネットは +z 方向に着磁。
"""

import math
import os
import tempfile
import time

MU0 = 4e-7 * math.pi

IRON_MATERIALS = ["1010 Steel", "1006 Steel", "1018 Steel", "1020 Steel", "1117 Steel"]
MAGNET_TYPES = {
    # mu_r, 導電率 [MS/m], 密度 [g/cm3], Br の妥当範囲
    "NdFeB": {"mur": 1.05, "cond": 0.667, "density": 7.5, "br_range": (0.9, 1.55)},
    "Ferrite": {"mur": 1.10, "cond": 0.0, "density": 4.9, "br_range": (0.2, 0.5)},
}
IRON_DENSITY = 7.87  # g/cm3
SAT_WARN = 1.5   # T  これ以上で「注意」
SAT_LIMIT = 1.8  # T  これ以上で「飽和」
XMAX_RATIO = 0.82

DEFAULTS = {
    "magnet_type": "NdFeB",
    "br": 1.3,            # T
    "magnet_r": 10.0,     # mm (φ20)
    "magnet_t": 5.0,      # mm
    "plate_r": 10.0,      # mm
    "plate_t": 3.0,       # mm
    "gap": 1.0,           # mm
    "yoke_wall": 2.0,     # mm
    "yoke_bottom": 2.5,   # mm
    "iron": "1010 Steel",
    "coil_h": 6.0,        # mm 巻き高さ
    "coil_turns": 60,     # 巻数
    "coil_offset": 0.0,   # mm 静止位置のずれ (+ = 上)
    "gap_mesh": 0.1,      # mm
}


# ---------------------------------------------------------------- 入力検証
def normalize(params):
    p = dict(DEFAULTS)
    for k, v in (params or {}).items():
        if k not in DEFAULTS:
            continue
        if isinstance(DEFAULTS[k], str):
            p[k] = str(v)
        elif k == "coil_turns":
            p[k] = int(round(float(v)))
        else:
            p[k] = float(v)
    return p


def validate(p):
    """不正な形状のエラーメッセージのリストを返す (空なら OK)。index.html と同じ規則。"""
    e = []
    num = ["br", "magnet_r", "magnet_t", "plate_r", "plate_t", "gap", "yoke_wall",
           "yoke_bottom", "coil_h", "coil_turns", "gap_mesh"]
    for k in num:
        v = p[k]
        if not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0:
            e.append(f"{k} は正の数にしてください")
    if not math.isfinite(p["coil_offset"]):
        e.append("coil_offset が数値ではありません")
    if e:
        return e
    if p["magnet_type"] not in MAGNET_TYPES:
        e.append("マグネット種類が不正です")
    else:
        lo, hi = MAGNET_TYPES[p["magnet_type"]]["br_range"]
        if not lo <= p["br"] <= hi:
            e.append(f"{p['magnet_type']} の Br は {lo}〜{hi} T の範囲にしてください")
    if p["iron"] not in IRON_MATERIALS:
        e.append("鉄材質が不正です")
    if p["magnet_r"] > p["plate_r"]:
        e.append("マグネット半径はプレート半径以下にしてください (内磁型)")
    if p["gap_mesh"] > p["gap"] / 2:
        e.append("ギャップのメッシュサイズはギャップ幅の 1/2 以下にしてください")
    if p["gap_mesh"] < 0.01:
        e.append("ギャップのメッシュサイズは 0.01 mm 以上にしてください")
    if p["coil_h"] > 50:
        e.append("巻き高さが大きすぎます (50 mm 以下)")
    if abs(p["coil_offset"]) > p["coil_h"] + p["plate_t"]:
        e.append("静止位置のずれが大きすぎます (|ずれ| ≤ 巻き高さ + プレート厚)")
    if p["coil_turns"] > 10000:
        e.append("巻数が大きすぎます")
    return e


# ---------------------------------------------------------------- 形状
def geometry(p):
    tb, hm, tp = p["yoke_bottom"], p["magnet_t"], p["plate_t"]
    Rm, Rp, g, tw = p["magnet_r"], p["plate_r"], p["gap"], p["yoke_wall"]
    zpb = tb + hm
    ztop = zpb + tp
    Ri = Rp + g
    Ro = Ri + tw
    zcoil = (zpb + ztop) / 2 + p["coil_offset"]
    # Bl(x) を求める変位範囲
    xr = round(max(4.0, 1.5 * max(p["coil_h"], tp)), 1)
    zs0 = max(tb + 0.05, zcoil - xr - p["coil_h"] / 2)
    zs1 = zcoil + xr + p["coil_h"] / 2
    zg1 = max(ztop + 1.0, zs1 + 0.5)   # 細分割ギャップ領域の上端
    size = max(Ro, zg1)
    return {
        "tb": tb, "zpb": zpb, "ztop": ztop, "Rm": Rm, "Rp": Rp, "Ri": Ri, "Ro": Ro,
        "rgap": Rp + g / 2, "zcoil": zcoil, "xr": xr, "zs0": zs0, "zs1": zs1,
        "zg1": zg1, "R_out": 5.0 * size, "zc": ztop / 2,
    }


def part_at(G, r, z):
    """点 (r, z) の部位名 (r は絶対値)"""
    r = abs(r)
    if 0 <= z <= G["tb"] and r <= G["Ro"]:
        return "yoke_bottom"
    if G["tb"] < z <= G["ztop"] and G["Ri"] <= r <= G["Ro"]:
        return "yoke_wall"
    if G["tb"] < z <= G["zpb"] and r <= G["Rm"]:
        return "magnet"
    if G["zpb"] < z <= G["ztop"] and r <= G["Rp"]:
        return "plate"
    return "air"


def masses(p, G):
    def cyl(r, h):
        return math.pi * r * r * h / 1000.0  # cm3
    v_yoke = cyl(G["Ro"], G["tb"]) + (cyl(G["Ro"], 1) - cyl(G["Ri"], 1)) * (G["ztop"] - G["tb"])
    v_plate = cyl(G["Rp"], p["plate_t"])
    v_mag = cyl(G["Rm"], p["magnet_t"])
    m = {
        "yoke": v_yoke * IRON_DENSITY,
        "plate": v_plate * IRON_DENSITY,
        "magnet": v_mag * MAGNET_TYPES[p["magnet_type"]]["density"],
    }
    m["total"] = sum(m.values())
    return {k: round(v, 2) for k, v in m.items()}


def _segments(G, with_axis=True):
    """全領域の境界線分を、交点で分割・重複除去して返す。"""
    tb, zpb, ztop = G["tb"], G["zpb"], G["ztop"]
    Rm, Rp, Ri, Ro, zg1 = G["Rm"], G["Rp"], G["Ri"], G["Ro"], G["zg1"]

    def poly(pts):
        return [(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]

    edges = []
    edges += poly([(0, 0), (Ro, 0), (Ro, ztop), (Ri, ztop), (Ri, tb), (0, tb)])     # ヨーク
    edges += poly([(0, tb), (Rm, tb), (Rm, zpb), (0, zpb)])                          # マグネット
    edges += poly([(0, zpb), (Rp, zpb), (Rp, ztop), (0, ztop)])                      # プレート
    edges += poly([(Rp, tb), (Ri, tb), (Ri, zg1), (Rp, zg1)])                        # ギャップ細分割領域
    if with_axis:
        edges.append(((0, G["zc"] - G["R_out"]), (0, G["zc"] + G["R_out"])))

    def key(pt):
        return (round(pt[0], 9), round(pt[1], 9))

    pts = {key(a) for e in edges for a in e}
    out = set()
    for a, b in edges:
        a, b = key(a), key(b)
        if a == b:
            continue
        dx, dy = b[0] - a[0], b[1] - a[1]
        L2 = dx * dx + dy * dy
        on = []
        for q in pts:
            cross = (q[0] - a[0]) * dy - (q[1] - a[1]) * dx
            if abs(cross) > 1e-9 * math.sqrt(L2):
                continue
            t = ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / L2
            if -1e-12 <= t <= 1 + 1e-12:
                on.append((t, q))
        on.sort()
        for (t1, q1), (t2, q2) in zip(on, on[1:]):
            if q1 != q2:
                out.add(tuple(sorted((q1, q2))))
    nodes = sorted({q for s in out for q in s})
    return nodes, sorted(out)


# ---------------------------------------------------------------- FEMM モデル
def build_and_solve(femm, p, G, workdir):
    femm.newdocument(0)
    femm.mi_probdef(0, "millimeters", "axi", 1e-8, 0, 30)

    femm.mi_getmaterial("Air")
    femm.mi_getmaterial(p["iron"])
    mt = MAGNET_TYPES[p["magnet_type"]]
    hc = p["br"] / (MU0 * mt["mur"])
    femm.mi_addmaterial("Magnet", mt["mur"], mt["mur"], hc, 0, mt["cond"], 0, 0, 1, 0, 0, 0)
    femm.mi_addboundprop("A0", 0, 0, 0, 0, 0, 0, 0, 0, 0)

    nodes, segs = _segments(G)
    for x, y in nodes:
        femm.mi_addnode(x, y)
    for (x1, y1), (x2, y2) in segs:
        femm.mi_addsegment(x1, y1, x2, y2)

    # 外側境界: A=0 の半円
    R, zc = G["R_out"], G["zc"]
    femm.mi_addarc(0, zc - R, 0, zc + R, 180, 2)
    femm.mi_selectarcsegment(R, zc)
    femm.mi_setarcsegmentprop(2, "A0", 0, 0)
    femm.mi_clearselected()

    iron_mesh = max(p["gap_mesh"], min(0.5, min(p["yoke_wall"], p["yoke_bottom"], p["plate_t"]) / 4))
    tb, zpb, ztop = G["tb"], G["zpb"], G["ztop"]

    def label(x, y, mat, automesh, mesh, magdir=0, group=0):
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop(mat, automesh, mesh, "", magdir, group, 0)
        femm.mi_clearselected()

    label((G["Ri"] + G["Ro"]) / 2, tb / 2, p["iron"], 0, iron_mesh, group=1)          # ヨーク
    label(G["Rm"] / 2, (tb + zpb) / 2, "Magnet", 0, max(iron_mesh, 0.5), 90, group=2)  # マグネット (+z)
    label(G["Rp"] / 2, (zpb + ztop) / 2, p["iron"], 0, iron_mesh, group=3)            # プレート
    label(G["rgap"], (tb + G["zg1"]) / 2, "Air", 0, p["gap_mesh"], group=4)           # ギャップ
    if G["Rm"] < G["Rp"]:
        label((G["Rm"] + G["Rp"]) / 2, (tb + zpb) / 2, "Air", 1, 0, group=5)          # プレート下の空気
    label(G["Ro"] + 0.3 * (R - G["Ro"]), zc, "Air", 1, 0, group=0)                    # 外側空気

    path = os.path.join(workdir, "speaker.fem")
    femm.mi_saveas(path)
    femm.mi_analyze(1)
    femm.mi_loadsolution()

    def getb(r, z):
        b = femm.mo_getb(r, z)
        return float(b[0]), float(b[1])

    return getb


def close_docs(femm):
    for fn in ("mo_close", "mi_close"):
        try:
            getattr(femm, fn)()
        except Exception:
            pass


# ---------------------------------------------------------------- 後処理 (FEMM/モック共通)
def postprocess(p, G, getb):
    out = {}
    rg = G["rgap"]

    # ギャップ中心の Br(z)
    dz = 0.1
    n = int(round((G["zs1"] - G["zs0"]) / dz)) + 1
    zs = [G["zs0"] + i * dz for i in range(n)]
    brs = [getb(rg, z)[0] for z in zs]
    out["br_z"] = {"r": round(rg, 4), "z": [round(z, 3) for z in zs], "br": [round(b, 5) for b in brs]}

    def br_at(z):
        if z < zs[0]:
            return 0.0  # ヨーク底より下 (物理的に到達不可)
        if z >= zs[-1]:
            return brs[-1]
        i = int((z - zs[0]) / dz)
        i = min(i, n - 2)
        t = (z - zs[i]) / dz
        return brs[i] * (1 - t) + brs[i + 1] * t

    # Bl(x): 1 回の解析結果からコイル位置をずらして計算
    h, N = p["coil_h"], p["coil_turns"]
    circ = 2 * math.pi * rg / 1000.0  # m

    def bl(x):
        z0 = G["zcoil"] + x - h / 2
        k = max(20, int(h / 0.05))
        s = 0.0
        for j in range(k + 1):
            w = 0.5 if j in (0, k) else 1.0
            s += w * br_at(z0 + h * j / k)
        mean_br = s / k
        return mean_br * N * circ

    dx = 0.1
    nx = int(round(G["xr"] / dx))
    xs = [i * dx for i in range(-nx, nx + 1)]
    bls = [bl(x) for x in xs]
    bl0 = bls[nx]
    out["bl_x"] = {"x": [round(x, 2) for x in xs], "bl": [round(b, 5) for b in bls]}

    def find_xmax(direction):
        thr = XMAX_RATIO * abs(bl0)
        prev_x, prev_b = 0.0, abs(bl0)
        for i in range(1, nx + 1):
            x = direction * i * dx
            b = abs(bls[nx + direction * i])
            if b < thr:
                t = (prev_b - thr) / (prev_b - b) if prev_b != b else 0
                return abs(prev_x + (x - prev_x) * t)
            prev_x, prev_b = x, b
        return None  # 範囲内では落ちない

    xp, xm = find_xmax(+1), find_xmax(-1)
    cands = [v for v in (xp, xm) if v is not None]
    xmax = min(cands) if cands else None

    # ギャップ磁束密度
    b_center = getb(rg, (G["zpb"] + G["ztop"]) / 2)
    k = 30
    b_avg = sum(getb(rg, G["zpb"] + p["plate_t"] * (j + 0.5) / k)[0] for j in range(k)) / k

    # |B| マップ
    margin = max(2.0, 0.25 * max(G["Ro"], G["ztop"]))
    r0, r1 = 0.0, G["Ro"] + margin
    z0, z1 = -margin, max(G["ztop"] + margin, G["zs1"])
    step = max(0.1, math.sqrt((r1 - r0) * (z1 - z0) / 4500.0))
    nr = int(math.ceil((r1 - r0) / step)) + 1
    nzm = int(math.ceil((z1 - z0) / step)) + 1
    mbr, mbz = [], []
    parts = {"yoke_bottom": [], "yoke_wall": [], "plate": [], "magnet": []}
    for j in range(nzm):
        z = z0 + j * step
        for i in range(nr):
            r = r0 + i * step
            br, bz = getb(r, z)
            mbr.append(round(br, 4))
            mbz.append(round(bz, 4))
            pt = part_at(G, r, z)
            if pt in parts:
                parts[pt].append(math.hypot(br, bz))
    out["map"] = {"r0": r0, "z0": round(z0, 4), "step": round(step, 5), "nr": nr, "nz": nzm,
                  "br": mbr, "bz": mbz}

    # 鉄の各部位 |B| と飽和判定 (角の特異点を避けるため 95 パーセンタイルで判定)
    names = {"yoke_bottom": "ヨーク底", "yoke_wall": "ヨーク側壁", "plate": "トッププレート"}
    iron = []
    for key in ("plate", "yoke_wall", "yoke_bottom"):
        vals = sorted(parts[key])
        if not vals:
            iron.append({"part": names[key], "key": key, "n": 0})
            continue
        b95 = vals[min(len(vals) - 1, int(0.95 * len(vals)))]
        status = "OK" if b95 < SAT_WARN else ("注意" if b95 < SAT_LIMIT else "飽和")
        iron.append({"part": names[key], "key": key, "n": len(vals),
                     "b_mean": round(sum(vals) / len(vals), 3), "b95": round(b95, 3),
                     "b_max": round(vals[-1], 3), "status": status})
    out["iron"] = iron

    wire_len = N * circ
    out["summary"] = {
        "b_gap_center": round(b_center[0], 4),
        "b_gap_avg": round(b_avg, 4),
        "bl0": round(bl0, 4),
        "bl_peak": round(max(bls, key=abs), 4),
        "xmax": None if xmax is None else round(xmax, 2),
        "xmax_plus": None if xp is None else round(xp, 2),
        "xmax_minus": None if xm is None else round(xm, 2),
        "coil_radius": round(rg, 3),
        "wire_length_m": round(wire_len, 3),
        "overhang": round((p["coil_h"] - p["plate_t"]) / 2, 3),
    }
    out["mass"] = masses(p, G)
    return out


# ---------------------------------------------------------------- 入口
def analyze(params, femm=None, keep_open=False):
    """解析を実行して結果 dict を返す。femm にはモジュール (pyfemm かモック) を渡せる。"""
    t0 = time.time()
    p = normalize(params)
    errs = validate(p)
    if errs:
        return {"ok": False, "errors": errs, "params": p}
    G = geometry(p)
    if femm is None:
        import femm  # pyfemm
    if hasattr(femm, "configure"):  # モックバックエンド用
        femm.configure(p, G)
    opened = False
    try:
        if not keep_open:
            femm.openfemm(1)
            opened = True
        with tempfile.TemporaryDirectory(prefix="spk_femm_") as wd:
            getb = build_and_solve(femm, p, G, wd)
            res = postprocess(p, G, getb)
            close_docs(femm)
    finally:
        if opened:
            try:
                femm.closefemm()
            except Exception:
                pass
    res.update({
        "ok": True,
        "backend": getattr(femm, "BACKEND_NAME", "femm"),
        "params": p,
        "geometry": {k: round(v, 5) for k, v in G.items()},
        "elapsed_s": round(time.time() - t0, 2),
    })
    return res


if __name__ == "__main__":
    import json
    import sys
    r = analyze(DEFAULTS)
    json.dump({k: v for k, v in r.items() if k != "map"}, sys.stdout, ensure_ascii=False, indent=1)
