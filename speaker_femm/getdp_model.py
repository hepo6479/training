"""GetDP + Gmsh 版の解析部 (macOS / Linux / Windows)

femm_model と同じモデル・同じ後処理で、FEM 解析だけを Gmsh (メッシュ) と
GetDP (ソルバ) に置き換える。analyze(params) は femm_model.analyze と同じ形の dict を返す。

    pip install gmsh            # メッシュ生成 (Python API)
    GetDP: https://getdp.info から実行ファイルを入手 (Linux は apt install getdp でも可)
           PATH に通すか、環境変数 GETDP に実行ファイルのパスを入れる

定式化: 軸対称のベクトルポテンシャル A_phi (BF_PerpendicularEdge + VolAxiSqu)。
鉄は FEMM ライブラリと同じ B-H カーブ (bh_data.py) の非線形材料として Newton 法で解く。
マグネットは H = (B - Br)/(mu0 mu_r) の線形材料。外側は A = 0 の半円 (と対称軸)。
"""

import math
import os
import shutil
import subprocess
import tempfile
import time
from bisect import bisect_right

import femm_model as fm
from bh_data import BH_CURVES

BACKEND_NAME = "getdp"
PHYS_IRON, PHYS_MAGNET, PHYS_AIR, PHYS_DIR = 101, 102, 103, 201


def find_getdp():
    exe = os.environ.get("GETDP") or shutil.which("getdp")
    if exe and os.path.exists(exe):
        return exe
    for c in ("/Applications/getdp/bin/getdp", "/usr/local/bin/getdp", "/opt/homebrew/bin/getdp",
              os.path.expanduser("~/getdp/bin/getdp"), r"C:\getdp\getdp.exe"):
        if os.path.exists(c):
            return c
    return None


def available():
    """使えなければ理由の文字列、使えれば None"""
    try:
        import gmsh  # noqa: F401
    except Exception as e:
        return f"gmsh (Python) が使えません: {e}"
    if not find_getdp():
        return "getdp の実行ファイルが見つかりません (PATH か環境変数 GETDP で指定)"
    return None


# ---------------------------------------------------------------- B-H
def b2_nu_table(iron):
    """GetDP の InterpolationLinear 用 (B^2, nu) 表。最終点より先は傾き mu0 で延長 (FEMM と同じ)"""
    pts = sorted(BH_CURVES[iron])
    bh = [pts[0]]
    for b, h in pts[1:]:
        if b > bh[-1][0] and h > bh[-1][1]:
            bh.append((b, h))
    bl, hl = bh[-1]
    for k in range(1, 9):
        b = bl + 0.5 * k * k
        bh.append((b, hl + (b - bl) / fm.MU0))
    tab = []
    for b, h in bh[1:]:
        tab.append((b * b, h / b))
    tab.insert(0, (0.0, tab[0][1]))
    return tab


# ---------------------------------------------------------------- メッシュ
def build_mesh(p, G, path):
    import gmsh
    # interruptible=False: サーバのワーカースレッド (メインスレッド以外) から呼べるようにする
    gmsh.initialize(["-noenv"], interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("speaker")
        occ = gmsh.model.occ
        R, zc = G["R_out"], G["zc"]
        tb, zpb, ztop = G["tb"], G["zpb"], G["ztop"]
        Rm, Rp, Ri, Ro, zg1 = G["Rm"], G["Rp"], G["Ri"], G["Ro"], G["zg1"]

        disk = occ.addDisk(0, zc, 0, R, R)
        box = occ.addRectangle(0, zc - R - 1, 0, R + 1, 2 * R + 2)
        half, _ = occ.intersect([(2, disk)], [(2, box)])
        pieces = {
            "iron": [occ.addRectangle(0, 0, 0, Ro, tb), occ.addRectangle(Ri, tb, 0, Ro - Ri, ztop - tb),
                     occ.addRectangle(0, zpb, 0, Rp, ztop - zpb)],
            "magnet": [occ.addRectangle(0, tb, 0, Rm, zpb - tb)],
            "gap": [occ.addRectangle(Rp, tb, 0, Ri - Rp, zg1 - tb)],
        }
        tools = [(2, t) for v in pieces.values() for t in v]
        _, omap = occ.fragment(half, tools)
        occ.synchronize()
        children = [set(t for d, t in m if d == 2) for m in omap]
        allsurf = set().union(*children)
        k = len(half)
        idx = k
        sets = {}
        for name, tags in pieces.items():
            sets[name] = set()
            for _ in tags:
                sets[name] |= children[idx]
                idx += 1
        iron, magnet = sets["iron"], sets["magnet"]
        air = allsurf - iron - magnet
        gmsh.model.addPhysicalGroup(2, sorted(iron), PHYS_IRON)
        gmsh.model.addPhysicalGroup(2, sorted(magnet), PHYS_MAGNET)
        gmsh.model.addPhysicalGroup(2, sorted(air), PHYS_AIR)
        bnd = gmsh.model.getBoundary([(2, t) for t in sorted(allsurf)], combined=True, oriented=False)
        gmsh.model.addPhysicalGroup(1, sorted({abs(t) for d, t in bnd}), PHYS_DIR)

        # メッシュサイズ: ギャップを細かく、鉄・マグネットは中程度、外側は粗く
        iron_mesh = max(p["gap_mesh"], min(0.5, min(p["yoke_wall"], p["yoke_bottom"], p["plate_t"]) / 4))
        big = R / 12
        f = gmsh.model.mesh.field
        f1 = f.add("Box")
        for kk, v in (("VIn", p["gap_mesh"]), ("VOut", big), ("XMin", Rp), ("XMax", Ri),
                      ("YMin", tb), ("YMax", zg1), ("Thickness", 2.0)):
            f.setNumber(f1, kk, v)
        f2 = f.add("Box")
        for kk, v in (("VIn", iron_mesh), ("VOut", big), ("XMin", 0), ("XMax", Ro + 1),
                      ("YMin", -1), ("YMax", max(ztop, zg1) + 1), ("Thickness", 0.6 * max(Ro, ztop))):
            f.setNumber(f2, kk, v)
        f3 = f.add("Min")
        f.setNumbers(f3, "FieldsList", [f1, f2])
        f.setAsBackgroundMesh(f3)
        for kk, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                      ("Mesh.MeshSizeFromCurvature", 0), ("Mesh.Algorithm", 6),
                      ("Mesh.MshFileVersion", 2.2)):
            gmsh.option.setNumber(kk, v)
        gmsh.model.mesh.generate(2)
        gmsh.write(path)

        # 後処理用に節点と三角形 (要素番号つき) を取り出す
        ntags, coords, _ = gmsh.model.mesh.getNodes()
        nodes = {int(t): (coords[3 * i], coords[3 * i + 1]) for i, t in enumerate(ntags)}
        tris = {}
        for reg, surfs in ((PHYS_IRON, iron), (PHYS_MAGNET, magnet), (PHYS_AIR, air)):
            for s in surfs:
                types, etags, enodes = gmsh.model.mesh.getElements(2, s)
                for ty, et, en in zip(types, etags, enodes):
                    if ty != 2:
                        continue
                    for i, e in enumerate(et):
                        tris[int(e)] = (reg, int(en[3 * i]), int(en[3 * i + 1]), int(en[3 * i + 2]))
        return nodes, tris
    finally:
        gmsh.finalize()


# ---------------------------------------------------------------- GetDP
PRO = r"""
Group {
  Iron = Region[%(iron)d]; Magnet = Region[%(magnet)d]; Air = Region[%(air)d];
  Vol = Region[{Iron, Magnet, Air}];
  Vol_NL = Region[{Iron}];
  Dir = Region[%(dir)d];
}
Function {
  mu0 = 4e-7*Pi;
  nu[Air] = 1/mu0;
  nu[Magnet] = 1/(mu0*%(mur).10g);
  br[Magnet] = Vector[0, %(br).10g, 0];
  Mat_b2() = {%(b2)s};
  Mat_nu() = {%(nu)s};
  Mat_b2_nu() = ListAlt[Mat_b2(), Mat_nu()];
  nu_fe[] = InterpolationLinear[SquNorm[$1]]{List[Mat_b2_nu]};
  dnudb2_fe[] = dInterpolationLinear[SquNorm[$1]]{List[Mat_b2_nu]};
  nu[Iron] = nu_fe[$1];
  dhdb_NL[Iron] = 2*dnudb2_fe[$1#1]*SquDyadicProduct[#1];
}
Constraint {
  { Name Dirichlet_a; Case { { Region Dir; Value 0; } } }
}
FunctionSpace {
  { Name Hcurl_a; Type Form1P;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge; Support Vol; Entity NodesOf[All]; }
    }
    Constraint { { NameOfCoef ae; EntityType NodesOf; NameOfConstraint Dirichlet_a; } }
  }
}
Jacobian { { Name JVol; Case { { Region All; Jacobian VolAxiSqu; } } } }
Integration {
  { Name I1; Case { { Type Gauss; Case { { GeoElement Triangle; NumberOfPoints 4; } } } } }
}
Formulation {
  { Name MagSta_a; Type FemEquation;
    Quantity { { Name a; Type Local; NameOfSpace Hcurl_a; } }
    Equation {
      Galerkin { [ nu[{d a}] * Dof{d a}, {d a} ]; In Vol; Jacobian JVol; Integration I1; }
      Galerkin { JacNL[ dhdb_NL[{d a}] * Dof{d a}, {d a} ]; In Vol_NL; Jacobian JVol; Integration I1; }
      Galerkin { [ -nu[] * br[], {d a} ]; In Magnet; Jacobian JVol; Integration I1; }
    }
  }
}
Resolution {
  { Name MagSta;
    System { { Name A; NameOfFormulation MagSta_a; } }
    Operation {
      InitSolution[A];
      IterativeLoop[60, 1e-7, 1] { GenerateJac[A]; SolveJac[A]; }
      SaveSolution[A];
    }
  }
}
PostProcessing {
  { Name MagSta; NameOfFormulation MagSta_a;
    Quantity { { Name b; Value { Local { [ {d a} ]; In Vol; Jacobian JVol; } } } }
  }
}
PostOperation {
  { Name Out; NameOfPostProcessing MagSta;
    Operation { Print[ b, OnElementsOf Vol, Depth 0, Format SimpleTable, File "b.txt" ]; }
  }
}
"""


def write_pro(p, path):
    mt = fm.MAGNET_TYPES[p["magnet_type"]]
    tab = b2_nu_table(p["iron"])
    with open(path, "w") as f:
        f.write(PRO % {
            "iron": PHYS_IRON, "magnet": PHYS_MAGNET, "air": PHYS_AIR, "dir": PHYS_DIR,
            "mur": mt["mur"], "br": p["br"],
            "b2": ", ".join("%.8g" % b2 for b2, _ in tab),
            "nu": ", ".join("%.8g" % nu for _, nu in tab),
        })


def run_getdp(workdir):
    exe = find_getdp()
    if not exe:
        raise RuntimeError("getdp が見つかりません")
    cmd = [exe, "speaker.pro", "-msh", "speaker.msh", "-solve", "MagSta", "-pos", "Out", "-v", "3"]
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=600)
    if r.returncode != 0 or not os.path.exists(os.path.join(workdir, "b.txt")):
        raise RuntimeError("GetDP の実行に失敗しました:\n" + (r.stdout + r.stderr)[-2000:])
    return r.stdout + r.stderr


def read_b_elements(path, tris, nodes):
    """b.txt (SimpleTable, Depth 0) は 1 行 1 要素で「重心 x y z  bx by bz」。重心で要素に対応づける。"""
    def key(x, y):
        return (round(x, 6), round(y, 6))
    cent = {}
    for eid, (_, a, b, c) in tris.items():
        cent[key((nodes[a][0] + nodes[b][0] + nodes[c][0]) / 3, (nodes[a][1] + nodes[b][1] + nodes[c][1]) / 3)] = eid
    bel = {}
    with open(path) as f:
        for line in f:
            w = line.split()
            if len(w) < 6:
                continue
            eid = cent.get(key(float(w[0]), float(w[1])))
            if eid is not None:
                bel[eid] = (float(w[3]), float(w[4]))
    missing = len(tris) - len(bel)
    if missing:
        raise RuntimeError(f"GetDP の出力と要素の対応がとれません ({missing} 要素)")
    return bel


# ---------------------------------------------------------------- 補間器
class Field:
    """要素ごとの B を、領域ごとに節点平均してから三角形内で線形補間する (FEMM の平滑化に相当)"""

    def __init__(self, nodes, tris, bel):
        self.nodes = nodes
        self.tris = tris          # eid -> (reg, n1, n2, n3)
        acc = {}
        for eid, (reg, a, b, c) in tris.items():
            (xa, ya), (xb, yb), (xc, yc) = nodes[a], nodes[b], nodes[c]
            area = abs((xb - xa) * (yc - ya) - (xc - xa) * (yb - ya)) / 2
            br, bz = bel[eid]
            for n in (a, b, c):
                s = acc.setdefault((reg, n), [0.0, 0.0, 0.0])
                s[0] += br * area
                s[1] += bz * area
                s[2] += area
        self.nb = {k: (v[0] / v[2], v[1] / v[2]) for k, v in acc.items()}
        # 検索用の一様バケット (構造の周辺だけ。外側は線形探索でよい)
        xs = [xy[0] for xy in nodes.values()]
        ys = [xy[1] for xy in nodes.values()]
        self.x0, self.y0 = min(xs), min(ys)
        self.h = max(max(xs) - self.x0, max(ys) - self.y0) / 400.0
        self.bk = {}
        for eid, (reg, a, b, c) in tris.items():
            px = [nodes[a][0], nodes[b][0], nodes[c][0]]
            py = [nodes[a][1], nodes[b][1], nodes[c][1]]
            i0, i1 = int((min(px) - self.x0) / self.h), int((max(px) - self.x0) / self.h)
            j0, j1 = int((min(py) - self.y0) / self.h), int((max(py) - self.y0) / self.h)
            if (i1 - i0 + 1) * (j1 - j0 + 1) > 400:
                self.bk.setdefault("big", []).append(eid)
                continue
            for i in range(i0, i1 + 1):
                for j in range(j0, j1 + 1):
                    self.bk.setdefault((i, j), []).append(eid)

    def _locate(self, x, y):
        key = (int((x - self.x0) / self.h), int((y - self.y0) / self.h))
        for eid in self.bk.get(key, []) + self.bk.get("big", []):
            reg, a, b, c = self.tris[eid]
            (xa, ya), (xb, yb), (xc, yc) = self.nodes[a], self.nodes[b], self.nodes[c]
            d = (yb - yc) * (xa - xc) + (xc - xb) * (ya - yc)
            l1 = ((yb - yc) * (x - xc) + (xc - xb) * (y - yc)) / d
            l2 = ((yc - ya) * (x - xc) + (xa - xc) * (y - yc)) / d
            l3 = 1 - l1 - l2
            if min(l1, l2, l3) >= -1e-9:
                return reg, (a, b, c), (l1, l2, l3)
        return None

    def getb(self, r, z):
        hit = self._locate(abs(r), z)
        if hit is None:
            return 0.0, 0.0
        reg, ns, ls = hit
        br = sum(self.nb[(reg, n)][0] * l for n, l in zip(ns, ls))
        bz = sum(self.nb[(reg, n)][1] * l for n, l in zip(ns, ls))
        return br, bz


# ---------------------------------------------------------------- 入口
def analyze(params, keep=None):
    t0 = time.time()
    p = fm.normalize(params)
    errs = fm.validate(p)
    if errs:
        return {"ok": False, "errors": errs, "params": p}
    G = fm.geometry(p)
    wd = keep or tempfile.mkdtemp(prefix="spk_getdp_")
    os.makedirs(wd, exist_ok=True)
    try:
        nodes, tris = build_mesh(p, G, os.path.join(wd, "speaker.msh"))
        t_mesh = time.time() - t0
        write_pro(p, os.path.join(wd, "speaker.pro"))
        log = run_getdp(wd)
        t_solve = time.time() - t0 - t_mesh
        bel = read_b_elements(os.path.join(wd, "b.txt"), tris, nodes)
        field = Field(nodes, tris, bel)
        res = fm.postprocess(p, G, field.getb)
    finally:
        if keep is None:
            shutil.rmtree(wd, ignore_errors=True)
    import re
    m = re.search(r"converged \((\d+) iterations", log)
    iters = int(m.group(1)) if m else None
    res.update({
        "ok": True, "backend": BACKEND_NAME, "params": p,
        "geometry": {k: round(v, 5) for k, v in G.items()},
        "mesh": {"nodes": len(nodes), "elements": len(tris), "t_mesh_s": round(t_mesh, 2),
                 "t_solve_s": round(t_solve, 2), "iterations": iters},
        "elapsed_s": round(time.time() - t0, 2),
    })
    return res


if __name__ == "__main__":
    import json
    import sys
    r = analyze(fm.DEFAULTS, keep=sys.argv[1] if len(sys.argv) > 1 else None)
    json.dump({k: v for k, v in r.items() if k not in ("map", "br_z", "bl_x")}, sys.stdout,
              ensure_ascii=False, indent=1)
