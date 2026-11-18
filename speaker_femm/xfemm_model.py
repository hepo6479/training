"""xfemm 版の解析部 (FEMM のソルバの移植版。macOS / Linux)

femm_model.build_and_solve が pyfemm に出す呼び出しを Lua スクリプトに書き出し、
xfemm の femmcli で実行する。B を読む点は femm_model.postprocess と同じ点を前もって集め、
Lua 側でまとめてファイルに書き出す。analyze(params) は femm_model.analyze と同じ形の dict を返す。

    femmcli: https://github.com/crobarcro/xfemm をビルド (cfemm で cmake → make femmcli-bin)
             PATH に置くか、環境変数 FEMMCLI に実行ファイルのパスを入れる

鉄は matlib.dat を読まず、bh_data.py (FEMM ライブラリと同じ B-H 点列) から定義するので、
xfemm の材料ライブラリの置き場所は気にしなくてよい。
"""

import os
import re
import shutil
import subprocess
import tempfile
import time

import femm_model as fm
from bh_data import BH_CURVES

BACKEND_NAME = "xfemm"


def find_femmcli():
    exe = os.environ.get("FEMMCLI") or shutil.which("femmcli")
    if exe and os.path.exists(exe):
        return exe
    for c in ("/usr/local/bin/femmcli", "/opt/homebrew/bin/femmcli",
              os.path.expanduser("~/xfemm/cfemm/bin/femmcli")):
        if os.path.exists(c):
            return c
    return None


def available():
    """使えなければ理由の文字列、使えれば None"""
    if not find_femmcli():
        return "femmcli (xfemm) が見つかりません (PATH か環境変数 FEMMCLI で指定)"
    return None


def _lua(v):
    if isinstance(v, str):
        return '"' + v.replace("\\", "/").replace('"', '\\"') + '"'
    return repr(float(v)) if isinstance(v, float) else str(v)


class LuaFemm:
    """pyfemm の代わりに Lua の行を記録する"""

    def __init__(self):
        self.lines = []

    def _call(self, name, *a):
        self.lines.append(f"{name}({', '.join(_lua(x) for x in a)})")

    def mi_getmaterial(self, name):
        if name == "Air":
            self._call("mi_addmaterial", "Air", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
            return
        # FEMM ライブラリの鉄と同じ: 非線形 B-H、導電率は静磁場では無関係
        self._call("mi_addmaterial", name, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0)
        for b, h in BH_CURVES[name]:
            self._call("mi_addbhpoint", name, float(b), float(h))

    def __getattr__(self, name):
        if name.startswith("mo_") or name in ("openfemm", "closefemm", "mi_close"):
            return lambda *a: None
        return lambda *a: self._call(name, *a)


def collect_points(p, G):
    """postprocess が B を読む点を、値を返さずに集める (点の並びは値に依らない)"""
    pts = []

    def getb(r, z):
        pts.append((r, z))
        return 0.0, 0.0
    fm.postprocess(p, G, getb)
    return pts


def solve(p, G, wd, femmcli=None):
    """解いて (点ごとの (Br, Bz) のリスト, 情報 dict) を返す"""
    femmcli = femmcli or find_femmcli()
    if not femmcli:
        raise RuntimeError("femmcli が見つかりません")
    pts = collect_points(p, G)
    lua = LuaFemm()
    fm.build_and_solve(lua, p, G, wd)          # mi_saveas, mi_analyze, mi_loadsolution まで
    out = os.path.join(wd, "b.txt").replace("\\", "/")
    body = "\n".join(lua.lines) + "\n"
    body += f'fh = openfile("{out}", "w")\n'
    body += 'write(fh, "NELEM ", mo_numelements(), "\\n")\n'
    body += "pts = {" + ",".join(f"{{{r!r},{z!r}}}" for r, z in pts) + "}\n"
    body += ("for i = 1, getn(pts) do\n"
             "  a, b1, b2 = mo_getpointvalues(pts[i][1], pts[i][2])\n"  # xfemm には mo_getb が無い
             "  if b1 == nil then b1 = 0; b2 = 0 end\n"
             '  write(fh, b1, " ", b2, "\\n")\n'
             "end\n"
             "closefile(fh)\n")
    script = os.path.join(wd, "model.lua")
    with open(script, "w") as f:
        f.write(body)
    r = subprocess.run([femmcli, f"--lua-script={script}"], cwd=wd, capture_output=True,
                       text=True, timeout=900)
    log = r.stdout + r.stderr
    if not os.path.exists(out):
        raise RuntimeError("xfemm の実行に失敗しました:\n" + log[-2000:])
    nelem, vals = None, []
    for line in open(out):
        w = line.split()
        if not w:
            continue
        if w[0] == "NELEM":
            nelem = int(float(w[1]))
        else:
            vals.append((float(w[0]), float(w[1])))
    if len(vals) != len(pts):
        raise RuntimeError(f"xfemm の出力点数が合いません ({len(vals)} / {len(pts)})\n" + log[-1000:])
    it = re.findall(r"Newton Iteration\((\d+)\)", log)
    return vals, {"elements": nelem, "iterations": int(it[-1]) if it else None}


def analyze(params, femmcli=None):
    t0 = time.time()
    p = fm.normalize(params)
    errs = fm.validate(p)
    if errs:
        return {"ok": False, "errors": errs, "params": p}
    G = fm.geometry(p)
    wd = tempfile.mkdtemp(prefix="spk_xfemm_")
    try:
        vals, info = solve(p, G, wd, femmcli)
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    it = iter(vals)
    res = fm.postprocess(p, G, lambda r, z: next(it))
    res.update({
        "ok": True, "backend": BACKEND_NAME, "params": p,
        "geometry": {k: round(v, 5) for k, v in G.items()},
        "mesh": info,
        "elapsed_s": round(time.time() - t0, 2),
    })
    return res


if __name__ == "__main__":
    import json
    import sys
    r = analyze(fm.DEFAULTS)
    json.dump({k: v for k, v in r.items() if k not in ("map", "br_z", "bl_x")}, sys.stdout,
              ensure_ascii=False, indent=1)
