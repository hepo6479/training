"""xfemm (FEMM のソルバの Linux/macOS 移植版) で femm_model と同じモデルを解き、GetDP 版と比べる。

femm_model.build_and_solve が pyfemm に出す呼び出しを Lua スクリプトに書き出し、
xfemm の femmcli で実行する。B を読む点も femm_model.postprocess と同じ点を使う。

    python xfemm_compare.py /path/to/femmcli [matlib.dat のあるディレクトリ]
"""

import json
import os
import subprocess
import sys
import tempfile

import femm_model as fm


def _lua(v):
    if isinstance(v, str):
        return '"' + v.replace("\\", "/") + '"'
    return repr(float(v)) if isinstance(v, float) else str(v)


class LuaFemm:
    """pyfemm の代わりに Lua の行を記録する"""
    BACKEND_NAME = "xfemm"

    def __init__(self):
        self.lines = []

    def __getattr__(self, name):
        if name.startswith("mo_") or name in ("openfemm", "closefemm", "mi_close"):
            return lambda *a: None
        return lambda *a: self.lines.append(f"{name}({', '.join(_lua(x) for x in a)})")


def collect_points(p, G):
    pts = []

    def getb(r, z):
        pts.append((r, z))
        return 0.0, 0.0
    fm.postprocess(p, G, getb)
    return pts


def run(femmcli, params, base_dir=None):
    p = fm.normalize(params)
    G = fm.geometry(p)
    pts = collect_points(p, G)
    wd = tempfile.mkdtemp(prefix="xfemm_")
    lua = LuaFemm()
    fm.build_and_solve(lua, p, G, wd)          # mi_saveas, mi_analyze, mi_loadsolution まで
    out = os.path.join(wd, "b.txt")
    body = "\n".join(lua.lines) + "\n"
    body += f'fh = openfile("{out}", "w")\n'
    body += "pts = {" + ",".join(f"{{{r!r},{z!r}}}" for r, z in pts) + "}\n"
    body += ("for i = 1, getn(pts) do\n"
             "  a, b1, b2 = mo_getpointvalues(pts[i][1], pts[i][2])\n"  # xfemm には mo_getb が無い
             '  if b1 == nil then b1 = 0; b2 = 0 end\n'
             '  write(fh, b1, " ", b2, "\\n")\n'
             "end\n"
             "closefile(fh)\n")
    script = os.path.join(wd, "model.lua")
    with open(script, "w") as f:
        f.write(body)
    cmd = [femmcli, "-q", f"--lua-script={script}"]
    if base_dir:
        cmd.append(f"--lua-base-dir={base_dir}")
    r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True, timeout=900)
    if not os.path.exists(out):
        raise RuntimeError(r.stdout[-2000:] + r.stderr[-2000:])
    vals = [tuple(map(float, l.split())) for l in open(out) if l.strip()]
    if len(vals) != len(pts):
        raise RuntimeError(f"点数が合いません {len(vals)} / {len(pts)}")
    it = iter(vals)
    res = fm.postprocess(p, G, lambda r, z: next(it))
    res["workdir"] = wd
    return res


if __name__ == "__main__":
    femmcli = sys.argv[1]
    base = sys.argv[2] if len(sys.argv) > 2 else None
    import getdp_model as gm
    cases = json.loads(sys.argv[3]) if len(sys.argv) > 3 else [{}]
    for c in cases:
        x = run(femmcli, c, base)
        g = gm.analyze(c)
        print(json.dumps(c, ensure_ascii=False))
        for k in ("b_gap_center", "b_gap_avg", "bl0", "xmax"):
            a, b = x["summary"][k], g["summary"][k]
            print(f"  {k:13s} xfemm={a}  getdp={b}  diff={(b - a) / a * 100:+.2f}%")
        for ix, ig in zip(x["iron"], g["iron"]):
            print(f"  {ix['part']}: mean xfemm={ix['b_mean']} getdp={ig['b_mean']}  b95 {ix['b95']} / {ig['b95']}")
