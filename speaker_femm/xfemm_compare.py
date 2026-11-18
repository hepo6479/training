"""xfemm 版と GetDP 版を同じ条件で解いて比べる。

    python xfemm_compare.py ['[{"gap": 0.8}, ...]']

femmcli は PATH か環境変数 FEMMCLI、getdp は PATH か環境変数 GETDP で指定する。
"""

import json
import sys

import getdp_model as gm
import xfemm_model as xm

if __name__ == "__main__":
    cases = json.loads(sys.argv[1]) if len(sys.argv) > 1 else [{}]
    for c in cases:
        x = xm.analyze(c)
        g = gm.analyze(c)
        print(json.dumps(c, ensure_ascii=False))
        for k in ("b_gap_center", "b_gap_avg", "bl0", "xmax"):
            a, b = x["summary"][k], g["summary"][k]
            print(f"  {k:13s} xfemm={a}  getdp={b}  diff={(b - a) / a * 100:+.2f}%")
        for ix, ig in zip(x["iron"], g["iron"]):
            print(f"  {ix['part']}: mean xfemm={ix['b_mean']} getdp={ig['b_mean']}  b95 {ix['b95']} / {ig['b95']}")
