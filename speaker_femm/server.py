"""ローカル Web サーバ (Python 標準ライブラリのみ)

    python server.py            # FEMM (pyfemm) で解析。http://127.0.0.1:8765/
    python server.py --mock     # FEMM 無しでモック解析 (UI 確認用)

FEMM は COM (ActiveX) で動くため、FEMM を触るのは専用ワーカースレッド 1 本だけにし、
そのスレッドで CoInitialize してから openfemm する (COM のアパートメントはスレッド単位)。
解析要求はロックで排他制御し、実行中に来た要求は 409 で断る。
"""

import argparse
import json
import os
import queue
import sys
import threading
import traceback
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

import femm_model

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")


class FemmWorker:
    def __init__(self, mock):
        self.mock = mock
        self.lock = threading.Lock()   # 排他制御: 同時に 1 解析だけ
        self.jobs = queue.Queue()
        self.femm = None
        self.backend_error = None
        t = threading.Thread(target=self._run, name="femm-worker", daemon=True)
        t.start()

    # ---- ワーカースレッド
    def _com_init(self):
        try:
            import pythoncom  # pywin32 (pyfemm の依存)
            pythoncom.CoInitialize()
            return pythoncom
        except ImportError:
            return None

    def _open(self):
        if self.mock:
            import mock_femm
            self.femm = mock_femm
            return
        import femm  # pyfemm
        femm.openfemm(1)  # 1 = ウィンドウを隠す
        self.femm = femm

    def _close(self):
        if self.femm is not None and not self.mock:
            try:
                self.femm.closefemm()
            except Exception:
                pass
        self.femm = None

    def _run(self):
        com = self._com_init()
        try:
            while True:
                params, box, done = self.jobs.get()
                try:
                    if self.femm is None:
                        self._open()
                    box["result"] = femm_model.analyze(params, femm=self.femm, keep_open=True)
                except Exception as e:  # FEMM が落ちた等: 次回開き直す
                    traceback.print_exc()
                    self._close()
                    box["result"] = {"ok": False, "errors": [f"解析エラー: {type(e).__name__}: {e}"]}
                finally:
                    done.set()
        finally:
            self._close()
            if com:
                com.CoUninitialize()

    # ---- HTTP スレッドから呼ぶ
    def analyze(self, params, timeout=600):
        if not self.lock.acquire(blocking=False):
            return None  # 実行中
        try:
            box, done = {}, threading.Event()
            self.jobs.put((params, box, done))
            if not done.wait(timeout):
                return {"ok": False, "errors": ["解析がタイムアウトしました"]}
            return box["result"]
        finally:
            self.lock.release()

    @property
    def busy(self):
        return self.lock.locked()


class Handler(SimpleHTTPRequestHandler):
    worker = None

    def __init__(self, *a, **kw):
        super().__init__(*a, directory=STATIC, **kw)

    def log_message(self, fmt, *args):
        if "/api/" in (self.path or ""):
            sys.stderr.write("%s %s\n" % (self.address_string(), fmt % args))

    def _json(self, code, obj):
        body = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path == "/api/config":
            return self._json(200, {
                "backend": "mock" if self.worker.mock else "femm",
                "busy": self.worker.busy,
                "defaults": femm_model.DEFAULTS,
                "iron_materials": femm_model.IRON_MATERIALS,
                "magnet_types": {k: {"br_range": v["br_range"]} for k, v in femm_model.MAGNET_TYPES.items()},
                "xmax_ratio": femm_model.XMAX_RATIO,
                "sat_warn": femm_model.SAT_WARN, "sat_limit": femm_model.SAT_LIMIT,
            })
        return super().do_GET()

    def do_POST(self):
        if self.path != "/api/analyze":
            return self._json(404, {"ok": False, "errors": ["not found"]})
        try:
            n = int(self.headers.get("Content-Length", 0))
            params = json.loads(self.rfile.read(n) or b"{}")
            if not isinstance(params, dict):
                raise ValueError("JSON object が必要です")
            p = femm_model.normalize(params)
        except Exception as e:
            return self._json(400, {"ok": False, "errors": [f"入力が不正です: {e}"]})
        errs = femm_model.validate(p)
        if errs:
            return self._json(400, {"ok": False, "errors": errs})
        res = self.worker.analyze(p)
        if res is None:
            return self._json(409, {"ok": False, "errors": ["別の解析を実行中です。終わってから再実行してください"]})
        return self._json(200 if res.get("ok") else 500, res)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--mock", action="store_true", help="FEMM を使わずモック解析で動かす")
    a = ap.parse_args()
    mock = a.mock
    if not mock:
        try:
            import femm  # noqa: F401
        except ImportError:
            print("pyfemm が見つかりません (pip install pyfemm)。--mock で UI だけ試せます。")
            sys.exit(1)
    Handler.worker = FemmWorker(mock)
    srv = ThreadingHTTPServer((a.host, a.port), Handler)
    print(f"{'MOCK' if mock else 'FEMM'} モードで起動: http://{a.host}:{a.port}/  (Ctrl+C で終了)")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
