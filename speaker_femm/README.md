# 内磁型スピーカ磁気回路解析ツール (FEMM 4.2 / Gmsh + GetDP)

カップ形ヨーク + 円板マグネット (+z 着磁) + トッププレートの軸対称静磁場解析を
FEMM または Gmsh + GetDP で実行し、ブラウザで結果を見るローカル Web ツール。

## 使い方 (macOS / Linux: Gmsh + GetDP)

```
pip3 install gmsh
# GetDP: https://getdp.info から macOS 用を入手して展開し、getdp を PATH に置く
#        (または export GETDP=/path/to/getdp)。Linux は apt install getdp でも可
cd speaker_femm
python3 server.py            # FEMM が無ければ自動で GetDP を使う (--backend getdp で明示)
```

macOS でダウンロードした getdp が「開発元を確認できない」で止まる場合は
`xattr -d com.apple.quarantine /path/to/getdp` を実行する。

## 使い方 (Windows: FEMM)

```
pip install pyfemm pywin32
cd speaker_femm
python server.py
```

ブラウザで http://127.0.0.1:8765/ を開き「解析」を押す。
FEMM 4.2 は `C:\femm42` に通常インストールされている前提 (pyfemm は COM 経由で起動)。

解析ソフトが無い環境では `python server.py --backend mock` で UI だけ動かせる
(簡易磁気回路から合成した **参考にならない** 値。画面に MOCK と表示される)。

## ファイル

| ファイル | 内容 |
|---|---|
| `femm_model.py` | `analyze(params)` がモデル作成→FEMM 解析→後処理して JSON 化可能な dict を返す |
| `server.py` | 標準ライブラリのみの HTTP サーバ。FEMM 専用ワーカースレッドで CoInitialize、ロックで排他制御 |
| `static/index.html` | 外部ライブラリなしの UI (入力・断面図 SVG・結果・履歴) |
| `getdp_model.py` | Gmsh でメッシュ、GetDP で解析する版。`analyze(params)` は femm_model と同じ形の結果を返す |
| `bh_data.py` | FEMM 材料ライブラリ (matlib.dat) から写した鉄の B-H カーブ |
| `mock_femm.py` | FEMM 無しで動かすための pyfemm 互換モック |
| `tests/` | モックで形状と後処理を検査 (`test_model.py`)、GetDP 版を解析解・基準値・メッシュ収束で検査 (`test_getdp.py`) |

## モデル

- 単位 mm、軸対称 (`mi_probdef(0,'millimeters','axi',1e-8,0,30)`)
- ヨーク側壁の上端はプレート上面と面一。ヨーク内径 = プレート半径 + ギャップ
- 鉄: FEMM ライブラリの非線形材料 (`1010 Steel` など)。マグネット: Hc = Br/(μ0·μr) の線形材料 (NdFeB μr 1.05 / フェライト μr 1.10)
- 外側境界: 半径 = 構造寸法の 5 倍の半円、A = 0
- ギャップ (プレート側面〜ヨーク内面) をコイル可動範囲の上まで別ブロックにしてメッシュを細かくする

## 評価

- ギャップ中心 (r = プレート半径 + ギャップ/2) の Br(z) を 0.1 mm 間隔で取得
- Bl(x) = 巻数 × 2πr × (コイル巻き幅内の Br 平均)。1 回の解析結果の Br(z) 上でコイル位置をずらして計算
- Xmax: Bl が静止位置の 82% に落ちる変位 (上下の小さい方)
- 鉄 (プレート・側壁・底) の |B|: 断面マップの格子点から平均・95%値・最大。判定は 95% 値で <1.5 T OK / <1.8 T 注意 / それ以上 飽和 (角の特異点を避けるため)
- 質量: 鉄 7.87、NdFeB 7.5、フェライト 4.9 g/cm³

## GetDP 版の定式化と検証

- 軸対称のベクトルポテンシャル A_φ (BF_PerpendicularEdge + VolAxiSqu)、鉄は Newton 法の非線形解析
- 鉄の B-H は FEMM と同じ点列 (`bh_data.py`)。最終点より先は傾き μ0 で延長
- 要素ごとの B を領域ごとに節点平均して線形補間 (FEMM の表示値の平滑化に相当)
- 検証: 空気中の円柱磁石 (鉄なし・μr = 1) の軸上 Bz が解析解と 1% 以内で一致
- メッシュ収束: ギャップのメッシュ 0.1 / 0.05 / 0.025 mm でギャップ B は 0.973 / 0.972 / 0.972 T

## 基準値の結果と未確認事項

基準値 (φ20×5, Br 1.3 T, プレート 3, ギャップ 1, 側壁 2, 底 2.5 mm) の GetDP 結果:

| 項目 | 値 |
|---|---|
| ギャップ B (中心) | 0.973 T |
| ギャップ B (プレート厚平均) | 0.951 T |
| Bl (60 回巻, 巻き高さ 6 mm) | 2.76 T·m |
| Xmax (82%) | ±1.75 mm |
| ヨーク側壁・底 | 95% 値 2.03 T で飽和 |

仕様の目安「約 0.92 T」とは約 5% 差がある。メッシュを細かくしても値は変わらないので、
差はモデルの前提 (プレート半径、側壁上端の高さ、マグネットの μr など) か FEMM 側の条件の違いと考えられる。
ヨークが飽和しているため、ギャップ B は鉄の B-H カーブとヨーク寸法に敏感
(側壁 3 / 底 3.5 mm にすると 1.01 T)。

FEMM での実解析はこの環境 (Linux) では未確認。Windows で `python femm_model.py` を実行すると
基準値の結果 (マップ以外) を表示するので、GetDP 版と比べられる。
