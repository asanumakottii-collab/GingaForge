# 全体の構成

## 方式とコマンド

| コマンド | 設定ファイル | 入力 | 出力（同梱設定） |
|---|---|---|---|
| `mellinger-ohp` | `mellinger-ohp.properties` | `source/mwpan2_RGB_3600.fits`、`source/hip_main.dat.gz`、StarNet2 | `output/mellinger-ohp/mellinger-ohp-plates.pdf` |
| `gaia-etch` | `gaia-etch.properties` | `source/gaia/` | `output/gaia-etch/inverted/gaia-etch-*.svg/pdf` |
| `gaia-ohp` | `gaia-ohp.properties` | `source/gaia/` | `output/gaia-ohp/inverted/gaia-ohp-*.svg/pdf` |

`gingaforge gaia-ohp` は `gingaforge gaia-etch -f gaia-ohp.properties` と同じです。
各コマンドのオプションは `gingaforge <コマンド> --help` で表示します。
`gingaforge` コマンドは `python -m gingaforge` でも実行できます。

ライブラリとしては次を使います。

- `from gingaforge.ohp import build` — Mellinger方式のOHP原盤PDF
- `from gingaforge.gaia_etch import generate` — Gaia方式の原盤SVG・PDF

## 担当領域

両方式とも、銀緯 −20度以上 +20度以下を次の4枚の原盤に分けます。
領域の定義は `src/gingaforge/galaxy_geometry.py` の `REGIONS` にまとめています。

| 原盤 | 銀経の範囲 | 中心銀経 |
|---|---|---|
| N0 | 0度以上60度未満 | 30度 |
| N1 | 60度以上120度未満 | 90度 |
| S0 | 300度以上360度未満 | 330度 |
| S1 | 240度以上300度未満 | 270度 |

- 銀経は0度以上360度未満に正規化し、各区間は下端を含み上端を含みません（境界上の星が重複しないため）。
- 銀経120度以上240度未満、および銀緯±20度の外は対象外です。
- N/Sは用紙の上段/下段の識別子で、赤緯・銀緯の正負を意味しません。
- 中心配置では、両方式の投影（レンズ中心を通る光線が平らな原板に当たる位置）の向きと座標が一致します。Gaia方式は偏心投影機にも対応し、Mellinger方式は中心配置のみです。

## 出力先

出力先は設定では指定せず、**設定ファイルと同じフォルダの `output/<設定ファイル名>/`** に保存します。
ファイル名の接頭辞も設定ファイル名です。Gaia方式はさらに `color.invert` で `normal/`（白地に黒穴）と `inverted/`（黒地に白穴）に分けます。

```
output/
  mellinger-ohp/            mellinger-ohp.properties の出力
  gaia-etch/inverted/       gaia-etch.properties の出力（color.invert = yes）
  gaia-etch/normal/         gaia-etch.properties を color.invert = no で生成したもの
  gaia-ohp/inverted/        gaia-ohp.properties の出力（color.invert = yes）
```

設定ごとにフォルダが分かれるため、ある設定の生成が別の設定の原盤を上書きすることはありません。
既存の原盤を残したまま試す場合は、設定ファイルを同じフォルダに別名でコピーし（例: `gaia-test.properties`）、
`-f gaia-test.properties` で指定してください。`output/gaia-test/` に出力されます。

入力ファイルの相対パスも、設定ファイルの所在フォルダが基準です。

設定ファイルの文法はJavaのプロパティファイルに準じます（`キー = 値`、`#` で始まる行はコメント）。

## コードの構成

```
pyproject.toml             パッケージ定義（gingaforge コマンド）
mellinger-ohp.properties   Mellinger方式の設定（OHPフィルム印刷）
gaia-etch.properties       Gaia方式の設定（金属エッチング）
gaia-ohp.properties        Gaia方式の設定（OHPフィルム印刷、天の川強調あり）
src/gingaforge/            ライブラリ本体
source/                    入力データ（FITS、Hipparcos星表、Gaiaキャッシュ）
tools/starnet2/            install-starnet で導入した StarNet2（gitの対象外、再配布しない）
output/                    生成物
tests/                     テスト
docs/                      ドキュメント
```

| モジュール（`src/gingaforge/`） | 役割 |
|---|---|
| `cli.py` | 共通コマンド `gingaforge` |
| `ohp.py` | Mellinger方式: 色の混合、星の除去（StarNet2・Hipparcos）、再投影、階調変換、OHP原盤PDF |
| `starnet_setup.py` | StarNet2 の検索、公式配布元からの導入、利用条件への同意の確認 |
| `gaia_etch.py` | Gaia原盤の生成 |
| `galaxy_geometry.py` | 共通の領域定義、レンズの公式、Gaiaの投影・逆投影 |
| `gaia_catalog.py` | Gaia DR3光量マップの取得・キャッシュ検証 |
| `galaxy_etching.py` | 光量の集計、穴の配置、開口率の画像（階調・点描）、天の川強調 |
| `models.py` / `mathvector.py` | 座標データ、ICRSと銀河座標の変換 |
| `plate_writer.py` | 原盤のSVG・PDF（穴はベクトル、画像は画素） |
| `config.py` | プロパティファイルの読込、リポジトリのルート |

## テスト

```sh
.venv/bin/python -m unittest discover -v
```

StarNet2 がない、または利用条件に同意していない環境では、StarNet2 を実際に動かすテストは飛ばされます。
StarNet2 の導入処理のテストは、代わりの小さなスクリプトと偽のフィードで行い、通信しません。
