<p align="center">
  <img src="docs/images/banner.webp" alt="Ginga Forge — 星空のフィルムに写した天の川" width="100%">
</p>

<p align="center">
  <strong>Gaiaの観測データから、プラネタリウムの天の川原盤をつくる。</strong><br>
  Generate Milky Way projection plates for planetarium projectors from Gaia DR3.
</p>

<p align="center">
  <a href="pyproject.toml"><img src="https://img.shields.io/badge/Python-3.11%2B-f4ead7?style=flat-square&labelColor=0a0e1a" alt="Python 3.11+"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/Code-GPL--2.0--or--later-f4ead7?style=flat-square&labelColor=0a0e1a" alt="Code: GPL-2.0-or-later"></a>
  <a href="NOTICE"><img src="https://img.shields.io/badge/Plates-CC%20BY--NC%203.0%20IGO-f4ead7?style=flat-square&labelColor=0a0e1a" alt="Plates: CC BY-NC 3.0 IGO"></a>
  <a href="https://www.cosmos.esa.int/web/gaia/dr3"><img src="https://img.shields.io/badge/Data-Gaia%20DR3-f4ead7?style=flat-square&labelColor=0a0e1a" alt="Data: Gaia DR3"></a>
</p>

<p align="center">
  <a href="#クイックスタート">クイックスタート</a> ·
  <a href="#2つの方式">2つの方式</a> ·
  <a href="#仕組み">仕組み</a> ·
  <a href="#ドキュメント">ドキュメント</a> ·
  <a href="#ライセンス">ライセンス</a>
</p>

---

**GingaForge** は、銀河投影機（ぎんとう）用の天の川原盤を生成するPythonライブラリです。
Gaia DR3の星の光を空の区画ごとに合計し、金属エッチングやOHPフィルム印刷でそのまま使えるSVG・PDFにします。

![天の川原盤のプレビュー。左からS1、S0、N0、N1](docs/images/plates.png)

<p align="center"><sub>同梱設定で生成した <code>gaia-ohp</code> の原盤（階調）。左から S1・S0・N0・N1 の順に、銀経240度から120度までの240度分を担当します。</sub></p>

## 特徴

- **観測データそのまま**: 星を間引かず、`7.5 < G <= 17` の全天体の光量を合計します。暗黒帯も残ります。
- **実際の光学系に合わせた形**: ドーム半径・投影機の位置・レンズの焦点距離から、平らな原板上の位置を計算します。
- **2種類の製版**: 金属エッチング用の穴と、OHPフィルム印刷用の階調・点描に対応します。
- **天の川強調**: 背景に相当する明るさを差し引き、ドームで天の川が際立つようにします（OHP用）。
- **1枚で比較**: 列ごとに設定を変えて並べ、1回のテスト印刷で描き方や明るさを比べられます。

## クイックスタート

Python 3.11以上を使います。このフォルダだけで実行できます。

```sh
# 1. インストール
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'

# 2. 原盤の生成（初回は Gaia の光量も取得する）
.venv/bin/gingaforge gaia-etch --download-gaia
.venv/bin/gingaforge gaia-ohp
```

生成物は `output/<設定ファイル名>/` に保存されます。

> [!IMPORTANT]
> PDFは **倍率100%・「用紙に合わせる」なし** で印刷してください。原盤の寸法が変わると、ドームに正しく投影できません。

## 2つの方式

| | `gaia-etch` | `gaia-ohp` |
|---|---|---|
| **用途** | 金属エッチング | OHPフィルム印刷 |
| **描き方** | 円形の穴 | 階調・点描の画像（穴も選択可） |
| **天の川強調**（同梱設定） | なし | あり |
| **設定ファイル** | [`gaia-etch.properties`](gaia-etch.properties) | [`gaia-ohp.properties`](gaia-ohp.properties) |
| **出力形式** | SVG・PDF | SVG・PDF |

`gaia-ohp` の描き方は `galaxy.ohp.render` で切り替えます。

![同じ領域を gray、dither、holes で描いた比較](docs/images/render-modes.png)

<p align="center"><sub>左から <code>gray</code>（階調）・<code>dither</code>（点描）・<code>holes</code>（穴）。<code>holes</code> は同梱設定で光量倍率を下げた比較用の列です。</sub></p>

| `render` | 描き方 | 向いている用途 |
|---|---|---|
| `gray` | 開口率の8ビット階調画像 | OHP印刷（インクジェット向き） |
| `dither` | 誤差拡散した点の画像 | OHP印刷（レーザー向き） |
| `holes` | 円形の穴 | 金属エッチング |

## 仕組み

```mermaid
flowchart LR
    A["Gaia DR3<br/>G帯の光量"] --> B["0.1度ごとに<br/>光量を集計"]
    B --> C["平らな原板へ<br/>投影"]
    C -->|gaia-etch| D["円形の穴"]
    C -->|gaia-ohp| E["天の川強調"]
    E --> F["階調・点描"]
    D --> G["SVG・PDF"]
    F --> G
```

銀緯±20度の範囲を、4枚の円形原盤に分けて担当します。銀経120度以上240度未満は対象外です。

| 原盤 | 銀経の範囲 | 原盤 | 銀経の範囲 |
|---|---|---|---|
| S1 | 240度以上300度未満 | N0 | 0度以上60度未満 |
| S0 | 300度以上360度未満 | N1 | 60度以上120度未満 |

作る原盤は設定ファイルの `plate.select` で選びます（同梱の設定は `N0, N1`）。

恒星原盤は別プロジェクトの[OTL](https://github.com/asanumakottii-collab/Orb-Transform-Library)が担当します。

## 使い方

寸法・明るさ・描き方などの設定は、同名の `.properties` ファイルで変えます。

| オプション | 内容 |
|---|---|
| `-f <設定ファイル>` | 別の設定ファイルを使う |
| `-svg` / `-pdf` | SVGのみ / PDFのみ出力 |
| `-mirrored` | SVG・PDFに加えて、裏面印刷用の左右反転PDFも出力 |
| `--download-gaia` | Gaiaの光量を取得し直してから生成 |

テストは次のコマンドで実行します。

```sh
.venv/bin/python -m unittest discover -v
```

## ドキュメント

- [全体の構成](docs/overview.md) — 担当領域、出力先、コードの構成
- [Gaia方式](docs/gaia.md) — 投影・穴の配置・描き方・天の川強調の仕組みと、すべての設定

## ライセンス

| 対象 | ライセンス |
|---|---|
| コード | [GPL-2.0-or-later](LICENSE) |
| Gaiaデータと、そこから生成した原盤（`output/`）・プレビュー画像 | [CC BY-NC 3.0 IGO](https://www.cosmos.esa.int/web/gaia-users/license)（クレジット: `ESA/Gaia/DPAC`） |

コードの一部は[OTL](https://github.com/asanumakottii-collab/Orb-Transform-Library)から移行して変更したものです。条件の詳細は[NOTICE](NOTICE)を参照してください。

Copyright (C) 2007, 2012–2014, 2026 東京大学地文研究会天文部
