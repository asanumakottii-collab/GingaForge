# GingaForge

銀河投影機（ぎんとう）用の天の川原盤を生成するライブラリです。

- **mellinger-ohp**: Mellinger全天画像 → OHPフィルム用の階調原盤（PDF）
- **gaia-etch**: Gaia DR3の光量 → 金属エッチング用の穴原盤（SVG・PDF）
- **gaia-ohp**: Gaia DR3の光量 → OHPフィルム印刷用の階調・点描原盤（SVG・PDF、天の川強調あり）

恒星原盤は別プロジェクトのOTLが担当します。このフォルダだけで実行できます。

## 使い方

Python 3.11以上を使います。

```sh
# 初回のみ
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/gingaforge install-starnet      # mellinger-ohp で使う StarNet2 を導入（利用条件への同意あり）

# 原盤の生成
.venv/bin/gingaforge mellinger-ohp
.venv/bin/gingaforge gaia-etch --download-gaia   # 初回は Gaia の光量を取得する
.venv/bin/gingaforge gaia-ohp

# テスト
.venv/bin/python -m unittest discover -v
```

設定は同名の `.properties` ファイル（`mellinger-ohp.properties` など）で変えます。
生成物は `output/<設定ファイル名>/` に保存します。

| オプション | 対象 | 内容 |
|---|---|---|
| `-f <設定ファイル>` | すべて | 別の設定ファイルを使う |
| `-svg` / `-pdf` | gaia-etch、gaia-ohp | SVGのみ / PDFのみ出力 |
| `-mirrored` | gaia-etch、gaia-ohp | 裏面印刷用の左右反転PDFのみ出力 |
| `--download-gaia` | gaia-etch、gaia-ohp | Gaiaの光量を取得し直してから生成 |

PDFは **倍率100%・「用紙に合わせる」なし** で印刷してください。

## ドキュメント

- [docs/overview.md](docs/overview.md) — 担当領域、出力先、コードの構成
- [docs/gaia.md](docs/gaia.md) — Gaia方式（gaia-etch・gaia-ohp）の仕組みと設定
- [docs/mellinger-ohp.md](docs/mellinger-ohp.md) — Mellinger方式の仕組みと設定、StarNet2の導入

## ライセンス

OTLから移行したコードはGPL-2.0-or-laterです。
Copyright (C) 2007, 2012–2014, 2026 東京大学地文研究会天文部。

Mellinger画像は、上映・公開配布の前に[作者の利用条件](https://milkywaysky.com/licenses.html)を確認し、必要な許諾を取得してください。
Gaiaデータ・Hipparcos星表・StarNet2の出典と利用条件は、各方式のドキュメントを参照してください。
