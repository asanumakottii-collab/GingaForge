# GingaForge

銀河投影機（ぎんとう）用の天の川原盤を生成するライブラリです。

- **gaia-etch**: Gaia DR3の光量 → 金属エッチング用の穴原盤（SVG・PDF）
- **gaia-ohp**: Gaia DR3の光量 → OHPフィルム印刷用の階調・点描原盤（SVG・PDF、天の川強調あり）

恒星原盤は別プロジェクトの[OTL](https://github.com/asanumakottii-collab/Orb-Transform-Library)が担当します。このフォルダだけで実行できます。

## 使い方

Python 3.11以上を使います。

```sh
# 初回のみ
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'

# 原盤の生成
.venv/bin/gingaforge gaia-etch --download-gaia   # 初回は Gaia の光量を取得する
.venv/bin/gingaforge gaia-ohp

# テスト
.venv/bin/python -m unittest discover -v
```

設定は同名の `.properties` ファイル（`gaia-etch.properties` など）で変えます。
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

## ライセンス

コードはGPL-2.0-or-laterです（[LICENSE](LICENSE)）。一部は[OTL](https://github.com/asanumakottii-collab/Orb-Transform-Library)から移行して変更したものです。
Copyright (C) 2007, 2012–2014, 2026 東京大学地文研究会天文部。

Gaiaデータと、そこから生成した `output/` の原盤にはGPLではなく、CC BY-NC 3.0 IGOが適用されます。
クレジットは `ESA/Gaia/DPAC` です。条件は[NOTICE](NOTICE)を参照してください。
