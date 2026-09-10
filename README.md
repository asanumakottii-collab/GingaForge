# GingaForge

天の川専用の原盤生成ツールです。Mellinger全天画像からOHPフィルム用の階調原盤を、
Gaia DR3の観測光量からエッチング用の円形穴原盤を生成します。
恒星原盤は別プロジェクトのOTLが担当します。このフォルダだけで実行でき、OTLのコードや星表への参照はありません。
旧 `galaxy.mode = legacy` の生成処理は廃止しました。

## 実行方法

Python 3.11以上を使用し、このフォルダで次を実行します。

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt

# OHPフィルム: 既存のFITSからPNG・2ページPDF・校正シートを生成
.venv/bin/python generate.py ohp

# Gaiaエッチング: 保存済みキャッシュからSVGと光量レポートを生成
.venv/bin/python generate.py gaia-etch

# 印刷用ベクトルPDF
.venv/bin/python generate.py gaia-etch -PDF

# キャッシュがない場合、または等級・集計ビンを変更した場合
.venv/bin/python generate.py gaia-etch --download-gaia

# 通信も原盤出力もせず、取得クエリを確認
.venv/bin/python generate.py gaia-etch --gaia-query

# 別の設定ファイルを使用
.venv/bin/python generate.py gaia-etch -f gaia-etch.properties

# OHP保存画像・PDFと、Gaia入出力・投影・穴間隔・統合コマンドを検証
.venv/bin/python -m unittest discover -v
```

| 方式 | 設定ファイル | 入力 | 主な出力 |
|---|---|---|---|
| `ohp` | `galaxyconfig.properties`（従来のまま） | `source/mwpan2_RGB_3600.fits` | `output/mellinger-ohp-plates.pdf`、PNG、校正PDF |
| `gaia-etch` | `gaia-etch.properties`（OTLから移行） | `source/gaia/`、または個別星CSV/CSV.gz | `output_inverse/galaxy/gaia-etch-*.svg/pdf`、レポート |

OHPの従来コマンド `python3 build_mellinger_ohp.py --config galaxyconfig.properties` も使用できます。
Gaiaも `python3 galaxy_transformer.py -f gaia-etch.properties` で直接実行できます。
生成方式はコマンドで選択するため、`galaxy.mode` の指定は不要です。
以前のGaia設定を持ち込む場合の `galaxy.mode = gaia-etch` は受け付けます。

両方式のN0/N1/S0/S1の銀経領域を `galaxy_geometry.py` に集約しています。
中心配置では円筒投影の向きと座標が一致します。Gaiaは偏心投影機にも対応し、OHPは中心配置のみです。
OHPの緑チャンネルとGaiaのG帯は異なる観測量のため、それぞれの光量処理・加工条件を使用します。

## コードの構成

| ファイル | 役割 |
|---|---|
| `generate.py` | OHP / Gaiaの共通コマンド |
| `build_mellinger_ohp.py` | FITSの再投影、階調変換、OHP用PNG・PDF |
| `galaxy_transformer.py` | Gaia原盤の生成コマンド |
| `galaxy_geometry.py` | 共通領域定義、Gaiaの円筒投影・逆投影 |
| `gaia_catalog.py` | Gaia DR3取得、キャッシュ検証、CSV入力 |
| `galaxy_etching.py` | 光量の集計、穴の配置、加工寸法・光量誤差の記録 |
| `models.py` / `mathvector.py` | 座標・穴データ、ICRSと銀河座標の変換 |
| `plate_writer.py` | エッチング原盤のSVG・ベクトルPDF |
| `config.py` | Gaia用プロパティファイルの読込 |

## Gaiaエッチング原盤の担当領域

`gaia-etch.properties` の `plate.frame-size` で、正方形原盤の一辺を出力上のmmで指定できます。OHP用の `galaxyconfig.properties` とは独立した設定です。例えば `100` なら100 × 100 mmの枠になります。`0` または省略時は用紙サイズと配置数から自動計算し、標準配置では138.5 × 138.5 mmです。`scale` は枠サイズには掛かりません。SVG / PDFに適用され、Gaiaの穴配置も指定した枠内に制限されます。

`GalaxyTransformer` は星表のICRS赤道座標を銀河座標へ変換し、**銀緯 −20度以上 +20度以下**の星を次の4枚に割り当てます。銀経は0度以上360度未満に正規化し、共有する境界上の星が重複しないよう、各区間の下端を含み上端を含めません。

| 原盤番号 | 銀経の範囲 | 中心銀経 |
|---|---|---|
| N0 | 0度以上60度未満 | 30度 |
| N1 | 60度以上120度未満 | 90度 |
| S0 | 300度以上360度未満 | 330度 |
| S1 | 240度以上300度未満 | 270度 |

銀経120度以上240度未満、および銀緯±20度の外は出力対象外です。N/Sは用紙の上段/下段の識別子で、赤緯・銀緯の正負を意味しません。`plate.column = 1`、`plate.row = 1` の標準設定では、`gaia-etch-0.svg` にN0/S0、`gaia-etch-1.svg` にN1/S1を配置し、計4枚の原盤を2用紙に出力します。PDFも同じ配置です。星がない領域も枠を出力します。

各原盤では中心銀経・銀緯0度への光線を原盤中央に合わせ、円筒投影で位置を求めます。`projector-horizontal` と `projector-vertical` は従来どおり赤道座標での距離で、投影機位置は `(0, ±projector-horizontal, ±projector-vertical)` です。横成分の符号は原盤番号0で正、1で負、縦成分はNで正、Sで負です。

座標変換の回転行列はHipparcosのICRS定義に基づきます（[ERFAの参照実装](https://github.com/liberfa/erfa/blob/master/src/icrs2g.c)）。

### Gaia-DR3からぎんとう原盤を作成

エッチング用にはGaia DR3の光量集計方式だけを使用します。

Gaiaモードは **Gaia DR3のG帯の観測光量だけ**を使います。Hipparcos・Tycho・RC3を重ねず、
V等級への変換や星間減光の除去もしません。暗黒帯の効果を残した、地球から見える恒星光を扱います。
Gaiaの混雑領域での欠測や限界等級より暗い恒星、散光星雲の光は補完しません。
これはG帯に基づく投影原稿で、肉眼の暗所視や写真の完全な再現ではありません。

初回の `--download-gaia` はESA Gaia Archiveに非同期ADQLクエリを送り、4領域それぞれの
`7.5 < G <= 17` の天体を銀経・銀緯0.1度のセルへ集計します。星をランダムに間引かず、
各セルで `SUM(10**(-0.4*(G - reference_magnitude)))` と天体数を取得します。
取得量を抑えるため個々の星をダウンロードせず、セル中心に合計光量を置いて再投影します。
`--gaia-query` で使用する4本のクエリを通信せず確認できます。
取得には時間がかかる場合があり、ジョブURLと状態を表示します。結果の上限超過（OVERFLOW）は拒否します。
途中で失敗した場合、同じコマンドを再実行すると、検証済みの領域ファイルは再利用します。

`galaxy.gaia.input` の既定値は設定ファイルからの相対パス `source/gaia` です。
クエリ・等級・領域・集計幅・チェックサムを保存し、通常実行時は通信せず検証して読み込みます。
入力・出力の相対パスは設定ファイルの所在フォルダが基準です。
Gaiaの等級範囲・取得ビン幅を変えた場合は `--download-gaia` で再取得してください。
個別星のCSV/CSV.gzをこの設定に指定することもできます。列は
`source_id,ra,dec,phot_g_mean_mag` が必須（赤経・赤緯はICRS、度）で、追加列は無視します。
不正値・重複IDはエラーにします。個別CSVは網羅性を確認できないため、その旨をレポートに記録します。
個別CSVの重複チェックはメモリを使用するため、大規模データでは集計キャッシュを推奨します。

| 設定 | 同梱設定の値 | 意味 |
|---|---:|---|
| `galaxy.etch.min-hole-diameter-mm` | 0.02 | 最終原盤での最小穴直径。すべての穴にこの径を使用 |
| `galaxy.etch.min-web-mm` | 0.2 | 穴の縁と縁の間に残す金属の最小幅 |
| `galaxy.etch.cell-size-mm` | 0.4 | 原盤で光量を合計するセル幅。穴径＋金属幅の整数倍に切り下げ |
| `galaxy.etch.flux-gain` | 0.05 | 全領域共通の光量倍率（穴面積の合計に掛ける）。微光星増加による飽和を抑える初期値 |
| `galaxy.etch.seed` | 0 | セル内の穴の選択を再現する乱数種 |
| `galaxy.gaia.bright-limit` | 7.5 | G等級の最輝側境界。この等級自体は含まない |
| `galaxy.gaia.faint-limit` | 17 | G等級の最微側境界。この等級自体を含む |
| `galaxy.gaia.sky-bin-deg` | 0.1 | 取得時の天球上の集計幅。0.05〜1度、60度を等分する値 |
| `galaxy.gaia.reference-magnitude` | 7.5 | `star-radius-7.5` の穴面積に対応させるG等級 |

加工寸法は未実測の仮値です。寸法は最終金属原盤でのmmで指定し、出力時は `scale` を掛けます。
`star-radius-7.5` は従来どおり**縮小前の製図上の半径**です。`scale` を変えるときはこの半径も
同じ倍率にしてください。SVG/PDFは「用紙に合わせる」を使わず、指定倍率で製版します。

光量から必要な穴面積を計算し、原盤上のセルごとに合計します。穴の候補点は原盤全体で共通の格子上に
置き、ピッチを穴径＋金属幅とするため、セル境界をまたいでも金属幅を確保します。
各セルの必要穴数の小数部分は蛇行順に次のセルへ繰り越し、容量制限後の原盤全体の量子化誤差を
最小穴0.5個分以内に抑えます。ゼロ光量セルに穴を足すことはありません。
配置可能な穴数を超えた光は別の場所へ移さず、飽和損失として記録します。
穴は**局所光量の標本で、明るい星も含め一対一の実際の恒星位置を表すものではありません**。
実効的な細かさは取得ビン・原盤セル・投影光学系で制限されます。

Gaiaモードは星座線を加工用原盤へ追加せず、入力待ちなしで生成します。
`output_inverse/galaxy/gaia-etch-0.svg` と `gaia-etch-1.svg`（`-PDF`指定時はPDF）に計4原盤を出力し、
接頭辞は `galaxy.gaia.file.prefix` で変更できます。
`gaia-etch-report.json` に入力の出典・設定・天体数・穴数・目標面積・出力面積・飽和損失・量子化誤差を記録します。
穴面積と投影光量の比例は均一照明・薄板を仮定しており、板厚、光源、エッチング精度、
投影時のぼけによる影響は別途試作・測定してください。最小寸法を満たすことは加工成功の保証ではありません。

GaiaモードのG等級境界は、恒星投影機のV等級境界と同一ではありません。
恒星投影機と同時に使う際は、色による違いも含めて重複・明るさを確認してください。

出典: [Gaia DR3](https://www.cosmos.esa.int/web/gaia/dr3)、
[ESA Archiveのプログラムアクセス](https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access)。
GaiaデータはESA/Gaia/DPACによります。発表等では[Gaiaの謝辞・引用案内](https://www.cosmos.esa.int/web/gaia-users/credits)に従ってください。

同梱の `gaia-etch.properties` は移行前の調整値を引き継いでいます。枠は60 mm四方、
基準穴半径は0.1 mm、`color.invert = no`（白地に黒穴）、出力先は `output_inverse/galaxy/` です。
`color.invert = yes` なら黒地に白穴になります。OHPの138.5 mm枠とは独立しています。
`output/galaxy/` と `output_inverse/galaxy/` に以前のGaia生成物も移しています。
既存ファイルを更新せず試す場合は、`output.directory` を別フォルダへ変更してください。

## Mellinger 全天画像による OHP 原盤

Axel Mellinger の科学用全天パノラマ FITS を銀河座標で再投影した、白黒階調の OHP 用原盤です。観測された天の川の濃淡と暗黒帯を利用しており、手描きの補完・生成画像・RC3 銀河の疑似星化は使用しません。

### 印刷するファイル

- `output/mellinger-ohp-plates.pdf`：A4、2ページ。1ページ目は N0 / S0、2ページ目は N1 / S1。黒い正方形は 138.5 mm 四方。
- `output/ohp-calibration.pdf`：同じ印刷条件で使う21段階の階調見本と、縦横 50 mm の寸法確認線。
- `output/preview.png`：4区画の比較用プレビュー。寸法を合わせた印刷には使用しません。
- `output/{N0,N1,S0,S1}-ohp.png`：区画別の8ビットグレースケール画像。PDFには物理寸法を正確に指定して配置しています。PNGだけを印刷する場合は下記の寸法を指定してください。

| 区画 | 銀経 l | 銀緯 b | 紙面の上 → 下 |
|---|---|---|---|
| N0 | 0–60° | −20–+20° | 60 → 0° |
| N1 | 60–120° | −20–+20° | 120 → 60° |
| S0 | 300–360° | −20–+20° | 360 → 300° |
| S1 | 240–300° | −20–+20° | 300 → 240° |

全区画で紙面の左が銀緯 −20°、右が +20°です。銀経0°と360°は同じ方向です。各画像の外周が指定範囲の境界に対応し、画素は境界の内側の画素中心をサンプリングします。境界の重複画素やフェザー処理は追加していません。

**PDFをA4・倍率100%・「用紙に合わせる」なしで印刷**してください。画像の幅は **36.397023 mm**、高さは **52.359878 mm** です（初期設定）。黒は遮光、白はインクを置かない透明部分、灰色は中間の透過を意図します。印刷時の白インクは不要です。ラベルと位置合わせマークは画像領域の外側にあります。取り付け時は画像領域以外を遮光してください。

紙面座標は既存OTLの計算方向に合わせてあり、追加の鏡像処理はしていません。フィルムの表裏・レンズによる左右反転は、投影機への取り付け方向に依存します。ラベルの銀経・銀緯方向を実機で照合してください。画像の外周を切り落とさず、黒い領域を含めて機構に合わせて切り出します。

### 投影式と光量の扱い

OTL由来の中心配置（ドーム半径6500 mm、水平・垂直偏心0 mm、原盤まで50 mm、倍率1）に対応します。中心配置ではドーム半径は式から消えます。

```
x = p × tan(b)
y = −p × (l − l_center)
```

角度はラジアン、xは紙面右向き、yは紙面下向き、pは原盤までの距離です。これを逆変換し、FITSのWCSで全天画像上の画素位置を求め、双線形補間します。平面への単なる長方形切り抜きではありません。元画像の強度を補間し、面積ヤコビアンを掛ける処理は行いません。レンズ・光源の周辺減光を補償する測定値は未取得です。

RGB FITSの中央のチャンネル（緑成分）を使用します。このG成分はGaiaのG等級ではなく、厳密な測光Vバンドや暗所視感度とも同一ではありません。色は白黒に統一されます。元FITSには強度の物理単位を示すBUNITがありません。

4区画共通の初期表示変換は、元強度 F に対して次の通りです。

```
u = clip((F − black-level) × exposure / (white-level − black-level), 0, 1)
v = asinh(u / asinh-softening) / asinh(1 / asinh-softening)
pixel = round(255 × v^(1 / display-gamma))
```

初期値は黒点0、白点1500、asinh-softening 0.05、exposure 1、gamma 1。白点1500は、4区画の元画像を等角度でサンプルした強度の約99パーセンタイルを参考に設定しました。区画ごとの自動露出・局所的なコントラスト補正は行いません。明るい星の一部は白で飽和します。`tone-mode = linear` にすると asinh の段階を省略できます。

**階調圧縮後の画素値は、空の表面輝度に比例するものではありません。** またプリンターの灰色値と実際のフィルム透過率も一致するとは限りません。初期値は実機調整前の表示設定です。テストシートを同じプリンター・フィルム・インク・設定で印刷し、黒の光漏れと灰色の透過を確認してください。測定した透過率に合わせた逆変換や、投影光学系の補正は本版には含まれません。

元画像の恒星も残しています。別の恒星投影機と同時に使うと、明るい星が重なる可能性があります。星除去は行っていません。白黒原盤のため、長時間露光写真の色そのものは再現しません。

### データと解像度

- 作者：Axel Mellinger
- 論文：[Mellinger (2009), PASP 121, 1180–1187](https://doi.org/10.1086/648480)、[arXiv:0908.4360](https://arxiv.org/abs/0908.4360)
- 原画像：[科学用 RGB FITS（3600 × 1800）](https://galaxy.phy.cmich.edu/~axel/mwpan2/mwpan2_RGB_3600.fits)
- 原画像の保存場所：`source/mwpan2_RGB_3600.fits`
- 元画像のサンプリング：0.1°/画素（6分角）。作者の最高解像度版とは異なります。
- 出力：1200 dpi の補間画像。画素数を増やしても、観測情報の解像度は増えません。
- `output/*-linear-G.npy`：表示用変換前の再投影済み32ビット浮動小数点データ。上下左右はPNGと同じ。
- `output/provenance.json`：原画像のSHA-256、URL、投影式、設定、区画ごとの強度統計。
- `output/source-fits-header.txt`：元のFITSヘッダー。

**利用条件：教育利用についても作者がライセンス案内を設けています。** ダウンロード可能であることは、上映・公開配布の許諾を意味しません。この作業では許諾取得・作者への連絡は行っていません。部の上映や外部配布の前に、[作者の利用条件](https://milkywaysky.com/licenses.html)を確認し、必要な許諾を取得してください。

### 再生成

このフォルダ内で実行します。

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python build_mellinger_ohp.py
.venv/bin/python -m unittest -v test_mellinger_ohp.py
```

元画像がない場合は、上記URLから別途取得して `source/mwpan2_RGB_3600.fits` に置いてください。スクリプトは自動取得しません。再生成では `output/` 内の同名生成物を更新します。

`galaxyconfig.properties` で原盤までの距離・出力倍率・解像度・濃淡を変更できます。原盤の物理寸法を変えたい場合は `plate-distance-mm` を変更します。`print-scale` を変えるとPDF上の画像寸法も変わりますが、黒枠と配置は固定です。指定した銀経区画は`galaxy_geometry.py` の `REGIONS` に記録しています。投影機が偏心する構成には未対応で、設定すると明示的にエラーにします。

入力・出力の相対パスは `--config` で指定した設定ファイルの所在フォルダを基準に解決します。

## コードのライセンス

OTLから移行したコードはGPL-2.0-or-laterです。原著作権表記およびライセンスヘッダーを保持しています。
Copyright (C) 2007, 2012–2014, 2026 東京大学地文研究会天文部。
Mellinger画像およびGaiaデータの出典・利用条件は、各方式の説明を参照してください。
