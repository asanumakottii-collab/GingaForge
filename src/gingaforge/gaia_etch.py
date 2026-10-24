# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""Gaia DR3のG帯光量からエッチング用の天の川原盤を生成する。既定でSVGとPDFの両方を出力する。

``column.count`` を2以上にすると、列ごとに設定を変えた原盤一式を用紙に左から並べる。
"""

import argparse
import math
from pathlib import Path
import sys

from .config import PROJECT_ROOT, load_properties, output_directory, positive_number
from .gaia_catalog import GaiaSettings, download_gaia, iter_flux_bins
from .galaxy_etching import EtchingBuilder, EtchingSettings
from .galaxy_geometry import GalaxyTransformer
from .plate_writer import PlateWriterPDF, PlateWriterSVG


def _init_galaxy_transformer(props):
    scale = positive_number(props, "scale", 1.)
    radius = positive_number(props, "star-radius-7.5", .125)
    sphere = positive_number(props, "dome-radius", 5000.)
    plate = positive_number(props, "plate-distance", 50.)
    horizontal = float(props.get("projector-horizontal", 300.))
    vertical = float(props.get("projector-vertical", 200.))
    if (not math.isfinite(horizontal) or not math.isfinite(vertical)
            or math.hypot(horizontal, vertical) >= sphere):
        raise ValueError("投影機位置は有限値で、ドームの内側にしてください。")
    return GalaxyTransformer(radius, sphere * scale, horizontal * scale, vertical * scale, plate * scale)


def _init_plate_writer(props, config_path, pdf=False, column_headings=((),), mirrored=False):
    column = int(props.get("plate.column", 1))
    row = int(props.get("plate.row", 1))
    # plate.frame-size は正方形枠の一辺、plate.frame-diameter は円形枠の直径。0 または省略時は自動計算。
    if "plate.frame-size" in props and "plate.frame-diameter" in props:
        raise ValueError("plate.frame-size（正方形）と plate.frame-diameter（円形）は、どちらか一方だけ書いてください。")
    circular = "plate.frame-diameter" in props
    frame_size = positive_number(props, "plate.frame-diameter" if circular else "plate.frame-size", 0.,
                                 allow_zero=True)
    invert_color = props.get("color.invert", "no").strip().lower() in {"y", "yes", "true"}
    # output/<設定ファイル名>/ の下で、白地に黒穴（normal）と黒地に白穴（inverted）を分ける。
    output_dir = output_directory(config_path, "inverted" if invert_color else "normal")
    args = (column, row, frame_size / 2, f"{Path(config_path).stem}-", invert_color, output_dir, column_headings)
    if mirrored:
        return PlateWriterPDF(*args, circular=circular, mirrored=True)
    return PlateWriterPDF(*args, circular=circular) if pdf else PlateWriterSVG(*args, circular=circular)


def _column_properties(props):
    """列ごとの ``(列で変えた設定, 列の設定)`` を左の列から返す。

    ``column.<列番号>.<キー> = 値`` は、その列だけ ``<キー>`` の値を変える。用紙の配置と白黒は全列共通。
    """
    try:
        count = int(props.get("column.count", 1))
    except ValueError as exc:
        raise ValueError("column.count: 1以上の整数を指定してください。") from exc
    if count < 1:
        raise ValueError("column.count: 1以上の整数を指定してください。")
    changes = [{} for _ in range(count)]
    for key, value in props.items():
        if not key.startswith("column.") or key == "column.count":
            continue
        number, _, name = key.removeprefix("column.").partition(".")
        if not number.isdecimal() or not name:
            raise ValueError(f"{key}: column.<列番号>.<キー> の形で書いてください。")
        if not 1 <= int(number) <= count:
            raise ValueError(f"{key}: 列番号は1〜{count}（column.count）にしてください。")
        if name.startswith(("plate.", "column.")) or name == "color.invert":
            raise ValueError(f"{key}: 用紙の配置・白黒（plate.*、color.invert）は全列共通です。")
        # 打ち間違えたキーが、列の見出しに書かれるだけで無視されないように。
        if name not in props:
            raise ValueError(f"{key}: 列で変える {name} は、列の指定なしでもこのファイルに書いてください。")
        if name in changes[int(number) - 1]:
            raise ValueError(f"{key}: 列{int(number)}の {name} が重複しています。")
        changes[int(number) - 1][name] = value
    return [(change, {**props, **change}) for change in changes]


def _column_headings(changes):
    """列の見出し。1列で設定を変えていなければ見出しを付けない（従来の原盤と同じ）。"""
    if len(changes) == 1 and not changes[0]:
        return [()]
    return [(f"#{number}", *([f"{key} = {value}" for key, value in change.items()] or ["(no changes)"]))
            for number, change in enumerate(changes, 1)]


def _build_column(gaia, etching, transformer, frame_radius, circular=False):
    """1列分の原盤一式の穴を配置し、``(穴, 穴の半径)`` を返します。"""
    regions, limit = GalaxyTransformer.LONGITUDE_REGIONS, GalaxyTransformer.LATITUDE_LIMIT
    print(f"円筒投影: ドーム半径 {transformer.sphere / etching.scale:g} mm / "
          f"原盤距離 {transformer.plate / etching.scale:g} mm / 倍率 {etching.scale:g}")
    builder = EtchingBuilder(transformer, etching, frame_radius, circular)
    print(f"Gaia G={gaia.bright:g}より暗く{gaia.faint:g}以下の光量を集計します。", flush=True)
    for sky_bin in iter_flux_bins(gaia, regions, limit):
        builder.add_bin(*sky_bin)
    if not builder.cells:
        raise ValueError("Gaia入力に担当領域・等級範囲の天体がありません。")
    if etching.milky_way:
        print("天の川強調: 背景の明るさを差し引いています。", flush=True)
    print(f"最終穴径 {etching.min_diameter_mm:g} mm、金属幅 {etching.min_web_mm:g} mm の穴を配置しています。", flush=True)
    holes, plates = builder.build()
    if builder.milky_way:
        milky_way = builder.milky_way
        print(f"天の川強調: 黒点 {milky_way['black_level']:.3g}、基準点 {milky_way['anchor_level']:.3g}"
              f"（G={gaia.reference:g}等星換算/平方度）、"
              f"{milky_way['dark_cells']:,}/{milky_way['cells']:,}セルを遮光", flush=True)
    if not len(holes):
        raise ValueError("出力できる穴がありません。光量倍率・最小穴径・用紙配置を確認してください。")
    for item in plates:
        requested = item["requested_area_mm2"]
        loss = item["capacity_loss_area_mm2"] / requested if requested else 0.
        print(f"{item['plate']}: {item['holes']:,}穴、容量不足の光量 {loss:.2%}")
        if loss > .01:
            print("  注意: 濃淡が飽和しています。flux-gain（galaxy.etch.* / galaxy.ohp.*）を下げるか、加工条件・集計セルを見直してください。")
    return holes, builder.radius


def generate(config_path, formats=("svg", "pdf"), download=False):
    """設定ファイルに従ってGaia原盤を formats（"svg" / "pdf" / "pdf-mirrored"）で生成し、出力フォルダを返します。

    "pdf-mirrored" は裏面印刷用に左右反転したPDFを ``<名前>-mirrored.pdf`` に出力する。
    """
    config_path = Path(config_path).resolve()
    props = load_properties(config_path)
    directory = config_path.parent
    # 時間のかかる集計の前に、全列の設定を確かめる。
    columns = []
    for number, (change, column_props) in enumerate(_column_properties(props), 1):
        try:
            columns.append((change, GaiaSettings.from_properties(column_props, directory),
                            EtchingSettings.from_properties(column_props), _init_galaxy_transformer(column_props)))
        except ValueError as exc:
            raise ValueError(f"列{number}: {exc}" if change else str(exc)) from exc
    caches = {}
    for number, (_, gaia, _, _) in enumerate(columns, 1):
        first, cached = caches.setdefault(gaia.path, (number, gaia))
        if cached != gaia:
            raise ValueError(f"列{first}と列{number}は galaxy.gaia.input が同じで、Gaiaの等級・集計幅が違います。"
                             "列ごとに galaxy.gaia.input も変えてください。")
    if download:
        for _, gaia in caches.values():
            download_gaia(gaia, GalaxyTransformer.LONGITUDE_REGIONS, GalaxyTransformer.LATITUDE_LIMIT)
    # 読み込み・配置に失敗したとき、既存原盤を途中結果で上書きしない。
    headings = _column_headings([change for change, *_ in columns])
    writers = [_init_plate_writer(props, config_path, fmt == "pdf", headings, fmt == "pdf-mirrored")
               for fmt in formats]
    results = []
    for number, (change, gaia, etching, transformer) in enumerate(columns, 1):
        if len(columns) > 1:
            print(f"列{number}: {', '.join(f'{k} = {v}' for k, v in change.items()) or '変更なし'}", flush=True)
        results.append(_build_column(gaia, etching, transformer, writers[0].r, writers[0].circular))
    for writer in writers:
        with writer:
            writer.write_frames(2)
            for column, (holes, radius) in enumerate(results):
                for dir_, index, x, y in holes:
                    writer.write_hole(dir_, index, x, y, radius, column)
    print(f"生成完了: {writers[0].output_dir}")
    return writers[0].output_dir


def main(argv=None, config="gaia-etch.properties", prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("-f", "--config", type=Path, default=PROJECT_ROOT / config)
    output_format = parser.add_mutually_exclusive_group()
    output_format.add_argument("-SVG", "-svg", "--svg", dest="formats", action="store_const", const=("svg",),
                               help="SVGのみ出力")
    output_format.add_argument("-PDF", "-pdf", "--pdf", dest="formats", action="store_const", const=("pdf",),
                               help="印刷用PDFのみ出力")
    output_format.add_argument("-mirrored", "--mirrored", dest="formats", action="store_const",
                               const=("pdf-mirrored",), help="裏面印刷用に左右反転したPDF（<名前>-mirrored.pdf）のみ出力")
    parser.set_defaults(formats=("svg", "pdf"))
    parser.add_argument("--download-gaia", action="store_true", help="光量マップを取得して生成")
    args = parser.parse_args(argv)
    try:
        generate(args.config, args.formats, args.download_gaia)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"Gaia原盤の生成に失敗しました: {exc}", file=sys.stderr)
        return 2
    return 0
