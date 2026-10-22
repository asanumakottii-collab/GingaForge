# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012,2014 東京大学地文研究会天文部
# Copyright (C) 2026 東京大学地文研究会天文部
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
"""Gaia DR3のG帯光量からエッチング用の天の川原盤（SVG / PDF）を生成する。"""

import argparse
import math
from pathlib import Path
import sys

from .config import PROJECT_ROOT, load_properties, output_directory, positive_number
from .gaia_catalog import GaiaSettings, download_gaia, iter_flux_samples
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
    print(f"円筒投影: ドーム半径 {sphere:g} mm / 原盤距離 {plate:g} mm / 倍率 {scale:g}")
    return GalaxyTransformer(radius, sphere * scale, horizontal * scale, vertical * scale, plate * scale)


def _init_plate_writer(props, config_directory, pdf=False):
    column = int(props.get("plate.column", 1))
    row = int(props.get("plate.row", 1))
    frame_size = positive_number(props, "plate.frame-size", 0., allow_zero=True)
    filename_prefix = props.get("galaxy.gaia.file.prefix", "gaia-etch-")
    invert_color = props.get("color.invert", "no").strip().lower() in {"y", "yes", "true"}
    # 白地に黒穴（normal）と黒地に白穴（inverted）は別フォルダに出力する。
    output_dir = output_directory(config_directory, "gaia-etch", "inverted" if invert_color else "normal")
    writer_class = PlateWriterPDF if pdf else PlateWriterSVG
    return writer_class(column, row, frame_size / 2, filename_prefix, invert_color, output_dir)


def generate(config_path, pdf=False, download=False):
    """設定ファイルに従ってGaia原盤を生成し、出力フォルダを返します。"""
    props = load_properties(config_path)
    directory = Path(config_path).resolve().parent
    gaia = GaiaSettings.from_properties(props, directory)
    etching = EtchingSettings.from_properties(props)
    regions, limit = GalaxyTransformer.LONGITUDE_REGIONS, GalaxyTransformer.LATITUDE_LIMIT
    transformer = _init_galaxy_transformer(props)
    if download:
        download_gaia(gaia, regions, limit)
    # 読み込み・配置に失敗したとき、既存原盤を途中結果で上書きしない。
    writer = _init_plate_writer(props, directory, pdf)
    builder = EtchingBuilder(transformer, etching, writer.r)
    print(f"Gaia G={gaia.bright:g}より暗く{gaia.faint:g}以下の光量を集計します。", flush=True)
    for position, flux in iter_flux_samples(gaia, regions, limit):
        builder.add_sample(position, flux)
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
    with writer:
        writer.write_frames(2)
        for dir_, index, x, y in holes:
            writer.write_hole(dir_, index, x, y, builder.radius)
    for item in plates:
        requested = item["requested_area_mm2"]
        loss = item["capacity_loss_area_mm2"] / requested if requested else 0.
        print(f"{item['plate']}: {item['holes']:,}穴、容量不足の光量 {loss:.2%}")
        if loss > .01:
            print("  注意: 濃淡が飽和しています。galaxy.etch.flux-gainを下げるか、加工条件・集計セルを見直してください。")
    print(f"生成完了: {writer.output_dir}")
    return writer.output_dir


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--config", type=Path, default=PROJECT_ROOT / "gaia-etch.properties")
    parser.add_argument("-PDF", "-pdf", "--pdf", dest="pdf", action="store_true", help="印刷用PDFを出力")
    parser.add_argument("--download-gaia", action="store_true", help="光量マップを取得して生成")
    args = parser.parse_args(argv)
    try:
        generate(args.config, args.pdf, args.download_gaia)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"Gaia原盤の生成に失敗しました: {exc}", file=sys.stderr)
        return 2
    return 0
