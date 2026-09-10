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
"""Gaia DR3のG帯光量からエッチング用の天の川原盤を生成する。"""

import argparse
import json
import math
import os
from pathlib import Path
import sys

from config import load_properties
from gaia_catalog import positive_number
from galaxy_geometry import GalaxyTransformer
from plate_writer import (DEFAULT_OUTPUT_DIR, PlateWriterPDF, PlateWriterSVG,
                          PlateWriterType, categorized_output_dir)

ROOT = Path(__file__).resolve().parent


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


def _init_plate_writer(props, writer_type, write_frames=True):
    column = int(props.get("plate.column", 1))
    row = int(props.get("plate.row", 1))
    frame_size = positive_number(props, "plate.frame-size", 0., allow_zero=True)
    filename_prefix = props.get("galaxy.gaia.file.prefix", "gaia-etch-")
    invert_color = props.get("color.invert", "no").strip().lower() in {"y", "yes", "true"}
    output_dir = categorized_output_dir(props.get("output.directory", DEFAULT_OUTPUT_DIR), "galaxy")
    writer_class = {PlateWriterType.SVG: PlateWriterSVG, PlateWriterType.PDF: PlateWriterPDF}[writer_type]
    writer = writer_class(column, row, frame_size / 2, False, filename_prefix, invert_color, output_dir)
    if write_frames:
        writer.write_frames(2)
    return writer


def _process_gaia(props, config_directory, writer_type, fetch=False, query_only=False):
    from gaia_catalog import (GaiaSettings, build_queries, download_gaia,
                              iter_flux_samples)
    from galaxy_etching import EtchingBuilder, EtchingSettings

    gaia = GaiaSettings.from_properties(props, config_directory)
    etching = EtchingSettings.from_properties(props)
    regions, limit = GalaxyTransformer.LONGITUDE_REGIONS, GalaxyTransformer.LATITUDE_LIMIT
    if query_only:
        for index, query in enumerate(build_queries(gaia, regions, limit)):
            print(f"-- region {index}\n{query};\n")
        return 0
    transformer = _init_galaxy_transformer(props)
    if fetch:
        download_gaia(gaia, regions, limit)
    # 読み込み・配置に失敗したとき、既存原盤を途中結果で上書きしない。
    writer = _init_plate_writer(props, writer_type, write_frames=False)
    try:
        builder = EtchingBuilder(transformer, etching, writer.r)
        input_stats = {}
        print(f"Gaia G={gaia.bright:g}より暗く{gaia.faint:g}以下の光量を集計します。", flush=True)
        print(f"最終穴径 {etching.min_diameter_mm:g} mm、金属幅 {etching.min_web_mm:g} mm（仮の加工条件）", flush=True)
        for position, flux, count in iter_flux_samples(gaia, regions, limit, input_stats):
            builder.add_sample(position, flux, count)
            if input_stats["rows"] % 100000 == 0:
                print(f"  {input_stats['rows']:,}レコード / {input_stats['stars']:,}天体", flush=True)
        if not sum(builder.counts.values()):
            raise ValueError("Gaia入力に担当領域・等級範囲の天体がありません。")
        print("最小穴径・金属幅を満たす穴を配置しています。", flush=True)
        holes, report = builder.build()
        if not holes:
            raise ValueError("出力できる穴がありません。光量倍率・最小穴径・用紙配置を確認してください。")
        report["generator"] = "GingaForge / Gaia DR3 etching"
        report["config"] = dict(props)
        report["frame_size_output_mm"] = 2 * writer.r
        report["input"] = input_stats
        report["geometry"] = {k: getattr(transformer, k) for k in (
            "radius", "sphere", "projector_horizontal", "projector_vertical", "plate")}
        writer.write_frames(2)
        for hole in holes:
            writer.write_star(hole)
    finally:
        writer.close()
    report_path = Path(writer.output_dir) / f"{writer.filename_prefix}report.json"
    partial = report_path.with_suffix(".json.part")
    partial.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, report_path)
    for item in report["plates"]:
        requested = item["requested_area_mm2"]
        loss = item["capacity_loss_area_mm2"] / requested if requested else 0.
        print(f"{item['plate']}: {item['source_stars']:,}天体 → {item['holes']:,}穴、容量不足の光量 {loss:.2%}")
        if loss > .01:
            print("  注意: 濃淡が飽和しています。galaxy.etch.flux-gainを下げるか、加工条件・集計セルを見直してください。")
    print(f"生成完了。加工条件と光量の記録: {report_path}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--config", type=Path, default=ROOT / "gaia-etch.properties")
    parser.add_argument("-PDF", "-pdf", "--pdf", dest="pdf", action="store_true", help="印刷用PDFを出力")
    parser.add_argument("--download-gaia", action="store_true", help="光量マップを取得して生成")
    parser.add_argument("--gaia-query", action="store_true", help="ADQLクエリだけを表示（通信なし）")
    parser.add_argument("-PS", "-ps", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.PS:
        print("PostScript出力は廃止されました。-PDFを使用してください。", file=sys.stderr)
        return 2
    try:
        props = load_properties(args.config)
        if props.get("galaxy.mode", "gaia-etch").strip().lower() != "gaia-etch":
            raise ValueError("このコマンドはgaia-etch専用です。galaxy.modeは省略するかgaia-etchを指定してください。")
        directory = args.config.resolve().parent
        output = Path(props.get("output.directory", DEFAULT_OUTPUT_DIR).strip() or DEFAULT_OUTPUT_DIR)
        props["output.directory"] = str(output if output.is_absolute() else directory / output)
        print("GingaForge — Gaia DR3 etching (OTL-derived, GPL-2.0-or-later)")
        return _process_gaia(props, directory, PlateWriterType.PDF if args.pdf else PlateWriterType.SVG,
                             args.download_gaia, args.gaia_query)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"Gaia原盤の生成に失敗しました: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
