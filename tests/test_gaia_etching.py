import base64
from collections import Counter
import csv
import gzip
from io import BytesIO
import json
import math
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
from io import StringIO
from unittest.mock import patch

from gingaforge.gaia_catalog import (GaiaSettings, build_queries, iter_flux_bins, iter_flux_samples,
                          _votable_to_csv, file_sha256)
from gingaforge.galaxy_etching import EtchingBuilder, EtchingSettings, dither_pixels
from gingaforge.gaia_etch import GalaxyTransformer, main
from gingaforge.models import SpherePosition
import numpy as np
from PIL import Image


def position(lon, lat=0.):
    return SpherePosition.from_galactic(lon, lat)


def lens_at_projector(horizontal=0., vertical=0., lens=50.):
    """ドームと投影機の位置を10^9倍にし、レンズ中心を投影機中心とみなせる変換（円筒投影の半径は約 lens mm）。"""
    k = 1e9
    return GalaxyTransformer(.1, 6500. * k, horizontal * k, vertical * k, 0., lens)


def write_cache(directory, cells, settings=None):
    """``{領域番号: [(銀経, 銀緯, 光量, 天体数), ...]}`` から検証可能なGaiaキャッシュを作る。"""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    settings = settings or GaiaSettings(directory)
    regions = GalaxyTransformer.LONGITUDE_REGIONS
    files = []
    for index in range(len(regions)):
        path = directory / f"region-{index}.csv.gz"
        with gzip.open(path, "wt", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["l", "b", "relative_flux", "source_count"])
            writer.writerows(cells.get(index, []))
        files.append({"name": path.name, "sha256": file_sha256(path)})
    manifest = {"complete": True, "files": files,
                "selection": settings.provenance(regions, GalaxyTransformer.LATITUDE_LIMIT)}
    (directory / "manifest.json").write_text(json.dumps(manifest))
    return directory


class EtchingTests(unittest.TestCase):
    def setUp(self):
        self.t = GalaxyTransformer(.1, 6500., 0., 0., 50., 19.)
        self.settings = EtchingSettings()

    def test_inverse_projection_with_offset(self):
        for horizontal, vertical in ((0., 0.), (355., 200.)):
            transformer = GalaxyTransformer(.1, 6500., horizontal, vertical, 50., 19.)
            for lon in (0., 30., 59.9, 60., 119., 240., 270., 300., 330., 359.9):
                for lat in (-20., 0., 20.):
                    sp = position(lon, lat)
                    pp = transformer.transform(sp)
                    back = transformer.inverse_transform_unit(pp.xmm, pp.ymm, pp.dir, pp.index)
                    actual_l, actual_b = back.to_galactic()
                    self.assertAlmostEqual(0., (actual_l - lon + 180) % 360 - 180, places=9)
                    self.assertAlmostEqual(lat, actual_b, places=9)

    def test_sixteen_faint_stars_combine_to_one_manufacturable_hole(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        # 基準穴の面積の1/100が最小穴。さらに1/16の星を16個。
        for _ in range(16):
            builder.add_sample(position(30.01, .01), .01 / 16)
        holes, plates = builder.build()
        self.assertEqual(1, len(holes))
        self.assertAlmostEqual(.02, 2 * builder.radius)
        self.assertEqual([1, 0, 0, 0], [p["holes"] for p in plates])

    def test_flux_sum_and_spacing_across_adjacent_cells(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        for ix in range(-2, 3):
            for iy in range(-2, 3):
                sp = self.t.inverse_transform_unit((ix + .5) * builder.cell,
                                                   (iy + .5) * builder.cell, 0, 0)
                builder.add_sample(sp, .023)
        holes, plates = builder.build()
        self.assertEqual(58, len(holes))  # 25 * 2.3 最小穴相当を四捨五入
        points = [(x, y) for _, _, x, y in holes]
        for i, first in enumerate(points):
            for second in points[i + 1:]:
                distance = math.hypot(first[0] - second[0], first[1] - second[1])
                self.assertGreaterEqual(distance - 2 * builder.radius + 1e-12, .02)
        self.assertAlmostEqual(0., plates[0]["capacity_loss_area_mm2"])

    def test_saturation_is_reported_and_not_spilled_to_dark_cells(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        builder.add_sample(position(30.01, .01), 100.)
        holes, plates = builder.build()
        self.assertEqual(25, len(holes))
        self.assertGreater(plates[0]["capacity_loss_area_mm2"], 0)
        self.assertEqual(1, len({(math.floor(x / builder.cell), math.floor(y / builder.cell))
                                 for _, _, x, y in holes}))

    def test_determinism_and_parameter_changes(self):
        def run(props):
            settings = EtchingSettings.from_properties(props)
            builder = EtchingBuilder(self.t, settings, 69.25)
            builder.add_sample(position(30.01, .01), .2)
            return builder, list(builder.build()[0])
        _, first = run({})
        _, second = run({})
        self.assertEqual(first, second)
        builder, wider = run({"galaxy.etch.min-hole-diameter-mm": ".04"})
        self.assertEqual(5, len(wider))
        self.assertTrue(math.isclose(builder.radius, .02))

    def test_final_dimensions_follow_scale(self):
        settings = EtchingSettings.from_properties({"scale": "2"})
        builder = EtchingBuilder(GalaxyTransformer(.2, 13000., 0., 0., 100., 38.), settings, 69.25)
        builder.add_sample(position(30.01, .01), .1)
        holes, _ = builder.build()
        self.assertEqual(10, len(holes))
        self.assertTrue(math.isclose(builder.radius / settings.scale, .01))
        self.assertAlmostEqual(.08, builder.pitch)  # (穴径 + 金属幅) * scale

    def test_outside_regions_and_physical_frame(self):
        builder = EtchingBuilder(self.t, self.settings, .1)
        builder.add_sample(position(150.), 1.)
        builder.add_sample(position(30., 21.), 1.)
        builder.add_sample(position(30., 10.), 1.)
        builder.add_sample(position(30.01, .01), .1)
        holes, plates = builder.build()
        self.assertEqual({(0, 0)}, {key[:2] for key in builder.cells})
        self.assertGreater(plates[0]["capacity_loss_area_mm2"], 0)
        self.assertTrue(all(max(abs(x), abs(y)) + builder.radius + .02 <= .1 for _, _, x, y in holes))

    def test_circular_frame_keeps_holes_inside_circle(self):
        # 直径10 mmの円。外接正方形の角に当たる光も足すが、円の外には穴を置かない。
        builder = EtchingBuilder(self.t, self.settings, 5.)
        for lon in range(25, 36):
            for lat in range(-5, 6):
                builder.add_sample(position(lon + .01, lat + .01), 1.)
        holes = list(builder.build()[0])
        margin = builder.radius + .02
        self.assertTrue(holes)
        self.assertTrue(all(math.hypot(x, y) + margin <= 5. for _, _, x, y in holes))

    def test_four_units_do_not_share_rounding_or_flux(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        for lon in (30.01, 90.01, 270.01, 330.01):
            builder.add_sample(position(lon, .01), .006)
        holes, plates = builder.build()
        self.assertEqual(4, len(holes))
        self.assertEqual([1] * 4, [p["holes"] for p in plates])

    @staticmethod
    def uniform_sky(builders, lon, lat, add_point=None):
        """銀経・銀緯の下端から0.1度のビンを200 x 60個敷き詰め、各ビンに光量1を置く。"""
        for i in range(200):
            for j in range(60):
                edges = (lon + .1 * i, lon + .1 * (i + 1), lat + .1 * j, lat + .1 * (j + 1))
                for builder in builders:
                    if builder is add_point:
                        builder.add_sample(position((edges[0] + edges[1]) / 2, (edges[2] + edges[3]) / 2), 1.)
                    else:
                        builder.add_bin(*edges, 1.)
        return 200 * 60

    @staticmethod
    def per_solid_angle(builder):
        """セルの光量を、セルが写る空の広さで割ったもの（平らな原板では端ほど空が原盤上で広がる）。"""
        return {key: value / builder._solid_angle_deg2(key[:2], *key[2:]) for key, value in builder.cells.items()}

    @staticmethod
    def largest_step_between_neighbours(cells):
        """光量のあるセル（端の3列・3行を除く）で、銀経方向に隣り合うセルの光量の比の最大のずれ。"""
        ixs, iys = sorted({k[2] for k in cells}), sorted({k[3] for k in cells})
        unit = next(iter(cells))[:2]
        worst = 0.
        for ix in ixs[3:-3]:
            column = [cells.get((*unit, ix, iy)) for iy in iys[3:-3]]
            worst = max(worst, *(abs(b / a - 1) for a, b in zip(column, column[1:]) if a and b))
        return worst

    def test_uniform_sky_bins_fill_cells_without_moire(self):
        # 0.1度のビンは原盤上で約0.087 mm。0.4 mmのセルには4個または5個ずつ入るので、
        # ビンの中心1点に光量を置くと、一様な空でも隣のセルと5/4倍違う縞ができる。
        settings = EtchingSettings(cell_mm=.4)
        area = EtchingBuilder(lens_at_projector(), settings, 69.25)
        point = EtchingBuilder(lens_at_projector(), settings, 69.25)
        count = self.uniform_sky((area, point), 20., -3., add_point=point)
        self.assertAlmostEqual(1., sum(area.cells.values()) / (count * area.reference_area), places=12)
        self.assertAlmostEqual(.25, self.largest_step_between_neighbours(point.cells), places=9)
        # セルの立体角は差分による近似なので、その誤差（1e-4程度）は残る。縞の0.25よりずっと小さい。
        self.assertLess(self.largest_step_between_neighbours(self.per_solid_angle(area)), 1e-3)
        # 銀緯方向も含め、内側のセルは写る空の広さあたり同じ光量。
        inner = [v for (_, _, ix, iy), v in self.per_solid_angle(area).items() if -5 <= ix <= 4 and -20 <= iy <= 19]
        self.assertLess(max(inner) / min(inner), 1.01)

    def test_uniform_sky_bins_with_projector_offset(self):
        transformer = lens_at_projector(355., 200.)
        area = EtchingBuilder(transformer, EtchingSettings(cell_mm=.4), 69.25)
        point = EtchingBuilder(transformer, EtchingSettings(cell_mm=.4), 69.25)
        count = self.uniform_sky((area, point), 80., 10., add_point=point)
        self.assertEqual({(0, 1)}, {key[:2] for key in area.cells})
        self.assertAlmostEqual(1., sum(area.cells.values()) / (count * area.reference_area), places=12)
        self.assertGreater(self.largest_step_between_neighbours(point.cells), .2)
        self.assertLess(self.largest_step_between_neighbours(self.per_solid_angle(area)), .005)

    def test_bin_is_split_by_overlap_area(self):
        builder = EtchingBuilder(lens_at_projector(), EtchingSettings(cell_mm=.4), 69.25)
        # 原盤中心（セルの角）をまたぐビン。レンズ中心がドーム中心なら、平らな原板では
        # x = L tan b / cos(l - 30度)、y = -L tan(l - 30度)。ビンの像は各辺の両端の平均の長方形とみなす。
        builder.add_bin(29.96, 30.06, -.04, .06, 1.)
        tan, cos = (lambda f: lambda degrees: f(math.radians(degrees)))(math.tan), (lambda d: math.cos(math.radians(d)))
        edge = lambda b: (tan(b) / cos(-.04) + tan(b) / cos(.06)) / 2
        x_fraction = {-1: -edge(-.04) / (edge(.06) - edge(-.04)), 0: edge(.06) / (edge(.06) - edge(-.04))}
        y_fraction = {-1: tan(.06) / (tan(.04) + tan(.06)), 0: tan(.04) / (tan(.04) + tan(.06))}
        expected = {(0, 0, ix, iy): x_fraction[ix] * y_fraction[iy] for ix in (-1, 0) for iy in (-1, 0)}
        self.assertEqual(set(expected), set(builder.cells))
        for key, fraction in expected.items():
            self.assertAlmostEqual(fraction, builder.cells[key] / builder.reference_area, places=9)

    def test_bins_outside_regions_degenerate_and_invalid(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        builder.add_bin(150., 150.1, 0., .1, 1.)
        builder.add_bin(30., 30.1, 20.5, 20.6, 1.)
        self.assertEqual({}, dict(builder.cells))
        # 銀緯20度ちょうどの天体だけを集計した幅0のビンは、点として1つのセルに入る。
        builder.add_bin(30., 30.1, 20., 20., 1.)
        self.assertEqual(1, len(builder.cells))
        self.assertAlmostEqual(builder.reference_area, sum(builder.cells.values()))
        with self.assertRaises(ValueError):
            builder.add_bin(30., 30.1, 0., .1, -1.)

    def test_invalid_manufacturing_settings_fail(self):
        for key, value in (("galaxy.etch.min-hole-diameter-mm", "0"),
                           ("galaxy.etch.min-web-mm", "-1"),
                           ("galaxy.etch.flux-gain", "nan"),
                           ("scale", "inf"), ("galaxy.etch.seed", "1.2"),
                           ("galaxy.etch.cell-size-mm", ".01"),
                           ("galaxy.etch.milky-way.black-percentile", "-1"),
                           ("galaxy.etch.milky-way.black-percentile", "95"),
                           ("galaxy.etch.milky-way.anchor-percentile", "101"),
                           ("galaxy.etch.milky-way.gamma", "0"),
                           ("galaxy.etch.milky-way.smoothing-mm", "-1"),
                           ("galaxy.etch.milky-way.smoothing-mm", "100")):
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                EtchingSettings.from_properties({key: value})

    def milky_way_builder(self, **settings):
        # 原盤上 ±8 mm の範囲を使うので、空の角度が小さく済む原板とレンズの距離50 mmで調べる。
        transformer = lens_at_projector()
        builder = EtchingBuilder(transformer, EtchingSettings(**settings), 69.25)
        # 銀緯方向(x)に明るさが落ちる帯と、外側の明るい星1つ。
        for ix in range(-20, 20):
            for iy in range(-10, 10):
                x = (ix + .5) * builder.cell
                flux = .02 + .1 * math.exp(-(x / 1.5) ** 2) + (.4 if (ix, iy) == (-18, 0) else 0.)
                builder.add_sample(transformer.inverse_transform_unit(x, (iy + .5) * builder.cell, 0, 0), flux)
        return builder

    def test_milky_way_enhancement_darkens_background_and_keeps_band(self):
        linear_builder = self.milky_way_builder()
        linear_holes, linear_plates = linear_builder.build()
        builder = self.milky_way_builder(milky_way=True)
        holes, plates = builder.build()
        column = lambda found: Counter(math.floor(x / linear_builder.cell) for _, _, x, _ in found)
        linear, enhanced = column(linear_holes), column(holes)
        outside = [*range(-20, -14), *range(15, 20)]  # 明るい星の列 -18 を含む
        self.assertGreater(sum(linear[ix] for ix in outside), 0)
        self.assertEqual(0, sum(enhanced[ix] for ix in outside))
        self.assertGreaterEqual(enhanced[-1] + enhanced[0], linear[-1] + linear[0])
        info = builder.milky_way
        self.assertEqual(5, info["smoothing_window_cells"])
        self.assertLess(info["black_level"], info["anchor_level"])
        self.assertGreater(info["dark_cells"], 0)
        plate = plates[0]
        self.assertAlmostEqual(linear_plates[0]["requested_area_mm2"], plate["linear_area_mm2"])
        self.assertLess(plate["requested_area_mm2"], plate["linear_area_mm2"])
        self.assertIsNone(linear_builder.milky_way)

    def test_ohp_prefix_is_same_settings_as_etch_prefix(self):
        etch = EtchingSettings.from_properties({"galaxy.etch.cell-size-mm": ".4", "galaxy.etch.milky-way": "yes"})
        ohp = EtchingSettings.from_properties({"galaxy.ohp.cell-size-mm": ".4", "galaxy.ohp.milky-way": "yes"})
        self.assertEqual(etch, ohp)
        with self.assertRaisesRegex(ValueError, "galaxy.ohp.seed"):
            EtchingSettings.from_properties({"galaxy.ohp.seed": "1.2"})
        with self.assertRaisesRegex(ValueError, "どちらか一方"):
            EtchingSettings.from_properties({"galaxy.etch.seed": "1", "galaxy.ohp.seed": "2"})

    def test_milky_way_enhancement_settings_and_uniform_input(self):
        settings = EtchingSettings.from_properties({"galaxy.etch.milky-way": "yes",
                                                    "galaxy.etch.milky-way.gamma": "2"})
        self.assertTrue(settings.milky_way)
        self.assertEqual(2., settings.milky_way_gamma)
        self.assertFalse(EtchingSettings.from_properties({}).milky_way)
        builder = EtchingBuilder(self.t, settings, 69.25)
        builder.add_sample(position(30.01, .01), .2)
        with self.assertRaises(ValueError):
            builder.build()

    def test_render_settings(self):
        self.assertEqual(("holes", 1200.), (self.settings.render, self.settings.raster_dpi))
        gray = EtchingSettings.from_properties({"galaxy.ohp.render": "Gray", "galaxy.ohp.raster-dpi": "600"})
        self.assertEqual(("gray", 600.), (gray.render, gray.raster_dpi))
        # 画像には穴の候補格子がないため、穴径＋金属幅（0.04 mm）より細かいセルも使える。
        fine = {"galaxy.ohp.cell-size-mm": ".03"}
        with self.assertRaisesRegex(ValueError, "cell-size-mm"):
            EtchingSettings.from_properties(fine)
        self.assertEqual(.03, EtchingSettings.from_properties({**fine, "galaxy.ohp.render": "dither"}).cell_mm)
        builder = EtchingBuilder(self.t, EtchingSettings(cell_mm=.13, render="gray"), 69.25)
        self.assertAlmostEqual(.13, builder.cell)
        for key, value in (("render", "stipple"), ("raster-dpi", "50"), ("raster-dpi", "4800"), ("raster-dpi", "x")):
            with self.subTest(key=key, value=value), self.assertRaisesRegex(ValueError, f"galaxy.ohp.{key}"):
                EtchingSettings.from_properties({f"galaxy.ohp.{key}": value})

    def raster_of_uniform_sky(self, render):
        """一様な空（銀経20〜40度、銀緯±3度）の N0 の画像と集計、内側の画素を返す。"""
        transformer = lens_at_projector()
        builder = EtchingBuilder(transformer, EtchingSettings(cell_mm=.13, flux_gain=.1, render=render,
                                                              raster_dpi=600), 69.25)
        self.uniform_sky((builder,), 20., -3.)
        images, plates = builder.build_raster()
        self.assertEqual([(0, 0)], [image[:2] for image in images])
        _, _, pixels, left, top, width, height = images[0]
        pixel = 25.4 / 600
        self.assertEqual((round(height / pixel), round(width / pixel)), pixels.shape)
        x = left + (np.arange(pixels.shape[1]) + .5) * pixel
        y = top + (np.arange(pixels.shape[0]) + .5) * pixel
        # データの端（銀経20・40度、中央の銀緯±3度。平らな原板では銀緯の端は両端ほど外に広がる）からセル3個分より内側。
        y_edge = abs(transformer.transform(position(40., 0.)).ymm)
        x_edge = abs(transformer.transform(position(30., 3.)).xmm)
        margin = 3 * builder.cell
        inner = np.ix_(np.abs(y) < y_edge - margin, np.abs(x) < x_edge - margin)
        return pixels, plates[0], pixels[inner]

    def test_gray_raster_of_uniform_sky_is_uniform_and_keeps_light(self):
        pixels, plate, inner = self.raster_of_uniform_sky("gray")
        # ビン（約0.087 mm）とセル（0.13 mm）の比が整数でなくても、縞は出ない（隣り合う画素の差は階調1段以内）。
        # 平らな原板では空が端ほど広く写るため、開口率は端に向かってなめらかに少し下がる。
        self.assertLessEqual(np.abs(np.diff(inner.astype(int), axis=0)).max(), 1)
        self.assertLessEqual(np.abs(np.diff(inner.astype(int), axis=1)).max(), 1)
        self.assertAlmostEqual(.4, inner.mean() / 255, delta=.05)
        self.assertAlmostEqual(1., plate["open_area_mm2"] / plate["requested_area_mm2"], delta=.005)
        self.assertAlmostEqual(plate["requested_area_mm2"], plate["linear_area_mm2"])
        self.assertEqual(0., plate["capacity_loss_area_mm2"])
        self.assertAlmostEqual(inner.max() / 255, plate["max_open_fraction"], delta=.01)

    def test_dither_raster_of_uniform_sky_keeps_light(self):
        pixels, plate, inner = self.raster_of_uniform_sky("dither")
        _, _, gray_inner = self.raster_of_uniform_sky("gray")
        self.assertEqual({0, 255}, set(np.unique(pixels).tolist()))
        self.assertAlmostEqual(gray_inner.mean() / 255, (inner == 255).mean(), delta=.005)
        self.assertAlmostEqual(1., plate["open_area_mm2"] / plate["requested_area_mm2"], delta=.005)

    def test_raster_saturation_is_clipped_and_not_spilled(self):
        builder = EtchingBuilder(self.t, EtchingSettings(cell_mm=.2, render="gray", raster_dpi=600), 69.25)
        builder.add_sample(position(30.05, .05), 1e4)
        builder.add_sample(position(31.05, .05), .01)
        images, plates = builder.build_raster()
        plate = plates[0]
        self.assertEqual(1., plate["max_open_fraction"])
        self.assertAlmostEqual(plate["requested_area_mm2"] - builder.cell ** 2 - builder.reference_area * .01,
                               plate["capacity_loss_area_mm2"])
        _, _, pixels, left, top, _, _ = images[0]
        pixel = 25.4 / 600
        # 明るいセルの中心（.1, -.1）。画素の中心はセルの中心と一致しないため、補間で少し暗くなる。
        row, col = int((-.1 - top) / pixel), int((.1 - left) / pixel)
        self.assertGreater(pixels[row, col], 200)
        self.assertEqual(pixels.max(), pixels[row - 1:row + 2, col - 1:col + 2].max())
        # 飽和した光は、暗いセルへ移さない。
        far = builder.transformer.transform(position(31.05, .05))
        row, col = int((far.ymm - top) / pixel), int((far.xmm - left) / pixel)
        self.assertLess(pixels[row, col], 20)

    def test_dither_pixels(self):
        quarter = dither_pixels(np.full((60, 80), .25))
        self.assertAlmostEqual(.25, (quarter == 255).mean(), delta=.005)
        np.testing.assert_array_equal(quarter, dither_pixels(np.full((60, 80), .25)))
        self.assertFalse(dither_pixels(np.zeros((5, 7))).any())
        self.assertTrue((dither_pixels(np.ones((5, 7))) == 255).all())
        # 暗い所でも、点の数は光量の合計に近い。8ビットに丸めると 0.001 は0になり、点が1つも打たれない。
        for value in (.003, .001):
            faint = dither_pixels(np.full((1000, 100), value))
            self.assertAlmostEqual(value * faint.size, (faint == 255).sum(), delta=1)
        # 光のない所（開口率0）へ誤差を配らないので、明るい四角の外に点がなく、四角の中の光も減らない。
        square = np.zeros((50, 60))
        square[10:30, 20:45] = .1
        dots = dither_pixels(square) == 255
        self.assertFalse(dots[square == 0].any())
        self.assertAlmostEqual(.1 * 20 * 25, dots.sum(), delta=1)


class GaiaInputTests(unittest.TestCase):
    def test_missing_cache_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, "--download-gaia"):
                list(iter_flux_samples(GaiaSettings(Path(directory)), ((0, 60),), 20))

    def test_queries_use_flux_and_disjoint_regions_without_sampling(self):
        queries = build_queries(GaiaSettings(Path("unused")), GalaxyTransformer.LONGITUDE_REGIONS, 20)
        self.assertEqual(4, len(queries))
        for query in queries:
            self.assertIn("SUM(POWER(10.0, -0.4", query)
            self.assertIn("COUNT(*)", query)
            self.assertNotIn("TOP", query)
        self.assertIn("l >= 300 AND l < 360", queries[2])

    def votable(self, status="OK"):
        return ('<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE>'
                '<INFO name="QUERY_STATUS" value="OK"/><TABLE>'
                + ''.join(f'<FIELD name="{name}"/>' for name in
                          ("lon_bin", "lat_bin", "source_count", "relative_flux"))
                + '<DATA><TABLEDATA><TR><TD>300.0</TD><TD>900.0</TD><TD>16</TD><TD>0.01</TD>'
                  '</TR></TABLEDATA></DATA></TABLE>'
                + f'<INFO name="QUERY_STATUS" value="{status}"/></RESOURCE></VOTABLE>')

    def test_votable_overflow_never_accepted_as_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            source, target = Path(directory) / "data.xml", Path(directory) / "data.csv.gz"
            source.write_text(self.votable("OVERFLOW"))
            with self.assertRaises(RuntimeError):
                _votable_to_csv(source, target, GaiaSettings(Path(directory)), (0, 60), 20)

    def test_verified_cache_and_selection_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source, target = directory / "data.xml", directory / "region-0.csv.gz"
            source.write_text(self.votable())
            settings = GaiaSettings(directory)
            self.assertEqual((1, 16), _votable_to_csv(source, target, settings, (0, 60), 20))
            manifest = {"complete": True, "selection": settings.provenance(((0, 60),), 20),
                        "files": [{"name": target.name, "sha256": file_sha256(target)}]}
            (directory / "manifest.json").write_text(json.dumps(manifest))
            samples = list(iter_flux_samples(settings, ((0, 60),), 20))
            self.assertEqual(1, len(samples))
            self.assertAlmostEqual(30.05, samples[0][0].to_galactic()[0])
            self.assertAlmostEqual(.01, samples[0][1])
            with self.assertRaises(ValueError):
                list(iter_flux_samples(GaiaSettings(directory, faint=18), ((0, 60),), 20))
            target.write_bytes(b"broken")
            with self.assertRaises(ValueError):
                list(iter_flux_samples(settings, ((0, 60),), 20))

    def test_cache_bins_restore_edges_including_clipped_latitude(self):
        regions = GalaxyTransformer.LONGITUDE_REGIONS
        cases = (
            # 0.1度：通常のビンと、銀緯20度ちょうどの天体だけを集計した幅0のビン。
            (.1, [(30.05, .05), (30.05, -19.95), (30.05, 20.)],
             [(30., 30.1, 0., .1), (30., 30.1, -20., -19.9), (30., 30.1, 20., 20.)]),
            # 0.12度：ビンの格子が銀緯±20度にそろわず、両端のビンが切り詰められる。
            (.12, [(30.06, .06), (30.06, -19.96), (30.06, 19.96)],
             [(30., 30.12, 0., .12), (30., 30.12, -20., -19.92), (30., 30.12, 19.92, 20.)]),
        )
        for step, centres, expected in cases:
            with self.subTest(step=step), tempfile.TemporaryDirectory() as directory:
                settings = GaiaSettings(Path(directory), sky_step=step)
                write_cache(directory, {0: [(l, b, .01, 1) for l, b in centres]}, settings)
                bins = list(iter_flux_bins(settings, regions, 20))
                self.assertEqual(len(expected), len(bins))
                for actual, edges in zip(bins, expected):
                    for value, edge in zip(actual, edges):
                        self.assertAlmostEqual(edge, value, places=9)
                    self.assertAlmostEqual(.01, actual[4])


class GaiaCommandTests(unittest.TestCase):
    def config_and_input(self, directory):
        directory = Path(directory)
        write_cache(directory / "gaia", {0: [(30.05, .05, .01, 1)]})
        config = directory / "galaxy.properties"
        config.write_text("galaxy.gaia.input=gaia\ncolor.invert=yes\n")
        return config

    def test_real_svg_pipeline_has_no_prompts(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            with (patch("builtins.input", side_effect=AssertionError("unexpected prompt")),
                  patch("sys.stdout", new_callable=StringIO)):
                self.assertEqual(0, main(["-svg", "-f", str(config)]))
            output = Path(directory) / "output" / "galaxy" / "inverted"
            pages = sorted(output.glob("*.svg"))
            self.assertEqual(2, len(pages))
            self.assertEqual(pages, sorted(output.iterdir()))
            holes = [c for c in ET.parse(pages[0]).getroot().findall("{http://www.w3.org/2000/svg}circle")
                     if not c.get("stroke")]
            self.assertGreater(len(holes), 0)
            self.assertTrue(all(math.isclose(float(c.attrib["r"][:-2]), .01) for c in holes))

    def test_missing_input_does_not_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            output = Path(directory) / "output" / "galaxy" / "inverted"
            output.mkdir(parents=True)
            sentinel = output / "galaxy-0.svg"
            sentinel.write_text("existing plate")
            (Path(directory) / "gaia" / "manifest.json").unlink()
            with patch("sys.stdout", new_callable=StringIO), patch("sys.stderr", new_callable=StringIO):
                self.assertEqual(2, main(["-f", str(config)]))
            self.assertEqual("existing plate", sentinel.read_text())

    def test_pdf_pipeline(self):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            with patch("sys.stdout", new_callable=StringIO):
                self.assertEqual(0, main(["-PDF", "-f", str(config)]))
            pages = sorted((Path(directory) / "output" / "galaxy" / "inverted").glob("*.pdf"))
            self.assertEqual(2, len(pages))
            self.assertEqual(["N0", "S0"], PdfReader(pages[0]).pages[0].extract_text().split())
            self.assertEqual(pages, sorted((Path(directory) / "output" / "galaxy" / "inverted").iterdir()))

    def column_config(self, directory, settings):
        directory = Path(directory)
        write_cache(directory / "gaia", {0: [(30.05, .05, .05, 1)], 3: [(270.05, .05, .05, 1)]})
        config = directory / "galaxy.properties"
        config.write_text("galaxy.gaia.input=gaia\ncolor.invert=yes\nplate.row=2\nplate.frame-diameter=22\n"
                          "galaxy.etch.flux-gain=1\n" + settings)
        return config

    def test_columns_side_by_side_with_changed_settings(self):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            config = self.column_config(directory, "column.count=2\ncolumn.2.galaxy.etch.flux-gain=2\n")
            with patch("sys.stdout", new_callable=StringIO) as stdout:
                self.assertEqual(0, main(["-f", str(config)]))
            self.assertIn("列2: galaxy.etch.flux-gain = 2", stdout.getvalue())
            output = Path(directory) / "output" / "galaxy" / "inverted"
            # 上段 N0・N1、下段 S0・S1 の1組を2列並べ、1用紙に収める。
            self.assertEqual(["galaxy-0.pdf", "galaxy-0.svg"], sorted(p.name for p in output.iterdir()))
            text = PdfReader(output / "galaxy-0.pdf").pages[0].extract_text()
            for label in ("#1", "(no changes)", "#2", "galaxy.etch.flux-gain = 2", "N0 #1", "S1 #2"):
                self.assertIn(label, text)
            svg = ET.parse(output / "galaxy-0.svg").getroot()
            namespace = "{http://www.w3.org/2000/svg}"
            texts = [element.text for element in svg.iter(namespace + "text")]
            self.assertEqual(["#1", "(no changes)", "#2", "galaxy.etch.flux-gain = 2"], texts[:4])
            self.assertEqual(8, sum(" #" in t for t in texts))
            # 左の列は x < 105 mm、右の列は x > 105 mm。光量倍率2の右の列は穴が約2倍。
            sides = Counter(float(c.attrib["cx"][:-2]) > 105 for c in svg.iter(namespace + "circle"))
            self.assertGreater(sides[False], 0)
            self.assertAlmostEqual(2., sides[True] / sides[False], delta=.2)

    def test_single_column_output_is_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            plain, counted = (self.column_config(Path(directory) / name, settings)
                              for name, settings in (("plain", ""), ("counted", "column.count=1\n")))
            with patch("sys.stdout", new_callable=StringIO):
                self.assertEqual(0, main(["-svg", "-f", str(plain)]))
                self.assertEqual(0, main(["-svg", "-f", str(counted)]))
            pages = [config.parent / "output" / "galaxy" / "inverted" / "galaxy-0.svg" for config in (plain, counted)]
            self.assertEqual(pages[0].read_text(), pages[1].read_text())
            texts = [e.text for e in ET.parse(pages[0]).getroot().iter("{http://www.w3.org/2000/svg}text")]
            self.assertEqual(["N0", "S0", "N1", "S1"], texts)

    def test_invalid_column_settings_fail_before_output(self):
        cases = (
            ("column.count=0\n", "column.count"),
            ("column.count=two\n", "column.count"),
            ("column.count=2\ncolumn.3.galaxy.etch.flux-gain=2\n", "1〜2"),
            ("column.count=2\ncolumn.x.galaxy.etch.flux-gain=2\n", "列番号"),
            ("column.count=2\ncolumn.2.galaxy.etch.flux-gian=2\n", "列の指定なし"),
            ("column.count=2\ncolumn.2.plate.frame-diameter=10\n", "全列共通"),
            ("column.count=2\ncolumn.2.color.invert=no\n", "全列共通"),
            ("column.count=2\ncolumn.2.galaxy.etch.flux-gain=2\ncolumn.02.galaxy.etch.flux-gain=3\n", "重複"),
            ("column.count=2\ncolumn.2.galaxy.etch.flux-gain=-1\n", "列2: galaxy.etch.flux-gain"),
            ("column.count=2\ngalaxy.gaia.faint-limit=17\ncolumn.2.galaxy.gaia.faint-limit=16\n",
             "galaxy.gaia.input"),
            # 原盤が列の幅（210 mm ÷ 列の数）に収まらない。
            ("column.count=10\n", "plate.frame-diameter"),
            ("column.count=4\nplate.frame-diameter=60\n", "plate.frame-diameter"),
            # 原盤は収まるが、列の見出しが最上段の原盤名と重なる。
            ("column.count=2\nplate.frame-diameter=60\ncolumn.2.galaxy.etch.flux-gain=2\n", "列の見出し"),
        )
        for settings, message in cases:
            with (self.subTest(settings=settings), tempfile.TemporaryDirectory() as directory,
                  patch("sys.stdout", new_callable=StringIO), patch("sys.stderr", new_callable=StringIO) as stderr):
                config = self.column_config(directory, settings)
                self.assertEqual(2, main(["-f", str(config)]))
                self.assertIn(message, stderr.getvalue())
                self.assertFalse((Path(directory) / "output").exists())

    @staticmethod
    def svg_images(path):
        """SVGの画像を ``(x mm, 画素の配列)`` で返す。"""
        found = []
        for element in ET.parse(path).getroot().iter("{http://www.w3.org/2000/svg}image"):
            data = element.attrib["{http://www.w3.org/1999/xlink}href"].removeprefix("data:image/png;base64,")
            found.append((float(element.attrib["x"][:-2]), np.asarray(Image.open(BytesIO(base64.b64decode(data))))))
        return found

    def test_columns_mix_holes_and_images(self):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            config = self.column_config(directory, "galaxy.etch.render=holes\ncolumn.count=2\n"
                                                   "column.2.galaxy.etch.render=gray\n")
            with patch("sys.stdout", new_callable=StringIO) as stdout:
                self.assertEqual(0, main(["-f", str(config)]))
                self.assertEqual(0, main(["-mirrored", "-f", str(config)]))
            self.assertIn("1200 dpi の階調の画像", stdout.getvalue())
            output = Path(directory) / "output" / "galaxy" / "inverted"
            svg = output / "galaxy-0.svg"
            holes = [c for c in ET.parse(svg).getroot().iter("{http://www.w3.org/2000/svg}circle")
                     if not c.get("stroke") and c.get("fill") == "white"]
            self.assertGreater(len(holes), 0)
            self.assertTrue(all(float(c.attrib["cx"][:-2]) < 105 for c in holes))
            # 光のある N0 と S1 だけに、右の列の画像がある。黒地に白で、光の所が明るい。
            images = self.svg_images(svg)
            self.assertEqual(2, len(images))
            self.assertTrue(all(x > 105 and pixels.max() > 0 and pixels.min() == 0 for x, pixels in images))
            for name in ("galaxy-0.pdf", "galaxy-0-mirrored.pdf"):
                self.assertEqual(2, len(PdfReader(output / name).pages[0].images))

    def test_normal_color_images_are_black_on_white(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.column_config(directory, "galaxy.etch.render=dither\n")
            config.write_text(config.read_text().replace("color.invert=yes", "color.invert=no"))
            with patch("sys.stdout", new_callable=StringIO):
                self.assertEqual(0, main(["-svg", "-f", str(config)]))
            images = self.svg_images(Path(directory) / "output" / "galaxy" / "normal" / "galaxy-0.svg")
            self.assertEqual(2, len(images))
            for _, pixels in images:
                self.assertEqual({0, 255}, set(np.unique(pixels).tolist()))
                self.assertGreater((pixels == 255).mean(), .9)

    def test_default_outputs_svg_and_pdf(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            with patch("sys.stdout", new_callable=StringIO):
                self.assertEqual(0, main(["-f", str(config)]))
            output = Path(directory) / "output" / "galaxy" / "inverted"
            self.assertEqual(["galaxy-0.pdf", "galaxy-0.svg", "galaxy-1.pdf", "galaxy-1.svg"],
                             sorted(p.name for p in output.iterdir()))


class BundledConfigTests(unittest.TestCase):
    def test_bundled_gaia_configs_are_valid(self):
        """同梱の設定ファイルが、生成の前の検査（列の設定・見出しの大きさを含む）を通る。"""
        from gingaforge.config import PROJECT_ROOT, load_properties
        from gingaforge.gaia_etch import (_column_headings, _column_properties, _init_galaxy_transformer,
                                          _init_plate_writer)
        expected = {
            "gaia-etch.properties": [("holes", .05)],
            "gaia-ohp.properties": [("gray", .4), ("dither", .4), ("holes", .05)],
        }
        for name, columns in expected.items():
            with self.subTest(name=name):
                path = PROJECT_ROOT / name
                props = load_properties(path)
                found = _column_properties(props)
                settings = [EtchingSettings.from_properties(column) for _, column in found]
                self.assertEqual(columns, [(s.render, s.flux_gain) for s in settings])
                for _, column in found:
                    GaiaSettings.from_properties(column, path.parent)
                    _init_galaxy_transformer(column)
                _init_plate_writer(props, path, column_headings=_column_headings([change for change, _ in found]))


if __name__ == "__main__":
    unittest.main()
