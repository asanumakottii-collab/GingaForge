from collections import Counter
import csv
import gzip
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
from gingaforge.galaxy_etching import EtchingBuilder, EtchingSettings
from gingaforge.gaia_etch import GalaxyTransformer, main
from gingaforge.models import SpherePosition


def position(lon, lat=0.):
    return SpherePosition.from_galactic(lon, lat)


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
        self.t = GalaxyTransformer(.1, 6500., 0., 0., 50.)
        self.settings = EtchingSettings()

    def test_inverse_projection_with_offset(self):
        for horizontal, vertical in ((0., 0.), (355., 200.)):
            transformer = GalaxyTransformer(.1, 6500., horizontal, vertical, 50.)
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
        builder = EtchingBuilder(GalaxyTransformer(.2, 13000., 0., 0., 100.), settings, 69.25)
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
    def largest_step_between_neighbours(builder):
        """光量のあるセル（端の3列・3行を除く）で、銀経方向に隣り合うセルの光量の比の最大のずれ。"""
        cells = builder.cells
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
        area = EtchingBuilder(self.t, settings, 69.25)
        point = EtchingBuilder(self.t, settings, 69.25)
        count = self.uniform_sky((area, point), 20., -3., add_point=point)
        self.assertAlmostEqual(1., sum(area.cells.values()) / (count * area.reference_area), places=12)
        self.assertAlmostEqual(.25, self.largest_step_between_neighbours(point), places=9)
        self.assertLess(self.largest_step_between_neighbours(area), 1e-9)
        # 銀緯方向も含め、内側のセルはほぼ同じ光量（tan b による伸びの差だけ）。
        inner = [v for (_, _, ix, iy), v in area.cells.items() if -5 <= ix <= 4 and -20 <= iy <= 19]
        self.assertLess(max(inner) / min(inner), 1.01)

    def test_uniform_sky_bins_with_projector_offset(self):
        transformer = GalaxyTransformer(.1, 6500., 355., 200., 50.)
        area = EtchingBuilder(transformer, EtchingSettings(cell_mm=.4), 69.25)
        point = EtchingBuilder(transformer, EtchingSettings(cell_mm=.4), 69.25)
        count = self.uniform_sky((area, point), 80., 10., add_point=point)
        self.assertEqual({(0, 1)}, {key[:2] for key in area.cells})
        self.assertAlmostEqual(1., sum(area.cells.values()) / (count * area.reference_area), places=12)
        self.assertGreater(self.largest_step_between_neighbours(point), .2)
        self.assertLess(self.largest_step_between_neighbours(area), .005)

    def test_bin_is_split_by_overlap_area(self):
        builder = EtchingBuilder(self.t, EtchingSettings(cell_mm=.4), 69.25)
        # 原盤中心（セルの角）をまたぐビン。x = p tan b、y = -p (l - 30度)。
        builder.add_bin(29.96, 30.06, -.04, .06, 1.)
        tan = lambda degrees: math.tan(math.radians(degrees))
        x_fraction = {-1: tan(.04) / (tan(.04) + tan(.06)), 0: tan(.06) / (tan(.04) + tan(.06))}
        y_fraction = {-1: .6, 0: .4}
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
        builder = EtchingBuilder(self.t, EtchingSettings(**settings), 69.25)
        # 銀緯方向(x)に明るさが落ちる帯と、外側の明るい星1つ。
        for ix in range(-20, 20):
            for iy in range(-10, 10):
                x = (ix + .5) * builder.cell
                flux = .02 + .1 * math.exp(-(x / 1.5) ** 2) + (.4 if (ix, iy) == (-18, 0) else 0.)
                builder.add_sample(self.t.inverse_transform_unit(x, (iy + .5) * builder.cell, 0, 0), flux)
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
            circles = ET.parse(pages[0]).getroot().findall("{http://www.w3.org/2000/svg}circle")
            self.assertGreater(len(circles), 0)
            self.assertTrue(all(math.isclose(float(c.attrib["r"][:-2]), .01) for c in circles))

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

    def test_default_outputs_svg_and_pdf(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            with patch("sys.stdout", new_callable=StringIO):
                self.assertEqual(0, main(["-f", str(config)]))
            output = Path(directory) / "output" / "galaxy" / "inverted"
            self.assertEqual(["galaxy-0.pdf", "galaxy-0.svg", "galaxy-1.pdf", "galaxy-1.svg"],
                             sorted(p.name for p in output.iterdir()))


if __name__ == "__main__":
    unittest.main()
