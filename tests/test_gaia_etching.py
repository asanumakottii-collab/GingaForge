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

from gingaforge.gaia_catalog import (GaiaSettings, build_queries, iter_flux_samples,
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

    def test_invalid_manufacturing_settings_fail(self):
        for key, value in (("galaxy.etch.min-hole-diameter-mm", "0"),
                           ("galaxy.etch.min-web-mm", "-1"),
                           ("galaxy.etch.flux-gain", "nan"),
                           ("scale", "inf"), ("galaxy.etch.seed", "1.2"),
                           ("galaxy.etch.cell-size-mm", ".01")):
            with self.subTest(key=key), self.assertRaises(ValueError):
                EtchingSettings.from_properties({key: value})


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


class GaiaCommandTests(unittest.TestCase):
    def config_and_input(self, directory):
        directory = Path(directory)
        write_cache(directory / "gaia", {0: [(30.05, .05, .01, 1)]})
        config = directory / "galaxy.properties"
        config.write_text(f"galaxy.gaia.input=gaia\n"
                          f"output.directory={directory / 'out'}\ncolor.invert=yes\n")
        return config

    def test_real_svg_pipeline_has_no_prompts(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            with (patch("builtins.input", side_effect=AssertionError("unexpected prompt")),
                  patch("sys.stdout", new_callable=StringIO)):
                self.assertEqual(0, main(["-f", str(config)]))
            output = Path(directory) / "out" / "galaxy"
            pages = sorted(output.glob("*.svg"))
            self.assertEqual(2, len(pages))
            self.assertEqual(pages, sorted(output.iterdir()))
            circles = ET.parse(pages[0]).getroot().findall("{http://www.w3.org/2000/svg}circle")
            self.assertGreater(len(circles), 0)
            self.assertTrue(all(math.isclose(float(c.attrib["r"][:-2]), .01) for c in circles))

    def test_missing_input_does_not_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            output = Path(directory) / "out" / "galaxy"
            output.mkdir(parents=True)
            sentinel = output / "gaia-etch-0.svg"
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
            pages = sorted((Path(directory) / "out" / "galaxy").glob("*.pdf"))
            self.assertEqual(2, len(pages))
            self.assertEqual(["N0", "S0"], PdfReader(pages[0]).pages[0].extract_text().split())


if __name__ == "__main__":
    unittest.main()
