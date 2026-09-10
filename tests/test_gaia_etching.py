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

from gaia_catalog import (GaiaSettings, build_queries, iter_flux_samples,
                          _votable_to_csv, file_sha256)
from galaxy_etching import EtchingBuilder, EtchingSettings
from galaxy_transformer import GalaxyTransformer, main
from models import SpherePosition


def position(lon, lat=0.):
    return SpherePosition.from_galactic(lon, lat)


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
        holes, report = builder.build()
        self.assertEqual(1, len(holes))
        self.assertAlmostEqual(.02, 2 * holes[0].rmm)
        self.assertEqual(16, report["plates"][0]["source_stars"])
        self.assertAlmostEqual(0., report["plates"][0]["quantization_error_area_mm2"])

    def test_flux_sum_and_spacing_across_adjacent_cells(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        for ix in range(-2, 3):
            for iy in range(-2, 3):
                sp = self.t.inverse_transform_unit((ix + .5) * builder.cell,
                                                   (iy + .5) * builder.cell, 0, 0)
                builder.add_sample(sp, .023, 20)
        holes, report = builder.build()
        self.assertEqual(58, len(holes))  # 25 * 2.3 最小穴相当を四捨五入
        for i, first in enumerate(holes):
            for second in holes[i + 1:]:
                distance = math.hypot(first.p.xmm - second.p.xmm, first.p.ymm - second.p.ymm)
                self.assertGreaterEqual(distance - first.rmm - second.rmm + 1e-12, .02)
        plate = report["plates"][0]
        self.assertLessEqual(abs(plate["quantization_error_area_mm2"]), builder.quantum_area / 2 + 1e-12)
        self.assertAlmostEqual(0., plate["capacity_loss_area_mm2"])

    def test_saturation_is_reported_and_not_spilled_to_dark_cells(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        builder.add_sample(position(30.01, .01), 100., 100000)
        holes, report = builder.build()
        self.assertEqual(25, len(holes))
        self.assertEqual(1, report["plates"][0]["saturated_cells"])
        self.assertGreater(report["plates"][0]["capacity_loss_area_mm2"], 0)
        self.assertEqual(1, len({(math.floor(h.p.xmm / builder.cell),
                                 math.floor(h.p.ymm / builder.cell)) for h in holes}))

    def test_determinism_and_parameter_changes(self):
        def run(props):
            settings = EtchingSettings.from_properties(props)
            builder = EtchingBuilder(self.t, settings, 69.25)
            builder.add_sample(position(30.01, .01), .2)
            return builder.build()
        first, _ = run({})
        second, _ = run({})
        self.assertEqual([(h.p.xmm, h.p.ymm, h.rmm) for h in first],
                         [(h.p.xmm, h.p.ymm, h.rmm) for h in second])
        wider, report = run({"galaxy.etch.min-hole-diameter-mm": ".04"})
        self.assertEqual(5, len(wider))
        self.assertTrue(all(math.isclose(h.rmm, .02) for h in wider))

    def test_final_dimensions_follow_scale(self):
        settings = EtchingSettings.from_properties({"scale": "2"})
        builder = EtchingBuilder(GalaxyTransformer(.2, 13000., 0., 0., 100.), settings, 69.25)
        builder.add_sample(position(30.01, .01), .1)
        holes, report = builder.build()
        self.assertEqual(10, len(holes))
        self.assertTrue(all(math.isclose(h.rmm / settings.scale, .01) for h in holes))
        self.assertAlmostEqual(.04, report["min_web_output_mm"])

    def test_outside_regions_and_physical_frame(self):
        builder = EtchingBuilder(self.t, self.settings, .1)
        builder.add_sample(position(150.), 1., 100)
        builder.add_sample(position(30., 21.), 1., 50)
        builder.add_sample(position(30., 10.), 1., 40)
        builder.add_sample(position(30.01, .01), .1, 20)
        holes, report = builder.build()
        self.assertEqual(150, report["outside_region_stars"])
        self.assertGreater(report["plates"][0]["capacity_loss_area_mm2"], 0)
        self.assertTrue(all(max(abs(h.p.xmm), abs(h.p.ymm)) + h.rmm + .02 <= .1 for h in holes))

    def test_four_units_do_not_share_rounding_or_flux(self):
        builder = EtchingBuilder(self.t, self.settings, 69.25)
        for lon in (30.01, 90.01, 270.01, 330.01):
            builder.add_sample(position(lon, .01), .006, 1)
        holes, report = builder.build()
        self.assertEqual(4, len(holes))
        self.assertEqual([1] * 4, [p["holes"] for p in report["plates"]])

    def test_invalid_manufacturing_settings_fail(self):
        for key, value in (("galaxy.etch.min-hole-diameter-mm", "0"),
                           ("galaxy.etch.min-web-mm", "-1"),
                           ("galaxy.etch.flux-gain", "nan"),
                           ("scale", "inf"), ("galaxy.etch.seed", "1.2"),
                           ("galaxy.etch.cell-size-mm", ".01")):
            with self.subTest(key=key), self.assertRaises(ValueError):
                EtchingSettings.from_properties({key: value})


class GaiaInputTests(unittest.TestCase):
    def write_stars(self, path, entries):
        with open(path, "w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["source_id", "ra", "dec", "phot_g_mean_mag"])
            writer.writerows(entries)

    def test_magnitude_filter_and_g_flux(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "stars.csv"
            self.write_stars(path, [(1, 10, 20, 7.5), (2, 10, 20, 12.5),
                                    (3, 10, 20, 17), (4, 10, 20, 18)])
            stats = {}
            samples = list(iter_flux_samples(GaiaSettings(path), ((0, 60),), 20, stats))
            self.assertEqual(2, len(samples))
            self.assertAlmostEqual(.01, samples[0][1])
            self.assertEqual(2, stats["filtered_stars"])

    def test_duplicate_nonfinite_and_bad_headers_are_rejected(self):
        for entries in ([(1, 10, 20, 10), (1, 10, 20, 10)],
                        [(1, 10, 20, "nan")], [(1, 360, 20, 10)]):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "stars.csv"
                self.write_stars(path, entries)
                with self.assertRaises(ValueError):
                    list(iter_flux_samples(GaiaSettings(path), ((0, 60),), 20, {}))

    def test_gzip_input_and_config_relative_path(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "stars.csv.gz"
            with gzip.open(path, "wt") as stream:
                stream.write("source_id,ra,dec,phot_g_mean_mag\n1,10,20,12.5\n")
            settings = GaiaSettings.from_properties({"galaxy.gaia.input": "stars.csv.gz"}, directory)
            samples = list(iter_flux_samples(settings, ((0, 60),), 20, {}))
            self.assertEqual(1, len(samples))

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
            stats = {}
            samples = list(iter_flux_samples(settings, ((0, 60),), 20, stats))
            self.assertEqual(16, stats["stars"])
            self.assertAlmostEqual(30.05, samples[0][0].to_galactic()[0])
            with self.assertRaises(ValueError):
                list(iter_flux_samples(GaiaSettings(directory, faint=18), ((0, 60),), 20, {}))
            target.write_bytes(b"broken")
            with self.assertRaises(ValueError):
                list(iter_flux_samples(settings, ((0, 60),), 20, {}))


class GaiaCommandTests(unittest.TestCase):
    def config_and_input(self, directory):
        directory = Path(directory)
        path = directory / "stars.csv"
        sp = position(30.01, .01)
        path.write_text(f"source_id,ra,dec,phot_g_mean_mag\n1,{sp.radeg % 360},{sp.dedeg},10\n")
        config = directory / "galaxy.properties"
        config.write_text(f"galaxy.mode=gaia-etch\ngalaxy.gaia.input=stars.csv\n"
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
            circles = ET.parse(pages[0]).getroot().findall("{http://www.w3.org/2000/svg}circle")
            self.assertGreater(len(circles), 0)
            self.assertTrue(all(math.isclose(float(c.attrib["r"][:-2]), .01) for c in circles))
            report = json.loads((output / "gaia-etch-report.json").read_text())
            self.assertEqual(1, report["input"]["stars"])
            self.assertEqual(len(circles), report["plates"][0]["holes"])

    def test_missing_input_does_not_overwrite_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            config = self.config_and_input(directory)
            output = Path(directory) / "out" / "galaxy"
            output.mkdir(parents=True)
            sentinel = output / "gaia-etch-0.svg"
            sentinel.write_text("existing plate")
            (Path(directory) / "stars.csv").unlink()
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
