import importlib.util
import math
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from gingaforge.gaia_etch import GalaxyTransformer, _init_galaxy_transformer, _init_plate_writer
from gingaforge.galaxy_geometry import lens_distance
from gingaforge.models import SpherePosition


def _position(longitude, latitude=0.):
    return SpherePosition.from_galactic(longitude, latitude)


class GalacticCoordinatesTests(unittest.TestCase):
    def test_canonical_galactic_pole_and_equator_node(self):
        # Hipparcosの定義角。行列の逆変換を使わずにICRS座標から検証する。
        position = SpherePosition()
        position.radeg, position.dedeg = 192.85948, 27.12825
        self.assertAlmostEqual(90., position.to_galactic()[1], places=10)
        position.radeg, position.dedeg = 282.85948, 0.
        longitude, latitude = position.to_galactic()
        self.assertAlmostEqual(32.93192, longitude, places=10)
        self.assertAlmostEqual(0., latitude, places=10)

        equatorial = _position(32.93192)
        self.assertAlmostEqual(282.85948, equatorial.radeg % 360., places=10)
        self.assertAlmostEqual(0., equatorial.dedeg, places=10)

    def test_galactic_coordinates_round_trip(self):
        for longitude in (-360., -120., 0., 60., 180., 359., 720.):
            for latitude in (-80., -20., 0., 20., 80.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    actual_l, actual_b = _position(longitude, latitude).to_galactic()
                    self.assertAlmostEqual(0., (actual_l - longitude + 180.) % 360. - 180., places=10)
                    self.assertAlmostEqual(latitude, actual_b, places=10)


class GalaxyTransformerTests(unittest.TestCase):
    def setUp(self):
        self.transformer = GalaxyTransformer(0.125, 6500., 355., 200., 50., 19.)

    def test_each_region_contains_both_galactic_latitude_signs(self):
        for longitude, unit in ((60., (0, 0)), (100., (0, 1)),
                                (20., (1, 0)), (340., (1, 1))):
            for latitude in (-19., 0., 19.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    position = _position(longitude, latitude)
                    self.assertEqual(unit, self.transformer.assigned_unit(position))
                    transformed = self.transformer.transform(position)
                    self.assertEqual(unit, (transformed.dir, transformed.index))
                    position.radeg += 360.
                    self.assertEqual(unit, self.transformer.assigned_unit(position))

    def test_longitude_boundaries_have_one_owner_and_gaps_are_excluded(self):
        for longitude, unit in (
            (0., (1, 0)), (40. - 1e-6, (1, 0)), (40., (0, 0)), (40. + 1e-6, (0, 0)),
            (80. - 1e-6, (0, 0)), (80., (0, 1)),
            (120. - 1e-6, (0, 1)), (120., None), (150., None),
            (180. - 1e-6, None), (180., None), (240., None), (300., None),
            (320. - 1e-6, None), (320., (1, 1)), (340., (1, 1)),
            (360. - 1e-6, (1, 1)), (360., (1, 0)), (-40., (1, 1)), (-40. - 1e-6, None),
        ):
            with self.subTest(longitude=longitude):
                self.assertEqual(unit, self.transformer.assigned_unit(_position(longitude)))

    def test_latitude_limits_are_inclusive_and_outside_is_excluded(self):
        for longitude in (20., 60., 100., 340.):
            for latitude in (-20., 20.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    self.assertIsNotNone(self.transformer.transform(_position(longitude, latitude)))
            for latitude in (-90., -20.000001, 20.000001, 90.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    self.assertIsNone(self.transformer.transform(_position(longitude, latitude)))

    def test_each_region_center_projects_to_plate_center_with_projector_offset(self):
        for longitude in (60., 100., 20., 340.):
            with self.subTest(longitude=longitude):
                position = self.transformer.transform(_position(longitude))
                self.assertAlmostEqual(0., position.xmm, places=10)
                self.assertAlmostEqual(0., position.ymm, places=10)

    def test_lens_distance_follows_lens_formula(self):
        for image, focus in ((6450., 19.), (6450., 50.), (100., 25.), (6.5e12, 50.)):
            with self.subTest(image=image, focus=focus):
                lens = lens_distance(image, focus)
                self.assertAlmostEqual(1., focus * (1 / lens + 1 / (image - lens)), places=12)
                self.assertLessEqual(lens, image / 2)
        for image, focus in ((99., 25.), (6450., 0.), (6450., -1.)):
            with self.subTest(image=image, focus=focus), self.assertRaisesRegex(ValueError, "焦点距離"):
                lens_distance(image, focus)

    def test_flat_plate_projection_scale_and_direction(self):
        # 中心配置。レンズ中心は原盤中心の方向に plate + lens だけ出ていて、原板は平ら。
        transformer = GalaxyTransformer(0.125, 6500., 0., 0., 50., 19.)
        lens = lens_distance(6500. - 50., 19.)
        offset = 50. + lens
        for center in (60., 100., 20., 340.):
            self.assertAlmostEqual(lens, transformer.lens_distance(*transformer.unit_of_galactic(center, 0.)), places=12)
            for delta_l, latitude in ((-10., -15.), (10., 15.), (30., 20.)):
                with self.subTest(center=center, delta_l=delta_l):
                    position = transformer.transform_unit(_position(center + delta_l, latitude),
                                                          *transformer.unit_of_galactic(center, 0.))
                    b, dl = math.radians(latitude), math.radians(delta_l)
                    forward = 6500. * math.cos(b) * math.cos(dl) - offset
                    east, north = 6500. * math.cos(b) * math.sin(dl), 6500. * math.sin(b)
                    self.assertAlmostEqual(lens * north / forward, position.xmm, places=10)
                    self.assertAlmostEqual(-lens * east / forward, position.ymm, places=10)

    def test_default_focal_length_gives_about_14_mm_plates(self):
        transformer = GalaxyTransformer(0.125, 6500., 0., 0., 50., 19.)
        self.assertAlmostEqual(19.0563, transformer.lens_distance(0, 0), places=4)
        # N0（銀経40〜80度）。銀経40度幅・銀緯40度幅の正方形の空。
        top, bottom = (transformer.transform_unit(_position(l), 0, 0) for l in (40., 80.))
        side = transformer.transform_unit(_position(60., 20.), 0, 0)
        corner = transformer.transform_unit(_position(40., 20.), 0, 0)
        self.assertAlmostEqual(14.03, top.ymm - bottom.ymm, places=2)
        # 平らな原板では、銀緯±20度の境界は両端ほど広がる。
        self.assertAlmostEqual(14.03, 2 * side.xmm, places=2)
        self.assertAlmostEqual(14.94, 2 * corner.xmm, places=2)
        self.assertLess(math.hypot(corner.xmm, corner.ymm), 15.)  # gaia-ohp の直径30 mmの枠に収まる

    def test_scale_shrinks_whole_geometry(self):
        for horizontal, vertical in ((0., 0.), (355., 200.)):
            full = GalaxyTransformer(0.125, 6500., horizontal, vertical, 50., 19.)
            scale = .25
            small = GalaxyTransformer(0.125, 6500. * scale, horizontal * scale, vertical * scale, 50. * scale, 19. * scale)
            for longitude, latitude in ((10., -15.), (95., 18.), (65., 3.), (333., -19.)):
                with self.subTest(horizontal=horizontal, longitude=longitude):
                    a, b = full.transform(_position(longitude, latitude)), small.transform(_position(longitude, latitude))
                    self.assertAlmostEqual(a.xmm * scale, b.xmm, places=10)
                    self.assertAlmostEqual(a.ymm * scale, b.ymm, places=10)

    def test_config_focal_length_is_scaled_like_otl(self):
        props = {"scale": ".5", "dome-radius": "6500", "plate-distance": "50", "focal-length": "19",
                 "projector-horizontal": "0", "projector-vertical": "0"}
        transformer = _init_galaxy_transformer(props)
        self.assertEqual(9.5, transformer.focus)
        self.assertAlmostEqual(.5 * lens_distance(6450., 19.), transformer.lens_distance(0, 0), places=12)
        with self.assertRaisesRegex(ValueError, "焦点距離"):
            _init_galaxy_transformer({**props, "focal-length": "2000"})

    def test_entire_band_fits_default_frames(self):
        for longitude in range(360):
            for latitude in (-20., 0., 20.):
                position = self.transformer.transform(_position(longitude, latitude))
                if position is not None:
                    with self.subTest(longitude=longitude, latitude=latitude):
                        self.assertTrue(math.isfinite(position.xmm))
                        self.assertTrue(math.isfinite(position.ymm))
                        self.assertLess(abs(position.xmm), 69.25)
                        self.assertLess(abs(position.ymm), 69.25)


class GalaxyOutputTests(unittest.TestCase):
    def test_configured_frame_diameter_in_svg_and_pdf(self):
        for pdf in (False, True):
            if pdf and not importlib.util.find_spec("reportlab"):
                continue
            # 既定は全4原盤を1列に並べるので、1枚の上限は 297 / 4 = 74.25 mm。自動の直径はそれより10 mm小さい。
            for configured, expected in ((None, 32.125), ("0", 32.125), ("60", 30.)):
                with self.subTest(pdf=pdf, configured=configured):
                    with tempfile.TemporaryDirectory() as directory:
                        props = {"scale": "2"}
                        if configured is not None:
                            props["plate.frame-diameter"] = configured
                        with _init_plate_writer(props, Path(directory, "gaia-etch.properties"), pdf) as writer:
                            self.assertEqual(expected, writer.r)
                            writer.write_frames()
                        if not pdf:
                            for page in Path(directory, "output", "gaia-etch", "normal").glob("*.svg"):
                                for frame in ET.parse(page).getroot().findall("{http://www.w3.org/2000/svg}circle"):
                                    self.assertEqual(f"{expected}mm", frame.get("r"))

    def test_configured_frame_diameter_draws_circle(self):
        with tempfile.TemporaryDirectory() as directory:
            props = {"plate.frame-diameter": "30", "color.invert": "yes"}
            with _init_plate_writer(props, Path(directory, "gaia-ohp.properties")) as writer:
                self.assertEqual(15., writer.r)
                writer.write_frames()
                # 円の内側、外接正方形の角（円の外側）、円の外側。
                writer.write_hole(0, 0, 10., -10., .01)
                writer.write_hole(0, 0, 14., 14., .01)
                writer.write_hole(0, 0, 15.01, 0., .01)
            root = ET.parse(Path(directory, "output", "gaia-ohp", "inverted", "gaia-ohp-0.svg")).getroot()
            namespace = "{http://www.w3.org/2000/svg}"
            self.assertEqual([], root.findall(namespace + "rect"))
            circles = root.findall(namespace + "circle")
            frames = [c for c in circles if c.get("stroke")]
            self.assertEqual(["15.0mm"] * 4, [c.get("r") for c in frames])
            self.assertEqual(1, len(circles) - len(frames))

    def test_config_name_and_color_invert_select_output(self):
        for value, folder in (("no", "normal"), ("yes", "inverted"), ("y", "inverted")):
            with self.subTest(value=value):
                writer = _init_plate_writer({"color.invert": value}, Path("configs", "gaia-ohp.properties"))
                self.assertEqual(Path("configs", "output", "gaia-ohp", folder), writer.output_dir)
                self.assertEqual("gaia-ohp-", writer.filename_prefix)

    def test_invalid_frame_diameter_is_rejected(self):
        for value in ("-1", "nan", "inf"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "plate.frame-diameter"):
                    _init_plate_writer({"plate.frame-diameter": value}, "gaia-etch.properties")

    def test_holes_outside_frame_are_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            with _init_plate_writer({}, Path(directory, "gaia-etch.properties")) as writer:
                writer.write_frames()
                writer.write_hole(0, 0, 0., -writer.r, .01)
                writer.write_hole(0, 0, writer.r + .01, 0., .01)
                writer.write_hole(1, 1, 0., math.nan, .01)
            holes = [c for page in Path(directory, "output", "gaia-etch", "normal").glob("*.svg")
                     for c in ET.parse(page).getroot().findall("{http://www.w3.org/2000/svg}circle") if not c.get("stroke")]
            self.assertEqual(1, len(holes))

    def test_svg_creates_all_four_frames_without_stars(self):
        with tempfile.TemporaryDirectory() as directory:
            with _init_plate_writer({}, Path(directory, "gaia-etch.properties")) as writer:
                writer.write_frames()
            pages = sorted(Path(directory, "output", "gaia-etch", "normal").glob("*.svg"))
            self.assertEqual(["gaia-etch-0.svg"], [p.name for p in pages])
            root = ET.parse(pages[0]).getroot()
            self.assertEqual(4, len(root.findall("{http://www.w3.org/2000/svg}circle")))
            # 既定は空でつながる順に、上から1枚ずつ。
            texts = root.findall("{http://www.w3.org/2000/svg}text")
            self.assertEqual(["N1", "N0", "S0", "S1"], [e.text for e in texts])
            self.assertEqual([37.125 - 32.125, 111.375 - 32.125, 185.625 - 32.125, 259.875 - 32.125],
                             [float(e.get("y")[:-2]) for e in texts])

    def test_select_order_and_per_column_split_pages(self):
        with tempfile.TemporaryDirectory() as directory:
            props = {"plate.select": "S0, N1, N0", "plate.per-column": "2", "plate.frame-diameter": "40"}
            with _init_plate_writer(props, Path(directory, "gaia-etch.properties")) as writer:
                self.assertEqual(148.5, writer.slot_height)
                writer.write_hole(1, 0, 0., 0., .5)
            pages = sorted(Path(directory, "output", "gaia-etch", "normal").glob("*.svg"))
            self.assertEqual(["gaia-etch-0.svg", "gaia-etch-1.svg"], [p.name for p in pages])
            labels = [[e.text for e in ET.parse(page).getroot().findall("{http://www.w3.org/2000/svg}text")]
                      for page in pages]
            self.assertEqual([["S0", "N1"], ["N0"]], labels)
            holes = [c for c in ET.parse(pages[0]).getroot().findall("{http://www.w3.org/2000/svg}circle")
                     if not c.get("stroke")]
            self.assertEqual([("105.0mm", "74.25mm")], [(c.get("cx"), c.get("cy")) for c in holes])

    def test_repeat_copies_the_order(self):
        with tempfile.TemporaryDirectory() as directory:
            props = {"plate.select": "N1, N0", "plate.repeat": "3", "plate.frame-diameter": "40"}
            with _init_plate_writer(props, Path(directory, "gaia-etch.properties")) as writer:
                writer.write_hole(0, 0, 0., 0., .5)
            root = ET.parse(Path(directory, "output", "gaia-etch", "normal", "gaia-etch-0.svg")).getroot()
            self.assertEqual(["N1", "N0"] * 3, [e.text for e in root.findall("{http://www.w3.org/2000/svg}text")])
            holes = [float(c.get("cy")[:-2]) for c in root.findall("{http://www.w3.org/2000/svg}circle")
                     if not c.get("stroke")]
            self.assertEqual([74.25, 173.25, 272.25], holes)

    def test_removed_layout_keys_are_rejected(self):
        for key in ("plate.column", "plate.row", "row.count"):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, f"{key} は廃止"):
                _init_plate_writer({key: "1"}, "gaia-etch.properties")
        for props, message in (({"plate.select": "N0, N0"}, "2回"), ({"plate.select": "N0, X1"}, "plate.select"),
                               ({"plate.per-column": "-1"}, "plate.per-column"),
                               ({"plate.repeat": "0"}, "plate.repeat")):
            with self.subTest(props=props), self.assertRaisesRegex(ValueError, message):
                _init_plate_writer(props, "gaia-etch.properties")

    @unittest.skipUnless(importlib.util.find_spec("reportlab") and importlib.util.find_spec("pypdf"),
                         "ReportLab and pypdf are required for PDF tests")
    def test_pdf_creates_all_four_frames_without_stars(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory:
            with _init_plate_writer({}, Path(directory, "gaia-etch.properties"), pdf=True) as writer:
                writer.write_frames()
            pages = sorted(Path(directory, "output", "gaia-etch", "normal").glob("*.pdf"))
            self.assertEqual(1, len(pages))
            pdf = PdfReader(pages[0])
            self.assertEqual(1, len(pdf.pages))
            self.assertEqual(["N1", "N0", "S0", "S1"], pdf.pages[0].extract_text().split())

    @unittest.skipUnless(importlib.util.find_spec("reportlab") and importlib.util.find_spec("pypdf"),
                         "ReportLab and pypdf are required for PDF tests")
    def test_mirrored_pdf_adds_flipped_second_page_to_same_file(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory:
            with _init_plate_writer({}, Path(directory, "gaia-etch.properties"), pdf=True, mirrored=True) as writer:
                writer.write_hole(0, 0, 10., 0., .5)
            folder = Path(directory, "output", "gaia-etch", "normal")
            self.assertEqual(["gaia-etch-0.pdf"], [p.name for p in folder.glob("*.pdf")])
            normal, mirrored = PdfReader(folder / "gaia-etch-0.pdf").pages
            self.assertEqual(["N1", "N0", "S0", "S1"], normal.extract_text().split())
            front = normal.get_contents().get_data().decode()
            back = mirrored.get_contents().get_data().decode()
            # 両ページは同じ描画内容を呼び出し、2ページ目にだけ x → 用紙幅 - x の変換が加わる。
            draw = re.search(r"/\S+ Do", front).group()
            flip = re.search(r"-1 0 0 1 ([\d.]+) 0 cm\n", back)
            self.assertAlmostEqual(float(mirrored.mediabox.width), float(flip.group(1)), places=3)
            self.assertNotIn("-1 0 0 1", front)
            self.assertIn(draw, back)
            self.assertLess(back.index(flip.group()), back.index(draw))


if __name__ == "__main__":
    unittest.main()
