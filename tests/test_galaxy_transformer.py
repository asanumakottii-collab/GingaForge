import importlib.util
import math
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from gingaforge.gaia_etch import GalaxyTransformer, _init_plate_writer
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
        self.transformer = GalaxyTransformer(0.125, 6500., 355., 200., 50.)

    def test_each_region_contains_both_galactic_latitude_signs(self):
        for longitude, unit in ((30., (0, 0)), (90., (0, 1)),
                                (330., (1, 0)), (270., (1, 1))):
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
            (0., (0, 0)), (60. - 1e-6, (0, 0)), (60., (0, 1)),
            (120. - 1e-6, (0, 1)), (120., None), (150., None),
            (180. - 1e-6, None), (180., None),
            (240. - 1e-6, None), (240., (1, 1)), (270., (1, 1)),
            (300. - 1e-6, (1, 1)), (300., (1, 0)), (300. + 1e-6, (1, 0)),
            (360. - 1e-6, (1, 0)), (360., (0, 0)), (-60., (1, 0)),
        ):
            with self.subTest(longitude=longitude):
                self.assertEqual(unit, self.transformer.assigned_unit(_position(longitude)))

    def test_latitude_limits_are_inclusive_and_outside_is_excluded(self):
        for longitude in (30., 90., 270., 300.):
            for latitude in (-20., 20.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    self.assertIsNotNone(self.transformer.transform(_position(longitude, latitude)))
            for latitude in (-90., -20.000001, 20.000001, 90.):
                with self.subTest(longitude=longitude, latitude=latitude):
                    self.assertIsNone(self.transformer.transform(_position(longitude, latitude)))

    def test_each_region_center_projects_to_plate_center_with_projector_offset(self):
        for longitude in (30., 90., 330., 270.):
            with self.subTest(longitude=longitude):
                position = self.transformer.transform(_position(longitude))
                self.assertAlmostEqual(0., position.xmm, places=10)
                self.assertAlmostEqual(0., position.ymm, places=10)

    def test_cylindrical_projection_scale_and_direction(self):
        transformer = GalaxyTransformer(0.125, 6500., 0., 0., 50.)
        for center in (30., 90., 330., 270.):
            for delta_l, latitude in ((-10., -15.), (10., 15.)):
                with self.subTest(center=center, delta_l=delta_l):
                    position = transformer.transform(_position(center + delta_l, latitude))
                    self.assertAlmostEqual(50. * math.tan(math.radians(latitude)), position.xmm, places=10)
                    self.assertAlmostEqual(-50. * math.radians(delta_l), position.ymm, places=10)

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
            for configured, expected in ((None, 69.25), ("0", 69.25), ("100", 50.)):
                with self.subTest(pdf=pdf, configured=configured):
                    with tempfile.TemporaryDirectory() as directory:
                        props = {"scale": "2"}
                        if configured is not None:
                            props["plate.frame-diameter"] = configured
                        with _init_plate_writer(props, Path(directory, "gaia-etch.properties"), pdf) as writer:
                            self.assertEqual(expected, writer.r)
                            writer.write_frames(2)
                        if not pdf:
                            for page in Path(directory, "output", "gaia-etch", "normal").glob("*.svg"):
                                for frame in ET.parse(page).getroot().findall("{http://www.w3.org/2000/svg}circle"):
                                    self.assertEqual(f"{expected}mm", frame.get("r"))

    def test_configured_frame_diameter_draws_circle(self):
        with tempfile.TemporaryDirectory() as directory:
            props = {"plate.frame-diameter": "30", "color.invert": "yes"}
            with _init_plate_writer(props, Path(directory, "gaia-ohp.properties")) as writer:
                self.assertEqual(15., writer.r)
                writer.write_frames(2)
                # 円の内側、外接正方形の角（円の外側）、円の外側。
                writer.write_hole(0, 0, 10., -10., .01)
                writer.write_hole(0, 0, 14., 14., .01)
                writer.write_hole(0, 0, 15.01, 0., .01)
            root = ET.parse(Path(directory, "output", "gaia-ohp", "inverted", "gaia-ohp-0.svg")).getroot()
            namespace = "{http://www.w3.org/2000/svg}"
            self.assertEqual([], root.findall(namespace + "rect"))
            circles = root.findall(namespace + "circle")
            frames = [c for c in circles if c.get("stroke")]
            self.assertEqual(["15.0mm"] * 2, [c.get("r") for c in frames])
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
                writer.write_frames(2)
                writer.write_hole(0, 0, 0., -writer.r, .01)
                writer.write_hole(0, 0, writer.r + .01, 0., .01)
                writer.write_hole(1, 1, 0., math.nan, .01)
            holes = [c for page in Path(directory, "output", "gaia-etch", "normal").glob("*.svg")
                     for c in ET.parse(page).getroot().findall("{http://www.w3.org/2000/svg}circle") if not c.get("stroke")]
            self.assertEqual(1, len(holes))

    def test_svg_creates_all_four_frames_without_stars(self):
        with tempfile.TemporaryDirectory() as directory:
            with _init_plate_writer({}, Path(directory, "gaia-etch.properties")) as writer:
                writer.write_frames(2)
            pages = sorted(Path(directory, "output", "gaia-etch", "normal").glob("*.svg"))
            self.assertEqual(["gaia-etch-0.svg", "gaia-etch-1.svg"], [p.name for p in pages])
            for page, labels in zip(pages, (["N0", "S0"], ["N1", "S1"])):
                root = ET.parse(page).getroot()
                self.assertEqual(2, len(root.findall("{http://www.w3.org/2000/svg}circle")))
                self.assertEqual(labels, [e.text for e in root.findall("{http://www.w3.org/2000/svg}text")])

    @unittest.skipUnless(importlib.util.find_spec("reportlab") and importlib.util.find_spec("pypdf"),
                         "ReportLab and pypdf are required for PDF tests")
    def test_pdf_creates_all_four_frames_without_stars(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory:
            with _init_plate_writer({}, Path(directory, "gaia-etch.properties"), pdf=True) as writer:
                writer.write_frames(2)
            pages = sorted(Path(directory, "output", "gaia-etch", "normal").glob("*.pdf"))
            self.assertEqual(2, len(pages))
            for path, labels in zip(pages, (["N0", "S0"], ["N1", "S1"])):
                pdf = PdfReader(path)
                self.assertEqual(1, len(pdf.pages))
                self.assertEqual(labels, pdf.pages[0].extract_text().split())

    @unittest.skipUnless(importlib.util.find_spec("reportlab") and importlib.util.find_spec("pypdf"),
                         "ReportLab and pypdf are required for PDF tests")
    def test_mirrored_pdf_is_flipped_horizontally_into_separate_file(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory:
            for mirrored in (False, True):
                with _init_plate_writer({}, Path(directory, "gaia-etch.properties"), pdf=True,
                                        mirrored=mirrored) as writer:
                    writer.write_hole(0, 0, 10., 0., .5)
            folder = Path(directory, "output", "gaia-etch", "normal")
            self.assertEqual(["gaia-etch-0-mirrored.pdf", "gaia-etch-0.pdf"],
                             sorted(p.name for p in folder.glob("*.pdf")))
            normal, mirrored = (PdfReader(folder / name).pages[0]
                                for name in ("gaia-etch-0.pdf", "gaia-etch-0-mirrored.pdf"))
            content = mirrored.get_contents().get_data().decode()
            # 描画内容は通常版と同じで、ページ全体に x → 用紙幅 - x の変換だけが加わる。
            flip = re.search(r"-1 0 0 1 ([\d.]+) 0 cm\n", content)
            self.assertAlmostEqual(float(mirrored.mediabox.width), float(flip.group(1)), places=3)
            self.assertEqual(normal.get_contents().get_data().decode(), content.replace(flip.group(0), "", 1))


if __name__ == "__main__":
    unittest.main()
