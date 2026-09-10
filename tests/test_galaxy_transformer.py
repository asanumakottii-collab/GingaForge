import importlib.util
import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from galaxy_transformer import GalaxyTransformer, _init_plate_writer, main
from models import SpherePosition
from plate_writer import PlateWriterType


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
    def test_configured_frame_size_in_svg_and_pdf(self):
        for writer_type in (PlateWriterType.SVG, PlateWriterType.PDF):
            if writer_type == PlateWriterType.PDF and not importlib.util.find_spec("reportlab"):
                continue
            for configured, expected in ((None, 69.25), ("0", 69.25), ("100", 50.)):
                with self.subTest(writer_type=writer_type, configured=configured):
                    with tempfile.TemporaryDirectory() as directory, patch("sys.stdout", new_callable=StringIO):
                        props = {"output.directory": directory, "scale": "2"}
                        if configured is not None:
                            props["plate.frame-size"] = configured
                        writer = _init_plate_writer(props, writer_type)
                        self.assertEqual(expected, writer.r)
                        self.assertTrue(writer._is_position_in_frame(expected, expected))
                        self.assertFalse(writer._is_position_in_frame(expected + .01, 0.))
                        writer.close()
                        if writer_type == PlateWriterType.SVG:
                            for page in Path(directory, "galaxy").glob("*.svg"):
                                for rect in ET.parse(page).getroot().findall("{http://www.w3.org/2000/svg}rect"):
                                    self.assertEqual(f"{expected * 2}mm", rect.get("width"))
                                    self.assertEqual(f"{expected * 2}mm", rect.get("height"))

    def test_invalid_frame_size_is_rejected(self):
        for value in ("-1", "nan", "inf"):
            with self.subTest(value=value), patch("sys.stdout", new_callable=StringIO):
                with self.assertRaisesRegex(ValueError, "plate.frame-size"):
                    _init_plate_writer({"plate.frame-size": value}, PlateWriterType.SVG)

    def test_svg_creates_all_four_frames_without_stars(self):
        with tempfile.TemporaryDirectory() as directory, patch("sys.stdout", new_callable=StringIO):
            writer = _init_plate_writer({"output.directory": directory}, PlateWriterType.SVG)
            writer.close()
            pages = sorted(Path(directory, "galaxy").glob("*.svg"))
            self.assertEqual(["gaia-etch-0.svg", "gaia-etch-1.svg"], [p.name for p in pages])
            for page, labels in zip(pages, (["N0", "S0"], ["N1", "S1"])):
                root = ET.parse(page).getroot()
                self.assertEqual(2, len(root.findall("{http://www.w3.org/2000/svg}rect")))
                self.assertEqual(labels, [e.text for e in root.findall("{http://www.w3.org/2000/svg}text")])

    @unittest.skipUnless(importlib.util.find_spec("reportlab") and importlib.util.find_spec("pypdf"),
                         "ReportLab and pypdf are required for PDF tests")
    def test_pdf_creates_all_four_frames_without_stars(self):
        from pypdf import PdfReader

        with tempfile.TemporaryDirectory() as directory, patch("sys.stdout", new_callable=StringIO):
            writer = _init_plate_writer({"output.directory": directory}, PlateWriterType.PDF)
            writer.close()
            pages = sorted(Path(directory, "galaxy").glob("*.pdf"))
            self.assertEqual(2, len(pages))
            for path, labels in zip(pages, (["N0", "S0"], ["N1", "S1"])):
                pdf = PdfReader(path)
                self.assertEqual(1, len(pdf.pages))
                self.assertEqual(labels, pdf.pages[0].extract_text().split())


if __name__ == "__main__":
    unittest.main()
