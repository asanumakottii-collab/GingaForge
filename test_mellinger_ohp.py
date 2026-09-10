import json
import unittest
import numpy as np
from astropy.io import fits
from PIL import Image
from pypdf import PdfReader
from build_mellinger_ohp import (ROOT, REGIONS, read_config, geometry, plate_to_sky,
                               celestial_wcs, bilinear, tone_map)


class MasterChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = read_config(ROOT / 'galaxyconfig.properties')
        cls.wcs = celestial_wcs(fits.getheader(ROOT / cls.cfg['source-fits']))

    def test_requested_regions_and_default_dimensions(self):
        self.assertEqual(REGIONS, {'N0': (0, 60), 'N1': (60, 120),
                                   'S0': (300, 360), 'S1': (240, 300)})
        _, _, _, _, width, height = geometry(self.cfg)
        self.assertAlmostEqual(width, 36.39702342662024)
        self.assertAlmostEqual(height, 52.35987755982988)

    def test_boundary_orientation_and_source_wcs(self):
        p, bmax, scale, _, width, height = geometry(self.cfg)
        for start, end in REGIONS.values():
            l, b = plate_to_sky(np.array([-width/2, 0, width/2]),
                                np.array([-height/2, 0, height/2]),
                                (start+end)/2, p, scale)
            np.testing.assert_allclose(l, [end % 360, (start+end)/2, start], atol=1e-10)
            np.testing.assert_allclose(b, [-bmax, 0, bmax], atol=1e-10)
            x, y = self.wcs.world_to_pixel_values(l, b)
            actual_l, actual_b = self.wcs.pixel_to_world_values(x, y)
            np.testing.assert_allclose((actual_l-l+180) % 360 - 180, 0, atol=1e-8)
            np.testing.assert_allclose(actual_b, b, atol=1e-8)
        # The source header puts Galactic center at this pixel (zero based).
        x0, y0 = self.wcs.world_to_pixel_values(0, 0)
        self.assertAlmostEqual(x0, 1799.950026799, places=6)
        self.assertAlmostEqual(y0, 899.950013387949, places=6)

    def test_saved_pixels_and_shared_tone_mapping(self):
        data = fits.getdata(ROOT / self.cfg['source-fits'])[1]
        p, _, scale, _, width, height = geometry(self.cfg)
        for name, (start, end) in REGIONS.items():
            linear = np.load(ROOT / f'output/{name}-linear-G.npy')
            rendered = np.array(Image.open(ROOT / f'output/{name}-ohp.png'))
            self.assertEqual(linear.shape, rendered.shape)
            self.assertTrue(np.all(np.isfinite(linear)))
            np.testing.assert_array_equal(tone_map(linear, self.cfg), rendered)
            ny, nx = linear.shape
            for row, col in [(0, 0), (ny//2, nx//2), (ny-1, nx-1), (ny//3, nx//4)]:
                x, y = (col+.5)/nx*width-width/2, (row+.5)/ny*height-height/2
                l, b = plate_to_sky(x, y, (start+end)/2, p, scale)
                sx, sy = self.wcs.world_to_pixel_values(l, b)
                expected = bilinear(data, sx, sy)
                self.assertAlmostEqual(float(linear[row, col])/float(expected), 1, places=6)

    def test_pdf_pages_labels_and_image_dimensions(self):
        reader = PdfReader(ROOT / 'output/mellinger-ohp-plates.pdf')
        self.assertEqual(len(reader.pages), 2)
        _, _, _, _, width, height = geometry(self.cfg)
        for page, names in zip(reader.pages, [('N0', 'S0'), ('N1', 'S1')]):
            self.assertAlmostEqual(float(page.mediabox.width), 210/25.4*72, places=3)
            self.assertAlmostEqual(float(page.mediabox.height), 297/25.4*72, places=3)
            for name in names:
                self.assertIn(name, page.extract_text())
            from pypdf.generic import ContentStream
            operations = ContentStream(page['/Contents'], reader).operations
            transforms = [values for values, op in operations if op == b'cm' and float(values[0]) > 10]
            self.assertEqual(len(transforms), 2)
            for values in transforms:
                self.assertAlmostEqual(float(values[0]), width/25.4*72, places=3)
                self.assertAlmostEqual(float(values[3]), height/25.4*72, places=3)
        self.assertEqual(len(PdfReader(ROOT / 'output/ohp-calibration.pdf').pages), 1)
        provenance = json.loads((ROOT / 'output/provenance.json').read_text())
        self.assertEqual(set(provenance['regions']), set(REGIONS))
        self.assertEqual(provenance['config'], dict(self.cfg))


if __name__ == '__main__':
    unittest.main()
