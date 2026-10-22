from pathlib import Path
import tempfile
import unittest

import numpy as np
from astropy.io import fits
from pypdf import PdfReader
from pypdf.generic import ContentStream

from gingaforge.config import PROJECT_ROOT as ROOT
from gingaforge.ohp import (REGIONS, read_config, geometry, plate_to_sky,
                            celestial_wcs, bilinear, reproject, tone_map, build)


class MasterChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = read_config(ROOT / 'mellinger-ohp.properties')
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

    def test_reprojection_samples_pixel_centers(self):
        cfg = read_config(ROOT / 'mellinger-ohp.properties')
        cfg['raster-dpi'] = '72'
        data = fits.getdata(ROOT / cfg['source-fits'])[1]
        p, _, scale, _, width, height = geometry(cfg)
        for start, end in REGIONS.values():
            linear = reproject(data, self.wcs, start, end, cfg)
            self.assertTrue(np.all(np.isfinite(linear)))
            ny, nx = linear.shape
            for row, col in [(0, 0), (ny//2, nx//2), (ny-1, nx-1), (ny//3, nx//4)]:
                x, y = (col+.5)/nx*width-width/2, (row+.5)/ny*height-height/2
                l, b = plate_to_sky(x, y, (start+end)/2, p, scale)
                sx, sy = self.wcs.world_to_pixel_values(l, b)
                self.assertAlmostEqual(float(linear[row, col])/float(bilinear(data, sx, sy)), 1, places=6)
            pixels = tone_map(linear, cfg)
            self.assertEqual(np.uint8, pixels.dtype)
            self.assertEqual(linear.shape, pixels.shape)

    def test_pdf_pages_labels_and_image_dimensions(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'ohp.properties'
            text = (ROOT / 'mellinger-ohp.properties').read_text()
            text = text.replace('source/mwpan2_RGB_3600.fits', str(ROOT / 'source/mwpan2_RGB_3600.fits'))
            config.write_text(text.replace('raster-dpi = 1200', 'raster-dpi = 72'))
            reader = PdfReader(build(config))
            self.assertEqual(len(reader.pages), 2)
            _, _, _, _, width, height = geometry(self.cfg)
            for page, names in zip(reader.pages, [('N0', 'S0'), ('N1', 'S1')]):
                self.assertAlmostEqual(float(page.mediabox.width), 210/25.4*72, places=3)
                self.assertAlmostEqual(float(page.mediabox.height), 297/25.4*72, places=3)
                for name in names:
                    self.assertIn(name, page.extract_text())
                operations = ContentStream(page['/Contents'], reader).operations
                transforms = [values for values, op in operations if op == b'cm' and float(values[0]) > 10]
                self.assertEqual(len(transforms), 2)
                for values in transforms:
                    self.assertAlmostEqual(float(values[0]), width/25.4*72, places=3)
                    self.assertAlmostEqual(float(values[3]), height/25.4*72, places=3)


if __name__ == '__main__':
    unittest.main()
