from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
from astropy.io import fits
from pypdf import PdfReader
from pypdf.generic import ContentStream

from gingaforge import ohp, starnet_setup
from gingaforge.config import PROJECT_ROOT as ROOT
from gingaforge.ohp import (REGIONS, read_config, geometry, plate_to_sky, relative_solid_angle,
                            celestial_wcs, bilinear, reproject, tone_map, build,
                            channel_mix, gaussian_blur, read_hipparcos, icrs_to_galactic,
                            star_disks, remove_projected_stars, starnet)


class MasterChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = read_config(ROOT / 'mellinger-ohp.properties')
        cls.wcs = celestial_wcs(fits.getheader(ROOT / cls.cfg['source-fits']))

    def test_requested_regions_and_default_dimensions(self):
        self.assertEqual(REGIONS, {'N0': (0, 60), 'N1': (60, 120),
                                   'S0': (300, 360), 'S1': (240, 300)})
        (lens, offset, sphere), _, scale, _, width, height = geometry(self.cfg)
        # focal-length-mm = 19 から、レンズの公式で原板とレンズの距離が決まる。
        self.assertEqual((6500., 1.), (sphere, scale))
        self.assertAlmostEqual(1., 19. * (1 / lens + 1 / (6500. - 50. - lens)), places=12)
        self.assertAlmostEqual(50. + lens, offset)
        self.assertAlmostEqual(lens, 19.0563, places=4)
        # 平らな原板で、四隅（銀経±30度・銀緯±20度）が写る位置で決まる長方形。
        self.assertAlmostEqual(width, 16.2297, places=4)
        self.assertAlmostEqual(height, 22.2954, places=4)

    def test_relative_solid_angle_matches_sky_area(self):
        p, _, scale, _, width, height = geometry(self.cfg)
        h = .01
        def sky_area(x, y):
            # 原盤上の h 四方が写る空の立体角を、銀経・銀緯の差から求める。
            l, b = plate_to_sky(np.array([x - h/2, x + h/2, x, x]), np.array([y, y, y - h/2, y + h/2]), 30., p, scale)
            dl, db = np.radians((l[[2, 3]] - 30 + 180) % 360 - 180), np.radians(b[[0, 1]])
            return abs((dl[1] - dl[0]) * (db[1] - db[0])) * np.cos(np.radians(plate_to_sky(x, y, 30., p, scale)[1]))
        center = sky_area(0., 0.)
        for x, y in ((0., 0.), (width/2, 0.), (0., height/2), (width/2, height/2)):
            with self.subTest(x=x, y=y):
                self.assertAlmostEqual(float(relative_solid_angle(x, y, p, scale)), sky_area(x, y) / center, places=3)
        # 四隅では、中央の約半分の空が原盤の同じ面積に写る。
        self.assertAlmostEqual(.53, float(relative_solid_angle(width/2, height/2, p, scale)), places=2)

    def test_print_scale_shrinks_whole_geometry(self):
        cfg = read_config(ROOT / 'mellinger-ohp.properties')
        _, _, _, _, width, height = geometry(cfg)
        cfg['print-scale'] = '0.5'
        _, _, _, _, half_width, half_height = geometry(cfg)
        self.assertAlmostEqual(width / 2, half_width)
        self.assertAlmostEqual(height / 2, half_height)
        cfg['focal-length-mm'] = '2000'
        with self.assertRaisesRegex(ValueError, '焦点距離'):
            geometry(cfg)

    def test_boundary_orientation_and_source_wcs(self):
        p, bmax, scale, _, width, height = geometry(self.cfg)
        for start, end in REGIONS.values():
            # 画像の四隅が担当領域の四隅。上が銀経の小さい側、右が銀緯のプラス側。
            l, b = plate_to_sky(np.array([-width/2, width/2, -width/2, width/2, 0]),
                                np.array([-height/2, -height/2, height/2, height/2, 0]),
                                (start+end)/2, p, scale)
            np.testing.assert_allclose((l - [end, end, start, start, (start+end)/2] + 180) % 360 - 180, 0, atol=1e-10)
            np.testing.assert_allclose(b, [-bmax, bmax, -bmax, bmax, 0], atol=1e-10)
            # 平らな原板では、辺の中央は担当領域の外になる（reproject() で黒にする）。
            _, b_side = plate_to_sky(np.array([width/2]), np.array([0.]), (start+end)/2, p, scale)
            self.assertGreater(b_side[0], bmax)
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
            for row, col in [(ny//2, nx//2), (ny//3, nx//4), (ny//2, nx-1), (0, nx//2)]:
                x, y = (col+.5)/nx*width-width/2, (row+.5)/ny*height-height/2
                l, b = plate_to_sky(x, y, (start+end)/2, p, scale)
                inside = (l - start) % 360 < end - start and abs(b) <= 20
                sx, sy = self.wcs.world_to_pixel_values(l, b)
                expected = float(bilinear(data, sx, sy) * relative_solid_angle(x, y, p, scale)) if inside else 0.
                self.assertAlmostEqual(float(linear[row, col]), expected, delta=abs(expected) * 1e-6)
            # 銀緯方向の辺の中央は担当領域の外なので黒（0）。
            self.assertEqual(0, linear[ny//2, nx-1])
            pixels = tone_map(linear, cfg)
            self.assertEqual(np.uint8, pixels.dtype)
            self.assertEqual(linear.shape, pixels.shape)

    def test_pdf_pages_labels_and_image_dimensions(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'ohp.properties'
            text = (ROOT / 'mellinger-ohp.properties').read_text()
            text = text.replace('source/mwpan2_RGB_3600.fits', str(ROOT / 'source/mwpan2_RGB_3600.fits'))
            # Layout only: star removal is covered by StarRemovalChecks and needs StarNet2.
            text = text.replace('star-removal = starnet', 'star-removal = none')
            config.write_text(text.replace('raster-dpi = 1200', 'raster-dpi = 72'))
            reader = PdfReader(build(config))
            self.assertEqual(len(reader.pages), 1)
            _, _, _, _, width, height = geometry(self.cfg)
            for page, names in zip(reader.pages, [('N0', 'S0', 'N1', 'S1')]):
                self.assertAlmostEqual(float(page.mediabox.width), 210/25.4*72, places=3)
                self.assertAlmostEqual(float(page.mediabox.height), 297/25.4*72, places=3)
                self.assertIn('B:G:R = 0.5:0.5:0', page.extract_text())
                for name in names:
                    self.assertIn(name, page.extract_text())
                operations = ContentStream(page['/Contents'], reader).operations
                transforms = [values for values, op in operations if op == b'cm' and float(values[0]) > 10]
                self.assertEqual(len(transforms), 4)
                for values in transforms:
                    self.assertAlmostEqual(float(values[0]), width/25.4*72, places=3)
                    self.assertAlmostEqual(float(values[3]), height/25.4*72, places=3)


class StarRemovalChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = read_config(ROOT / 'mellinger-ohp.properties')
        cls.wcs = celestial_wcs(fits.getheader(ROOT / cls.cfg['source-fits']))

    def test_channel_mix_uses_blue_green_red_order(self):
        cube = np.stack([np.full((2, 2), 1.), np.full((2, 2), 3.), np.full((2, 2), 100.)])
        np.testing.assert_allclose(channel_mix(cube, self.cfg), 2.)
        cfg = read_config(ROOT / 'mellinger-ohp.properties')
        cfg['channel-weights'] = '0, 2, 0'
        np.testing.assert_allclose(channel_mix(cube, cfg), 3.)
        cfg['channel-weights'] = '0.5, 0.5'
        with self.assertRaises(ValueError):
            channel_mix(cube, cfg)

    def test_hipparcos_selection_matches_otl(self):
        stars = read_hipparcos(ROOT / self.cfg['star-catalog'], 7.5)
        # OTL projects every Hipparcos star with V <= star.minimum (7.5, inclusive) on one of the two projectors.
        self.assertEqual(len(stars), 25713)
        self.assertEqual(stars[:, 2].max(), 7.5)
        vega = stars[np.argmin(np.hypot(stars[:, 0] - 279.2341, stars[:, 1] - 38.7830))]
        self.assertAlmostEqual(vega[2], 0.03)
        longitude, latitude = icrs_to_galactic(vega[0], vega[1])
        self.assertAlmostEqual(longitude, 67.448, places=2)
        self.assertAlmostEqual(latitude, 19.237, places=2)

    def test_blur_and_disks_wrap_in_longitude(self):
        image = np.zeros((40, 60))
        image[20, 0] = 1.
        blurred = gaussian_blur(image, 2.)
        self.assertAlmostEqual(blurred.sum(), 1.)
        self.assertGreater(blurred[20, 59], 0.)
        mask = star_disks(image.shape, [0.], [20.], [2.])
        self.assertTrue(mask[20, 58] and mask[20, 2] and not mask[20, 3])

    def test_projected_stars_removed_and_fainter_light_kept_as_glow(self):
        vega = [279.2341, 38.7830, 0.03]
        (vx,), (vy,) = self.wcs.world_to_pixel_values(*map(np.atleast_1d, icrs_to_galactic(*vega[:2])))
        vx, vy = round(float(vx)), round(float(vy))
        image = np.full((1800, 3600), 100.)
        image[vy, vx] = 5000.        # projected by the 1st-magnitude projector
        image[300, 1000] = 300.      # a fainter star: excess 200 is below the cap
        image[1500, 2000] = 2100.    # a fainter spike: excess 2000 is capped at faint-star-cap
        starless = np.full(image.shape, 100.)
        starless[vy - 3:vy + 4, vx - 3:vx + 4] += 50.   # halo left by StarNet2 around a bright star
        with mock.patch.object(ohp, 'starnet', return_value=starless), \
                mock.patch.object(ohp, 'read_hipparcos', return_value=np.array([vega])):
            result = remove_projected_stars(image, self.wcs, self.cfg, ROOT)
        self.assertAlmostEqual(result[vy, vx], 100., delta=.5)
        np.testing.assert_allclose(result[vy - 8:vy + 9, vx - 8:vx + 9], 100., atol=.5)
        cap = self.cfg.getfloat('faint-star-cap')
        self.assertAlmostEqual((result - 100.).sum(), 200. + cap, delta=1.)
        self.assertLess(result[300, 1000] - 100., 200. / 10)   # spread over the 0.3 deg glow

    @unittest.skipUnless(starnet_setup.ready(), 'StarNet2 is not installed or its license is not accepted')
    def test_starnet_separates_stars_and_keeps_shape(self):
        image = channel_mix(fits.getdata(ROOT / self.cfg['source-fits'])[:, 644:1156, 1288:2312], self.cfg)
        starless = starnet(image, self.cfg)
        self.assertEqual(starless.shape, image.shape)
        self.assertTrue(np.all(np.isfinite(starless)))
        stars = np.clip(image, 0, ohp.STARNET_RANGE) - starless
        self.assertGreater(np.percentile(stars, 99.9), 1000.)
        self.assertLess(abs(np.median(stars)), 0.1 * np.median(image))


if __name__ == '__main__':
    unittest.main()
