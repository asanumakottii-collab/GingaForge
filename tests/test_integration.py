import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from gingaforge import starnet_setup
from gingaforge.galaxy_geometry import GalaxyTransformer, REGIONS
from gingaforge.models import SpherePosition
from tests.test_gaia_etching import write_cache

ROOT = Path(__file__).resolve().parents[1]


class IntegrationTests(unittest.TestCase):
    def run_command(self, *args, cwd, timeout=60):
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
        return subprocess.run([sys.executable, "-m", "gingaforge", *args],
                              cwd=cwd, env=env, text=True, capture_output=True, timeout=timeout)

    def test_shared_projection_matches_ohp_in_all_four_regions(self):
        from gingaforge.ohp import plate_to_sky
        for distance, scale in ((50., 1.), (70., 2.)):
            transformer = GalaxyTransformer(.1, 6500. * scale, 0., 0., distance * scale)
            for region, (start, end) in enumerate(REGIONS.values()):
                center = (start + end) / 2
                for longitude in (start, center, end):
                    for latitude in (-20., 0., 20.):
                        # Explicit unit retains the upper boundary for orientation checks.
                        p = transformer.transform_unit(SpherePosition.from_galactic(longitude, latitude),
                                                       *divmod(region, 2))
                        actual_l, actual_b = plate_to_sky(p.xmm, p.ymm, center, distance, scale)
                        self.assertAlmostEqual(0., (actual_l - longitude + 180) % 360 - 180, places=9)
                        self.assertAlmostEqual(latitude, actual_b, places=9)

    def test_gaia_config_relative_io(self):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configs = root / 'configs'
            write_cache(configs / 'cache', {0: [(30.05, .05, .01, 1)]})
            config = configs / 'custom.properties'
            config.write_text('galaxy.gaia.input=cache\nplate.frame-size=60\ncolor.invert=y\n')
            result = self.run_command('gaia-etch', '--pdf', '-f', str(config), cwd=root)
            self.assertEqual(0, result.returncode, result.stderr)
            output = configs / 'output/custom/inverted'
            self.assertFalse((root / 'output').exists())
            self.assertEqual(['custom-0.pdf', 'custom-1.pdf'], sorted(p.name for p in output.iterdir()))
            self.assertEqual(['N0', 'S0'], PdfReader(output / 'custom-0.pdf').pages[0].extract_text().split())

    def test_ohp_command_generates_from_config_relative_paths(self):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configs = root / 'configs'
            configs.mkdir()
            (configs / 'sky.fits').symlink_to(ROOT / 'source/mwpan2_RGB_3600.fits')
            (configs / 'stars.dat.gz').symlink_to(ROOT / 'source/hip_main.dat.gz')
            config = configs / 'ohp.properties'
            text = (ROOT / 'mellinger-ohp.properties').read_text()
            text = text.replace('source/mwpan2_RGB_3600.fits', 'sky.fits').replace('raster-dpi = 1200', 'raster-dpi = 72')
            text = text.replace('source/hip_main.dat.gz', 'stars.dat.gz')
            if starnet_setup.ready():
                # Without 2x upsampling StarNet2 takes about half as long; quality is not checked here.
                text = text.replace('starnet-upsample = yes', 'starnet-upsample = no')
            else:
                text = text.replace('star-removal = starnet', 'star-removal = none')
            config.write_text(text)
            result = self.run_command('mellinger-ohp', '--config', str(config), cwd=root, timeout=300)
            self.assertEqual(0, result.returncode, result.stderr)
            output = configs / 'output/ohp'
            self.assertFalse((root / 'output').exists())
            self.assertEqual(['ohp-plates.pdf'], [p.name for p in output.iterdir()])
            self.assertEqual(2, len(PdfReader(output / 'ohp-plates.pdf').pages))

    def test_cli_help_and_invalid_options(self):
        with tempfile.TemporaryDirectory() as directory:
            for args in (('--help',), ('mellinger-ohp', '--help'), ('gaia-etch', '--help'), ('gaia-ohp', '--help'), ('install-starnet', '--help')):
                result = self.run_command(*args, cwd=directory)
                self.assertEqual(0, result.returncode, result.stderr)
                self.assertIn('usage:', result.stdout)
            for args in (('legacy',), ('ohp',), ('gaia-etch', '--unknown'), ('gaia-etch', '-f'), ('gaia-ohp', '--unknown'), ('gaia-etch', '-svg', '-pdf')):
                self.assertEqual(2, self.run_command(*args, cwd=directory).returncode)


if __name__ == '__main__':
    unittest.main()
