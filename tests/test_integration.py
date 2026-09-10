import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from io import StringIO

from galaxy_geometry import GalaxyTransformer, REGIONS
from galaxy_transformer import main as gaia_main
from models import SpherePosition

ROOT = Path(__file__).resolve().parents[1]


class IntegrationTests(unittest.TestCase):
    def run_command(self, *args, cwd):
        return subprocess.run([sys.executable, str(ROOT / 'generate.py'), *args],
                              cwd=cwd, text=True, capture_output=True, timeout=30)

    def test_shared_projection_matches_ohp_in_all_four_regions(self):
        from build_mellinger_ohp import plate_to_sky
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

    def test_default_query_is_standalone_and_has_no_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_command('gaia-etch', '--gaia-query', cwd=directory)
            self.assertEqual(0, result.returncode, result.stderr)
            self.assertEqual(4, result.stdout.count('FROM gaiadr3.gaia_source'))
            self.assertEqual([], list(Path(directory).iterdir()))

    def test_gaia_config_relative_io_and_no_mode_setting(self):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configs = root / 'configs'
            configs.mkdir()
            point = SpherePosition.from_galactic(30.01, .01)
            (configs / 'stars.csv').write_text(
                f'source_id,ra,dec,phot_g_mean_mag\n1,{point.radeg % 360},{point.dedeg},10\n')
            config = configs / 'custom.properties'
            config.write_text('galaxy.gaia.input=stars.csv\noutput.directory=masters\n'
                              'galaxy.gaia.file.prefix=custom-\nplate.frame-size=60\ncolor.invert=y\n')
            result = self.run_command('gaia-etch', '--pdf', '-f', str(config), cwd=root)
            self.assertEqual(0, result.returncode, result.stderr)
            output = configs / 'masters/galaxy'
            self.assertFalse((root / 'masters').exists())
            self.assertEqual(2, len(list(output.glob('*.pdf'))))
            self.assertEqual(['N0', 'S0'], PdfReader(output / 'custom-0.pdf').pages[0].extract_text().split())
            report = json.loads((output / 'custom-report.json').read_text())
            self.assertEqual(60., report['frame_size_output_mm'])
            self.assertEqual(1, report['input']['stars'])

    def test_legacy_rejected_even_with_download_or_query_flag(self):
        with tempfile.TemporaryDirectory() as directory:
            config = Path(directory) / 'legacy.properties'
            config.write_text('galaxy.mode=legacy\n')
            for flags in ([], ['--download-gaia'], ['--gaia-query']):
                with (patch('gaia_catalog.download_gaia', side_effect=AssertionError('unexpected network')),
                      patch('builtins.input', side_effect=AssertionError('unexpected prompt')),
                      patch('sys.stdout', new_callable=StringIO),
                      patch('sys.stderr', new_callable=StringIO)):
                    self.assertEqual(2, gaia_main(['-f', str(config), *flags]))
            self.assertEqual([config], list(Path(directory).iterdir()))

    def test_ohp_command_generates_from_config_relative_paths(self):
        from PIL import Image
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configs = root / 'configs'
            configs.mkdir()
            (configs / 'sky.fits').symlink_to(ROOT / 'source/mwpan2_RGB_3600.fits')
            config = configs / 'ohp.properties'
            text = (ROOT / 'galaxyconfig.properties').read_text()
            text = text.replace('source/mwpan2_RGB_3600.fits', 'sky.fits').replace('raster-dpi = 1200', 'raster-dpi = 72')
            config.write_text(text)
            result = self.run_command('ohp', '--config', str(config), cwd=root)
            self.assertEqual(0, result.returncode, result.stderr)
            output = configs / 'output'
            self.assertFalse((root / 'output').exists())
            self.assertEqual(2, len(PdfReader(output / 'mellinger-ohp-plates.pdf').pages))
            self.assertEqual(1, len(PdfReader(output / 'ohp-calibration.pdf').pages))
            for name in REGIONS:
                with Image.open(output / f'{name}-ohp.png') as image:
                    self.assertEqual((103, 148), image.size)

    def test_cli_help_and_invalid_options(self):
        with tempfile.TemporaryDirectory() as directory:
            for args in (('--help',), ('ohp', '--help'), ('gaia-etch', '--help')):
                result = self.run_command(*args, cwd=directory)
                self.assertEqual(0, result.returncode, result.stderr)
                self.assertIn('usage:', result.stdout)
            for args in (('legacy',), ('gaia-etch', '--unknown'), ('gaia-etch', '-f'), ('gaia-etch', '-PS')):
                self.assertEqual(2, self.run_command(*args, cwd=directory).returncode)

if __name__ == '__main__':
    unittest.main()
