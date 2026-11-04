import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from tests.test_gaia_etching import write_cache

ROOT = Path(__file__).resolve().parents[1]


class IntegrationTests(unittest.TestCase):
    def run_command(self, *args, cwd, timeout=60):
        env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
        return subprocess.run([sys.executable, "-m", "gingaforge", *args],
                              cwd=cwd, env=env, text=True, capture_output=True, timeout=timeout)

    def test_gaia_config_relative_io(self):
        from pypdf import PdfReader
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            configs = root / 'configs'
            write_cache(configs / 'cache', {0: [(30.05, .05, .01, 1)]})
            config = configs / 'custom.properties'
            config.write_text('galaxy.gaia.input=cache\nplate.frame-diameter=60\ncolor.invert=y\n')
            result = self.run_command('gaia-etch', '--pdf', '-f', str(config), cwd=root)
            self.assertEqual(0, result.returncode, result.stderr)
            output = configs / 'output/custom/inverted'
            self.assertFalse((root / 'output').exists())
            self.assertEqual(['custom-0.pdf', 'custom-1.pdf'], sorted(p.name for p in output.iterdir()))
            self.assertEqual(['N0', 'S0'], PdfReader(output / 'custom-0.pdf').pages[0].extract_text().split())

    def test_cli_help_and_invalid_options(self):
        with tempfile.TemporaryDirectory() as directory:
            for args in (('--help',), ('gaia-etch', '--help'), ('gaia-ohp', '--help')):
                result = self.run_command(*args, cwd=directory)
                self.assertEqual(0, result.returncode, result.stderr)
                self.assertIn('usage:', result.stdout)
            for args in (('legacy',), ('ohp',), ('mellinger-ohp',), ('install-starnet',), ('gaia-etch', '--unknown'), ('gaia-etch', '-f'), ('gaia-ohp', '--unknown'), ('gaia-etch', '-svg', '-pdf')):
                self.assertEqual(2, self.run_command(*args, cwd=directory).returncode)


if __name__ == '__main__':
    unittest.main()
