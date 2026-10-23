import json
from pathlib import Path
import shutil
import stat
import tempfile
import unittest
from unittest import mock
import zipfile

from gingaforge import starnet_setup
from gingaforge.gaia_catalog import file_sha256

# A stand-in package: a shell script that answers --version like StarNet2. No StarNet2 files or outputs are used.
SCRIPT = "#!/bin/sh\necho 'starnet2  version: 9.9.9'\n"
LICENSE = "TEST LICENSE AGREEMENT\n"
PLATFORM_KEY = starnet_setup.platform_key  # setUp replaces the module attribute


def make_zip(path, members=()):
    with zipfile.ZipFile(path, "w") as archive:
        executable = zipfile.ZipInfo("starnet2_test_cli/starnet2")
        executable.external_attr = (stat.S_IFREG | 0o755) << 16
        archive.writestr(executable, SCRIPT)
        archive.writestr("starnet2_test_cli/LICENSE.txt", LICENSE)
        archive.writestr("starnet2_test_cli/lib/libfoo.1.dylib", "library")
        link = zipfile.ZipInfo("starnet2_test_cli/lib/libfoo.dylib")
        link.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(link, "libfoo.1.dylib")
        for name, data in members:
            archive.writestr(name, data)


class StarNetSetupChecks(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.root = self.directory / "tools"
        self.archive = self.directory / "starnet2_test_cli.zip"
        make_zip(self.archive)
        self.write_feed(file_sha256(self.archive))
        for patch in (mock.patch.object(starnet_setup, "FEED_URL", (self.directory / "latest.json").as_uri()),
                      mock.patch.object(starnet_setup, "OFFICIAL_PREFIXES", ("file://",)),
                      mock.patch.object(starnet_setup, "platform_key", return_value="macos-arm64")):
            patch.start()
            self.addCleanup(patch.stop)

    def write_feed(self, sha256, url=None):
        archive = {"file_name": self.archive.name, "url": url or self.archive.as_uri(),
                   "sha256": sha256, "size_bytes": self.archive.stat().st_size}
        feed = {"schema": starnet_setup.FEED_SCHEMA, "tools": {"starnet2": {"latest": {"macos-arm64": {
            "version": "9.9.9", "build": "0001", "packages": {"archive": archive}}}}}}
        (self.directory / "latest.json").write_text(json.dumps(feed))

    def test_platform_keys_follow_the_official_feed(self):
        for system, machine, key in [("Darwin", "arm64", "macos-arm64"), ("Darwin", "x86_64", "macos-x64"),
                                     ("Linux", "x86_64", "linux-x64"), ("Windows", "AMD64", "windows-x64")]:
            with mock.patch("platform.system", return_value=system), mock.patch("platform.machine", return_value=machine):
                self.assertEqual(key, PLATFORM_KEY())
        with mock.patch("platform.system", return_value="Linux"), mock.patch("platform.machine", return_value="aarch64"):
            with self.assertRaises(RuntimeError):
                PLATFORM_KEY()

    def test_download_is_verified_and_installed_after_acceptance(self):
        executable = starnet_setup.install(accept=True, root=self.root)
        self.assertEqual(self.root / "starnet2_test_cli" / "starnet2", executable)
        self.assertTrue(executable.stat().st_mode & stat.S_IXUSR)
        self.assertTrue((executable.parent / "lib/libfoo.dylib").is_symlink())
        self.assertEqual("9.9.9", starnet_setup.version(executable))
        self.assertEqual(executable, starnet_setup.installed_executable(self.root))
        # The copy under tools/starnet2/ is preferred over a starnet2 on PATH.
        with mock.patch("shutil.which", return_value="/usr/local/bin/starnet2"):
            self.assertEqual(executable, starnet_setup.locate("", self.root))
        self.assertTrue(starnet_setup.ready("", self.root))

    def test_checksum_mismatch_installs_nothing(self):
        self.write_feed("0" * 64)
        with self.assertRaises(ValueError):
            starnet_setup.install(accept=True, root=self.root)
        self.assertIsNone(starnet_setup.installed_executable(self.root))
        self.assertEqual([], list(self.root.iterdir()))

    def test_unofficial_download_url_is_rejected(self):
        self.write_feed(file_sha256(self.archive), url="https://example.com/starnet2.zip")
        with mock.patch.object(starnet_setup, "OFFICIAL_PREFIXES", ("https://download.starnetastro.com/",)):
            with self.assertRaises(ValueError):
                starnet_setup.latest_archive()

    def test_license_is_accepted_once_per_text(self):
        with self.assertRaises(PermissionError):
            starnet_setup.install(interactive=False, root=self.root)
        self.assertIsNone(starnet_setup.installed_executable(self.root))
        with mock.patch("builtins.input", return_value="y"), mock.patch("builtins.print"):
            executable = starnet_setup.install(interactive=True, root=self.root)
        starnet_setup.ensure_license(executable, interactive=False, root=self.root)
        (executable.parent / "LICENSE.txt").write_text(LICENSE + "Revised.\n")
        self.assertFalse(starnet_setup.ready("", self.root))
        with self.assertRaises(PermissionError):
            starnet_setup.ensure_license(executable, interactive=False, root=self.root)
        with mock.patch("builtins.input", return_value="n"), mock.patch("builtins.print"):
            with self.assertRaises(PermissionError):
                starnet_setup.ensure_license(executable, interactive=True, root=self.root)
        starnet_setup.ensure_license(executable, accept=True, root=self.root)
        self.assertTrue(starnet_setup.ready("", self.root))

    def test_missing_starnet_offers_installation_only_when_interactive(self):
        with mock.patch("shutil.which", return_value=None):
            with self.assertRaisesRegex(FileNotFoundError, "install-starnet"):
                starnet_setup.prepare("", interactive=False, root=self.root)
            with mock.patch("builtins.input", side_effect=["y", "y"]), mock.patch("builtins.print"):
                executable = starnet_setup.prepare("", interactive=True, root=self.root)
        self.assertEqual(executable, starnet_setup.installed_executable(self.root))

    def test_install_from_a_manually_downloaded_folder(self):
        unpacked = self.directory / "unpacked"
        starnet_setup.extract(self.archive, unpacked)
        executable = starnet_setup.install(unpacked / "starnet2_test_cli", accept=True, root=self.root)
        self.assertEqual("9.9.9", starnet_setup.version(executable))
        self.assertTrue((unpacked / "starnet2_test_cli" / "starnet2").is_file())

    def test_zip_entries_outside_the_target_are_rejected(self):
        for name, data, link in [("../evil.txt", "x", False), ("starnet2_test_cli/escape", "../../evil", True)]:
            bad = self.directory / "bad.zip"
            with zipfile.ZipFile(bad, "w") as archive:
                info = zipfile.ZipInfo(name)
                if link:
                    info.external_attr = (stat.S_IFLNK | 0o777) << 16
                archive.writestr(info, data)
            with self.assertRaises(ValueError):
                starnet_setup.extract(bad, self.directory / "target")
            self.assertFalse((self.directory / "evil.txt").exists())
        shutil.rmtree(self.directory / "target", ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
